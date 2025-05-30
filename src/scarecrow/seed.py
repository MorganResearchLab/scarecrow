# -*- coding: utf-8 -*-
"""
#!/usr/bin/env python3
@author: David Wragg
"""

import os
import pysam
import logging
import random
from argparse import ArgumentTypeError, RawTextHelpFormatter
from collections import defaultdict, namedtuple
from functools import lru_cache
from typing import List, Dict, Set, Tuple
from scarecrow import __version__
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import (
    generate_random_string,
    reverse_complement,
    parse_seed_arguments,
)
from scarecrow.encode import BarcodeMatcherAhoCorasick

# Define a namedtuple for match results
MatchResult = namedtuple(
    "MatchResult",
    ["barcode", "whitelist", "orientation", "start", "end", "mismatches", "distance"],
)


class BarcodeMatcher:
    def __init__(self, barcode_sequences: Dict[str, Set[str]], mismatches: int = 0):
        """
        Base class for barcode matching strategies.
        """
        self.mismatches = mismatches
        self.barcode_sequences = barcode_sequences

    def find_matches(self, sequence: str) -> List[MatchResult]:
        """
        Find all matching barcodes in a sequence.
        Implemented by subclasses.
        """
        raise NotImplementedError


class BarcodeMatcherOptimized(BarcodeMatcher):
    def __init__(self, barcode_sequences: Dict[str, Set[str]], mismatches: int = 0):
        """
        Initialize optimized barcode matcher with pre-computed lookup tables.
        """
        super().__init__(barcode_sequences, mismatches)
        self.matchers = {}
        self.logger = setup_worker_logger()
        self.logger.info("Initializing set based barcode matcher")

        for whitelist_key, sequences in barcode_sequences.items():
            exact_matches = set(sequences)
            mismatch_lookup = (
                self._create_mismatch_lookup(sequences) if mismatches > 0 else None
            )
            self.matchers[whitelist_key] = {
                "exact": exact_matches,
                "mismatch": mismatch_lookup,
                "length": len(next(iter(sequences))),
            }

    def _create_mismatch_lookup(self, sequences: Set[str]) -> Dict[str, str]:
        """Create lookup table for sequences with allowed mismatches."""
        lookup = {}
        for seq in sequences:
            lookup[seq] = (seq, 0)  # (original sequence, mismatch count)
            if self.mismatches > 0:
                for i in range(len(seq)):
                    for base in "ACGTN":
                        if base != seq[i]:
                            variant = seq[:i] + base + seq[i + 1 :]
                            if variant not in lookup:
                                lookup[variant] = (seq, 1)
        return lookup

    @lru_cache(maxsize=1024)
    def _reverse_complement(self, sequence: str) -> str:
        """Cached reverse complement computation."""
        return reverse_complement(sequence)

    def find_matches(
        self,
        sequence: List[tuple[str, int]],
        whitelist_key: str,
        orientation: str,
        original_start: int,
    ) -> List[Dict]:
        """
        Find all matching barcodes in a sequence.
        Returns list of matches with positions and orientations.
        """
        matches = []
        sequence, start = sequence[0]
        seq_len = len(sequence)

        for whitelist_key, matcher in self.matchers.items():
            barcode_len = matcher["length"]

            # Scan sequence for matches in both orientations
            for start in range(seq_len - barcode_len + 1):
                candidate = sequence[start : start + barcode_len]

                # Reverse sequence if required
                if orientation == "reverse":
                    candidate = self._reverse_complement(candidate)

                # Check for matches
                match = self._check_match(candidate, matcher)
                if match:
                    matches.append(
                        MatchResult(
                            barcode = match[0],
                            whitelist = whitelist_key,
                            orientation = orientation,
                            start = start + 1,  # 1-based position
                            end = start + barcode_len,
                            mismatches = match[1],
                            distance = 0,
                        )
                    )

        return sorted(matches, key=lambda x: x.start)

    def _check_match(self, candidate: str, matcher: Dict) -> Tuple[str, int]:
        """Check if candidate matches any barcode in the matcher."""
        if candidate in matcher["exact"]:
            return (candidate, 0)
        if (
            self.mismatches > 0
            and matcher["mismatch"]
            and candidate in matcher["mismatch"]
        ):
            return matcher["mismatch"][candidate]
        return None


class SequenceFrequencyAnalyzer:
    def __init__(self):
        """Initialize analyzer to track nucleotide frequencies at each position"""
        self.position_counts = defaultdict(lambda: defaultdict(int))
        self.total_reads = 0

    def add_sequence(self, sequence: str) -> None:
        """
        Add a sequence to the frequency analysis

        Args:
            sequence (str): DNA sequence to analyze
        """
        self.total_reads += 1
        for pos, base in enumerate(sequence):
            self.position_counts[pos][base] += 1

    def get_frequencies(self) -> Dict[int, Dict[str, float]]:
        """
        Calculate frequencies for each nucleotide at each position

        Returns:
            Dict[int, Dict[str, float]]: Position -> {base -> frequency}
        """
        frequencies = {}
        for pos in sorted(self.position_counts.keys()):
            frequencies[pos] = {}
            for base in "ACGTN":
                count = self.position_counts[pos].get(base, 0)
                frequencies[pos][base] = count / self.total_reads
        return frequencies

    def find_conserved_runs(self, min_freq: float, min_length: int) -> List[Dict]:
        """
        Find runs of conserved bases meeting frequency and length thresholds

        Args:
            min_freq (float): Minimum frequency threshold (0-1)
            min_length (int): Minimum length of conserved run

        Returns:
            List[Dict]: List of conserved runs with start, end, sequence, and frequencies
        """
        frequencies = self.get_frequencies()
        conserved_runs = []
        current_run = None

        for pos in sorted(frequencies.keys()):
            # Find most frequent base at this position
            pos_freqs = frequencies[pos]
            max_base = max(pos_freqs.items(), key=lambda x: x[1])
            base, freq = max_base

            # Check if this position is conserved
            if freq >= min_freq:
                if current_run is None:
                    # Start new run
                    current_run = {"start": pos, "bases": [base], "freqs": [freq]}
                else:
                    # Extend current run
                    current_run["bases"].append(base)
                    current_run["freqs"].append(freq)
            else:
                # Check if we should save the current run
                if current_run and len(current_run["bases"]) >= min_length:
                    conserved_runs.append(self._finalize_run(current_run))
                current_run = None

        # Check final run
        if current_run and len(current_run["bases"]) >= min_length:
            conserved_runs.append(self._finalize_run(current_run))

        return conserved_runs

    def _finalize_run(self, run: Dict) -> Dict:
        """Helper to finalize a conserved run with summary statistics"""
        return {
            "start": run["start"] + 1,  # Convert to 1-based position
            "end": run["start"] + len(run["bases"]),
            "length": len(run["bases"]),
            "sequence": "".join(run["bases"]),
            "median_freq": sorted(run["freqs"])[len(run["freqs"]) // 2],
        }

    def write_frequencies(self, output_file: str, read_label: str) -> None:
        """
        Write frequencies to output file

        Args:
            output_file (str): Path to output file
            read_label (str): Label for the read (e.g., 'read1' or 'read2')
        """
        frequencies = self.get_frequencies()
        with open(output_file, "w") as f:
            # Write header
            f.write("read\tposition\tA\tC\tG\tT\tN\n")

            # Write frequencies for each position
            for pos in sorted(frequencies.keys()):
                freq = frequencies[pos]
                f.write(f"{read_label}\t{pos + 1}\t")  # 1-based position
                f.write(f"{freq.get('A', 0):.4f}\t")
                f.write(f"{freq.get('C', 0):.4f}\t")
                f.write(f"{freq.get('G', 0):.4f}\t")
                f.write(f"{freq.get('T', 0):.4f}\t")
                f.write(f"{freq.get('N', 0):.4f}\n")


def parser_seed(parser):
    subparser = parser.add_parser(
        "seed",
        description="""
Search reads in FASTQ files for barcodes in whitelists, record exact matches and conserved sequences.

Example:

scarecrow seed --fastqs R1.fastq.gz R2.fastq.gz\n\t--barcodes BC1:BC1.txt BC2:BC2.txt BC3:BC3.txt\n\t--out barcode_counts.csv 
---
""",
        help="Search fastq reads for barcodes",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument(
        "-f", "--fastqs", nargs="+", required=True, help="Pair of FASTQ files"
    )
    subparser.add_argument(
        "-n",
        "--num_reads",
        metavar = "<int>",
        help = ("Number of read pairs to sample [100000]"),
        type = int,
        default = 100000,
    )
    subparser.add_argument(
        "-u",
        "--upper_read_count",
        metavar = "<int>",
        help = ("Upper read count of reads to subsample [10000000]"),
        type = int,
        default = 10000000,
    )
    subparser.add_argument(
        "-r",
        "--random_seed",
        metavar = "<int>",
        help = ("Random seed for sampling read pairs [1234]"),
        type = int,
        default = 1234,
    )
    subparser.add_argument(
        "-c",
        "--barcodes",
        metavar = "<string>",
        nargs = "+",
        required = True,
        help = "Barcode whitelist files in format <barcode_name>:<whitelist_name>:<whitelist_file>\n\t(e.g. BC1:v1:barcodes1.txt BC2:v2:barcodes2.txt ...)",
    )
    subparser.add_argument(
        "-o",
        "--out",
        metavar = "<file>",
        help = (
            "CSV file to write barcode counts to, also serves as prefix for conserved.tsv and frequencies.tsv files."
        ),
        type = str,
        default = "./barcode_counts.csv",
    )
    subparser.add_argument(
        "-p",
        "--pickle",
        metavar = "<file>",
        help = ("Path to compressed pickle file to use for barcode matching [None]"),
        type = str,
        default = None,
    )
    subparser.add_argument(
        "-k",
        "--kmer_length",
        metavar = "<int>",
        help = ("K-mer length for building k-mer index for approximate matching"),
        type = int,
        default = None,
    )
    subparser.add_argument(
        "-b",
        "--batch_size",
        metavar = "<int>",
        help = ("Number of read pairs per batch to process at a time [10000]"),
        type = int,
        default = 10000,
    )
    subparser.add_argument(
        "-l",
        "--linker_base_frequency",
        metavar = "<float>",
        help = ("Nucleotide frequency indicative of linker (fixed) sequence [0.75]"),
        type = float,
        default = 0.75,
    )
    subparser.add_argument(
        "-m",
        "--linker_min_length",
        metavar = "<int>",
        help = (
            "Minimum run of bases exceeding linker_base_frequency to suggest linker sequence [10]"
        ),
        type = int,
        default = 10,
    )
    subparser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose output [false]"
    )
    return subparser


def validate_seed_args(parser, args):
    """
    Validate arguments
    """

    # Type checking for all arguments
    try:
        # Check fastqs - should be list of strings
        if not isinstance(args.fastqs, list) or not all(isinstance(f, str) for f in args.fastqs):
            raise TypeError("--fastqs must be a list of file paths")
        
        # Check that all FASTQ files exist and are readable
        for fastq in args.fastqs:
            if not os.path.exists(fastq):
                raise ArgumentTypeError(f"FASTQ file does not exist: {fastq}")
            if not os.access(fastq, os.R_OK):
                raise ArgumentTypeError(f"FASTQ file is not readable: {fastq}")

        # Check numeric parameters
        if not isinstance(args.num_reads, int):
            raise TypeError("--num_reads must be a positive integer")
        if not isinstance(args.upper_read_count, int):
            raise TypeError("--upper_read_count must be a positive integer")
        if not isinstance(args.random_seed, int):
            raise TypeError("--random_seed must be an integer")
        if not isinstance(args.batch_size, int) or args.batch_size <= 0:
            raise TypeError("--batch_size must be a positive integer")
        if not isinstance(args.linker_min_length, int) or args.linker_min_length <= 0:
            raise TypeError("--linker_min_length must be a positive integer")
        
        # Check float parameters
        if not isinstance(args.linker_base_frequency, float) or not (0 <= args.linker_base_frequency <= 1):
            raise TypeError("--linker_base_frequency must be a float between 0 and 1")
        
        # Check barcodes format
        if not isinstance(args.barcodes, list) or not all(isinstance(bc, str) for bc in args.barcodes):
            raise TypeError("--barcodes must be a list of strings in format 'name:version:file'")
        
        # Check output file is string
        if not isinstance(args.out, str):
            raise TypeError("--out must be a string file path")

        # Check if output directory exists
        output_dir = os.path.dirname(os.path.abspath(args.out)) or '.'
        if not os.path.exists(output_dir):
            raise ArgumentTypeError(f"Output directory does not exist: {output_dir}")
        if not os.access(output_dir, os.W_OK):
            raise ArgumentTypeError(f"No write permissions for output directory: {output_dir}")

        # Check optional parameters
        if args.pickle is not None and not isinstance(args.pickle, str):
            raise TypeError("--pickle must be a string file path or None")
        if args.kmer_length is not None and (not isinstance(args.kmer_length, int) or args.kmer_length <= 0):
            raise TypeError("--kmer_length must be a positive integer or None")
        
            
    except TypeError as e:
        parser.error(str(e))

    # Global logger setup
    logfile = f"./scarecrow_seed_{generate_random_string()}.log"
    logger = setup_logger(logfile)
    logger.info(f"scarecrow version {__version__}")
    logger.info(f"logfile: '{logfile}'")

    run_seed(
        fastqs = [f for f in args.fastqs],
        num_reads = args.num_reads,
        random_seed = args.random_seed,
        upper_read_count = args.upper_read_count,
        barcodes = args.barcodes,
        pickle_file = args.pickle,
        kmer_length = args.kmer_length,
        output_file = args.out,
        batches = args.batch_size,
        linker_base_frequency = args.linker_base_frequency,
        linker_min_length = args.linker_min_length,
        verbose = args.verbose,
    )


@log_errors
def process_read_batch(
    read_batch: List[Tuple[int, str, pysam.FastxRecord]],
    verbose: bool = False,
    matcher: BarcodeMatcherOptimized = None,
    whitelist_key: str = None,
) -> List[Dict]:
    """
    Process a batch of reads efficiently
    """
    logger = logging.getLogger("scarecrow")
    results = []

    for file_index, file_name, read in read_batch:
        for orientation in ["forward", "reverse"]:
            read_info = {
                "file_index": file_index,
                "file_name": file_name,
                "header": read.name,
                "regions": matcher.find_matches(
                    [(read.sequence, 1)], whitelist_key, orientation, 1
                ),
                "seqlen": len(read.sequence),
            }            
            results.append(read_info)
            if verbose:
                logger.info(
                    f"File index:{file_index}\t{file_name}\tRead\t{read.name}\n{read.sequence}\n{read_info['regions']}"
                )

    # Filter data to ensure no more than one barcode recorded for a given read at a given position and orientation
    filtered_data = filter_data(results)

    return filtered_data


# Function to filter data for unique combinations
def filter_data(data):
    filtered_data = []

    for entry in data:
        file_index = entry["file_index"]
        file_name = entry["file_name"]
        header = entry["header"]
        regions = entry["regions"]
        seqlen = entry["seqlen"]

        seen = set()
        filtered_regions = []

        for region in regions:
            # Create a key with the unique combination
            key = (file_index, file_name, header, region.whitelist, region.orientation, region.start)

            # If this combination hasn't been seen before, add it to the result
            if key not in seen:
                seen.add(key)
                filtered_regions.append(region)

        # Store the filtered regions back into the same structure
        filtered_data.append(
            {
                "file_index": file_index,
                "file_name": file_name,
                "header": header,
                "regions": filtered_regions,
                "seqlen": seqlen,
            }
        )

    return filtered_data


@log_errors
def run_seed(
    fastqs: List[str] = None,
    num_reads: int = 0,
    random_seed: int = 1234,
    upper_read_count: int = 10000000,
    barcodes: List[str] = None,
    pickle_file: str = None,
    kmer_length: int = None,
    output_file: str = None,
    batches: int = 10000,
    linker_base_frequency: float = 0.75,
    linker_min_length: int = 10,
    verbose: bool = False,
) -> None:
    """
    Optimized implementation of seed functionality
    """
    logger = logging.getLogger("scarecrow")
    random.seed(random_seed)

    # Parse barcode whitelists
    expected_barcodes = parse_seed_arguments(barcodes)
    whitelist_key = list(expected_barcodes.keys())[0]

    # Initialize matcher based on the chosen method
    if pickle_file:
        matcher = BarcodeMatcherAhoCorasick(
            barcode_sequences = {k: set(v) for k, v in expected_barcodes.items()},
            pickle_file = pickle_file,
            kmer_length = kmer_length,
            mismatches = 0,
        )
    else:
        matcher = BarcodeMatcherOptimized(
            barcode_sequences = {k: set(v) for k, v in expected_barcodes.items()},
            mismatches = 0,
        )

    # Initialize sequence analyzers for each read
    analyzers = [SequenceFrequencyAnalyzer() for _ in fastqs]

    # If subsetting FASTQ, first get total read count
    if num_reads > 0:
        if upper_read_count == 0:
            logger.info(
                "For subsetting, getting total read count to generate random indices for sampling."
            )
            upper_read_count = sum(1 for _ in pysam.FastxFile(fastqs[0]))
        num_reads = min(num_reads, upper_read_count)
        sample_indices = set(random.sample(range(upper_read_count), num_reads))
        logger.info(
            f"Selected {len(sample_indices)} random reads out of an upper limit of {upper_read_count} reads"
        )
        logger.info(f"Random seed used: {random_seed}")

    # Process files with minimal overhead
    fastq_handles = [pysam.FastxFile(fastq) for fastq in fastqs]
    with open(output_file, "w") as out_csv:
        # Write header
        out_csv.write(
            "file_index\tfile\tread_name\tseqlen\tbarcode_whitelist\tbarcode\torientation\tstart\tend\tmismatches\n"
        )

        # Create batches efficiently
        current_batch = []
        for idx, reads in enumerate(zip(*fastq_handles)):
            if num_reads > 0:
                # Current read index exceeds max index in subsample so exit
                if idx > max(sample_indices):
                    break
                # Current read index not in subsample index so skip to next read
                if idx not in sample_indices:
                    continue

            # Add sequences to analyzers and track file index
            for file_index, read in enumerate(reads):
                analyzers[file_index].add_sequence(read.sequence)
                current_batch.append((file_index, os.path.basename(fastqs[file_index]), read)) 

            if len(current_batch) >= batches:
                # Process batch
                results = process_read_batch(current_batch, verbose, matcher, whitelist_key)

                # Write results
                for result in results:
                    write_batch_results(result, out_csv)

                current_batch = []

        # Process remaining reads
        if current_batch:
            results = process_read_batch(current_batch, verbose, matcher, whitelist_key)
            for result in results:
                write_batch_results(result, out_csv)

    logger.info(f"Barcode seeds written to\n\t'{output_file}'")

    # Write frequency analysis results
    freq_output_base = os.path.splitext(output_file)[0]
    for file_index, analyzer in enumerate(analyzers):
        analyzer.write_frequencies(
            f"{freq_output_base}_file{file_index}_frequencies.tsv", f"file{file_index}"
        )
    logger.info(
        f"Sequence base frequency results written to:\n\t'{freq_output_base}_file*_frequencies.tsv'"
    )

    # Find and write conserved runs
    with open(f"{freq_output_base}_conserved.tsv", "w") as f:
        f.write("file_index\tstart\tend\tlength\tsequence\tmedian_frequency\n")

        # Analyze each file
        for file_index, analyzer in enumerate(analyzers):
            runs = analyzer.find_conserved_runs(
                linker_base_frequency, linker_min_length
            )
            for run in runs:
                f.write(
                    f"{file_index}\t{run['start']}\t{run['end']}\t{run['length']}\t{run['sequence']}\t{run['median_freq']:.4f}\n"
                )
                logger.info(
                    f"Possible conserved sequence in file {file_index}:\n\t'{run['start']}-{run['end']}'\t'{run['sequence']}'"
                )

    logger.info(
        f"Conserved sequence analysis results written to:\n\t'{freq_output_base}_conserved.tsv'"
    )


    logger.info("Finished!")


def write_batch_results(read_info: Dict, output_handler) -> None:
    """
    Write batch results to output file
    """
    file_index = read_info["file_index"]
    file_name = read_info["file_name"]
    header = read_info["header"]
    seqlen = read_info["seqlen"]

    for region in read_info["regions"]:
        output_handler.write(
            f"{file_index}\t{file_name}\t{header}\t{seqlen}\t{region.whitelist}\t{region.barcode}\t"
            f"{region.orientation}\t{region.start}\t{region.end}\t{region.mismatches}\n"
        )

def setup_worker_logger(log_file: str = None):
    """Configure logger for worker processes with file output"""
    logger = logging.getLogger("scarecrow")
    if not logger.handlers:  # Only add handlers if none exist
        # Create formatters
        formatter = logging.Formatter(
            "%(asctime)s - %(processName)s - %(name)s - %(levelname)s - %(message)s"
        )

        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

        # File handler if log_file provided
        if log_file:
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

        logger.setLevel(logging.INFO)
    return logger
