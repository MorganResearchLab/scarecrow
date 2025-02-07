#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
"""

import os
import pysam
import logging
from argparse import RawTextHelpFormatter
from collections import defaultdict
from functools import lru_cache
from typing import List, Dict, Set, Tuple
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string, reverse_complement

class BarcodeMatcherOptimized:
    def __init__(self, barcode_sequences: Dict[str, Set[str]], mismatches: int = 1):
        """
        Initialize optimized barcode matcher with pre-computed lookup tables
        """
        self.mismatches = mismatches
        self.matchers = {}
        
        for whitelist_key, sequences in barcode_sequences.items():
            exact_matches = set(sequences)
            mismatch_lookup = self._create_mismatch_lookup(sequences) if mismatches > 0 else None
            self.matchers[whitelist_key] = {
                'exact': exact_matches,
                'mismatch': mismatch_lookup,
                'length': len(next(iter(sequences)))
            }

    def _create_mismatch_lookup(self, sequences: Set[str]) -> Dict[str, str]:
        """Create lookup table for sequences with allowed mismatches"""
        lookup = {}
        for seq in sequences:
            lookup[seq] = (seq, 0)  # (original sequence, mismatch count)
            if self.mismatches > 0:
                for i in range(len(seq)):
                    for base in 'ACGTN':
                        if base != seq[i]:
                            variant = seq[:i] + base + seq[i+1:]
                            if variant not in lookup:
                                lookup[variant] = (seq, 1)
        return lookup

    @lru_cache(maxsize=1024)
    def _reverse_complement(self, sequence: str) -> str:
        """Cached reverse complement computation"""
        return reverse_complement(sequence)

    def find_matches(self, sequence: str) -> List[Dict]:
        """
        Find all matching barcodes in a sequence
        Returns list of matches with positions and orientations
        """
        matches = []
        seq_len = len(sequence)
        
        for whitelist_key, matcher in self.matchers.items():
            barcode_len = matcher['length']
            
            # Scan sequence for matches in both orientations
            for start in range(seq_len - barcode_len + 1):
                candidate = sequence[start:start + barcode_len]
                
                # Check forward orientation
                match = self._check_match(candidate, matcher, whitelist_key, 'forward')
                if match:
                    matches.append({
                        'barcode': match[0],
                        'whitelist': whitelist_key,
                        'orientation': 'forward',
                        'start': start + 1,  # 1-based position
                        'end': start + barcode_len,
                        'mismatches': match[1]
                    })
                
                # Check reverse orientation
                rev_candidate = self._reverse_complement(candidate)
                match = self._check_match(rev_candidate, matcher, whitelist_key, 'reverse')
                if match:
                    matches.append({
                        'barcode': match[0],
                        'whitelist': whitelist_key,
                        'orientation': 'reverse',
                        'start': start + 1,
                        'end': start + barcode_len,
                        'mismatches': match[1]
                    })
        
        return sorted(matches, key=lambda x: x['start'])

    def _check_match(self, candidate: str, matcher: Dict, whitelist: str, orientation: str) -> Tuple[str, int]:
        """Check if candidate matches any barcode in the matcher"""
        if candidate in matcher['exact']:
            return (candidate, 0)
        if self.mismatches > 0 and matcher['mismatch'] and candidate in matcher['mismatch']:
            return matcher['mismatch'][candidate]
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
            for base in 'ACGTN':
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
                    current_run = {
                        'start': pos,
                        'bases': [base],
                        'freqs': [freq]
                    }
                else:
                    # Extend current run
                    current_run['bases'].append(base)
                    current_run['freqs'].append(freq)
            else:
                # Check if we should save the current run
                if current_run and len(current_run['bases']) >= min_length:
                    conserved_runs.append(self._finalize_run(current_run))
                current_run = None
        
        # Check final run
        if current_run and len(current_run['bases']) >= min_length:
            conserved_runs.append(self._finalize_run(current_run))
        
        return conserved_runs
    
    def _finalize_run(self, run: Dict) -> Dict:
        """Helper to finalize a conserved run with summary statistics"""
        return {
            'start': run['start'] + 1,  # Convert to 1-based position
            'end': run['start'] + len(run['bases']),
            'length': len(run['bases']),
            'sequence': ''.join(run['bases']),
            'median_freq': sorted(run['freqs'])[len(run['freqs'])//2]
        }
    
    def write_frequencies(self, output_file: str, read_label: str) -> None:
        """
        Write frequencies to output file
        
        Args:
            output_file (str): Path to output file
            read_label (str): Label for the read (e.g., 'read1' or 'read2')
        """
        frequencies = self.get_frequencies()
        with open(output_file, 'w') as f:
            # Write header
            f.write(f"read\tposition\tA\tC\tG\tT\tN\n")
            
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
Search fastq reads for barcodes in whitelists

Example:

scarecrow seed --fastqs R1.fastq.gz R2.fastq.gz\n\t--barcodes BC1:BC1.txt BC2:BC2.txt BC3:BC3.txt\n\t--out barcode_counts.csv 
---
""",
        help="Search fastq reads for barcodes",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument(
        "-f", "--fastqs", 
        nargs="+", 
        help="Pair of FASTQ files")
    subparser.add_argument(
        "-n", "--num_reads", 
        metavar="<int>",
        help=("Number of read pairs to sample [100000]"),
        type=int,
        default=100000)
    subparser.add_argument(
        "-u", "--upper_read_count", 
        metavar="<int>",
        help=("Upper read count of reads to subsample [10000000]"),
        type=int,
        default=10000000)
    subparser.add_argument(
        "-r", "--random_seed", 
        metavar="<int>",
        help=("Random seed for sampling read pairs [1234]"),
        type=int,
        default=1234)
    subparser.add_argument(
        "-c", "--barcodes",
        metavar="<string>",
        nargs='+', 
        help='Barcode whitelist files in format <barcode_name>:<whitelist_name>:<whitelist_file>\n\t(e.g. BC1:v1:barcodes1.txt BC2:v2:barcodes2.txt ...)',
    )
    subparser.add_argument(
        "-o", "--out",
        metavar="<file>",
        help=("CSV file to write barcode counts to, also serves as prefix for conserved.tsv and frequencies.tsv files."),
        type=str,
        default="./barcode_counts.csv",
    )
    subparser.add_argument(
        "-b", "--batch_size",
        metavar="<int>",
        help=("Number of read pairs per batch to process at a time [10000]"),
        type=int,
        default=10000,
    )
    subparser.add_argument(
        "-l", "--linker_base_frequency",
        metavar="<float>",
        help=("Nucleotide frequency indicative of linker (fixed) sequence [0.75]"),
        type=float,
        default=0.75,
    )
    subparser.add_argument(
        "-m", "--linker_min_length",
        metavar="<int>",
        help=("Minimum run of bases exceeding linker_base_frequency to suggest linker sequence [10]"),
        type=int,
        default=10,
    )
    subparser.add_argument(
        "-v", "--verbose",
        action='store_true',
        help='Enable verbose output [false]'
    )
    return subparser

def validate_seed_args(parser, args):
    run_seed(fastqs = [f for f in args.fastqs],
             num_reads = args.num_reads,
             random_seed = args.random_seed,
             upper_read_count = args.upper_read_count,
             barcodes = args.barcodes, 
             output_file = args.out, 
             batches = args.batch_size, 
             linker_base_frequency = args.linker_base_frequency,
             linker_min_length = args.linker_min_length,
             verbose = args.verbose)
    
@log_errors
def parse_seed_arguments(barcode_args):
    """
    Parse seed arguments from command line.
    
    Args:
        barcode_args (List[str]): List of barcode arguments in format 'KEY:WHITELIST:FILE'
    
    Returns:
        Dict[str, List[str]]: Dictionary of barcodes with keys as region identifiers
    """
    logger = logging.getLogger('scarecrow')

    expected_barcodes = {}
    
    for arg in barcode_args:
        try:
            # Split the argument into key, whitelist and file path
            key, label, file_path = arg.split(':')
            
            # Read barcodes from the file
            try:
                barcodes = read_barcode_file(file_path)
            except FileNotFoundError:
                raise Exception(f"Barcode file not found: {file_path}")

            expected_barcodes[key,label] = barcodes
            logger.info(f"Loaded {len(barcodes)} barcodes for barcode '{key}' from whitelist '{label}' file '{file_path}'")
                                    
        except ValueError:
            logger.error(f"Invalid barcode argument format: {arg}. Use 'KEY:FILE'")
    
    return expected_barcodes

@log_errors
def read_barcode_file(file_path):
    """
    Read barcode sequences from a text file.
    
    Args:
        file_path (str): Path to the barcode file
    
    Returns:
        List[str]: List of unique barcode sequences
    """
    logger = logging.getLogger('scarecrow')

    try:
        with open(file_path, 'r') as f:
            # Read lines, strip whitespace, remove empty lines
            barcodes = [line.strip() for line in f if line.strip()]
        
        # Remove duplicates while preserving order
        unique_barcodes = list(dict.fromkeys(barcodes))
        
        if not unique_barcodes:
            logger.warning(f"No barcodes found in file: {file_path}")
        
        return unique_barcodes
    
    except FileNotFoundError:
        logger.error(f"Barcode file not found: {file_path}")
        sys.exit(1)
        return []
    except Exception as e:
        logger.error(f"Error reading barcode file {file_path}: {e}")
        sys.exit(1)

@log_errors
def process_read_batch(read_batch: List[Tuple], 
                       verbose: bool = False, 
                       matcher: BarcodeMatcherOptimized = None) -> List[Dict]:
    """
    Process a batch of reads efficiently
    """
    logger = logging.getLogger('scarecrow')
    results = []
    
    for read1, read2 in read_batch:
        read_pair_info = {
            'read1': {
                'header': read1.name,
                'regions': matcher.find_matches(read1.sequence)
            },
            'read2': {
                'header': read2.name,
                'regions': matcher.find_matches(read2.sequence)
            }
        }
        results.append(read_pair_info)
        if verbose:
            logger.info(f"Read 1\t{read1.name}\n{read1.sequence}\n{read_pair_info['read1']['regions']}")
            logger.info(f"Read 2\t{read2.name}\n{read2.sequence}\n{read_pair_info['read2']['regions']}")
    
    return results

@log_errors
def run_seed(
    fastqs: List[str] = None,
    num_reads: int = 0,
    random_seed: int = 1234,
    upper_read_count: int = 10000000,
    barcodes: List[str] = None,
    output_file: str = None,
    batches: int = 10000,
    linker_base_frequency: float = 0.75,
    linker_min_length: int = 10,
    verbose: bool = False
) -> None:
    """
    Optimized implementation of seed functionality
    """
    import random
    random.seed(random_seed)

    # Setup logging
    logfile = f'./scarecrow_seed_{generate_random_string()}.log'
    logger = setup_logger(logfile)
    logger.info(f"logfile: '{logfile}'")

    # Parse barcode whitelists
    expected_barcodes = parse_seed_arguments(barcodes)
    
    # Initialize optimized matcher
    matcher = BarcodeMatcherOptimized(
        barcode_sequences={k: set(v) for k, v in expected_barcodes.items()},
        mismatches=1
    )
    
    # Initialize sequence analyzers for each read
    read1_analyzer = SequenceFrequencyAnalyzer()
    read2_analyzer = SequenceFrequencyAnalyzer()

    # If subsetting FASTQ, first get total read count
    if num_reads > 0:
        if upper_read_count == 0:
            logger.info(f"For subsetting, getting total read count to generate random indices for sampling.")
            upper_read_count = sum(1 for _ in pysam.FastxFile(fastqs[0]))
        sample_indices = set(random.sample(range(upper_read_count), min(num_reads, upper_read_count)))
        logger.info(f"Selected {len(sample_indices)} random reads out of an upper limit of {upper_read_count} reads")
        logger.info(f"Random seed used: {random_seed}")

    # Process files with minimal overhead
    with pysam.FastxFile(fastqs[0]) as r1, \
         pysam.FastxFile(fastqs[1]) as r2, \
         open(output_file, 'w') as out_csv:
        
        # Write header
        out_csv.write("read\tname\tbarcode_whitelist\tbarcode\torientation\tstart\tend\tmismatches\n")
        
        # Create batches efficiently
        current_batch = []        
        for idx, (read1, read2) in enumerate(zip(r1,r2)):
            if num_reads > 0: 
                # Current read index exceeds max index in subsample so exit
                if idx > max(sample_indices):
                    break
                # Current read index not in subsample index so skip to next read
                if idx not in sample_indices:
                    continue

            # Add sequences to analyzers
            read1_analyzer.add_sequence(read1.sequence)
            read2_analyzer.add_sequence(read2.sequence)

            current_batch.append((read1, read2))
            
            if len(current_batch) >= batches:

                # Process batch
                results = process_read_batch(current_batch, verbose, matcher)
                
                # Write results
                for result in results:
                    write_batch_results(result, out_csv)
                
                current_batch = []
        
        # Process remaining reads
        if current_batch:
            results = process_read_batch(current_batch, verbose, matcher)
            for result in results:
                write_batch_results(result, out_csv)
        
    # Write frequency analysis results
    freq_output_base = os.path.splitext(output_file)[0]
    read1_analyzer.write_frequencies(f"{freq_output_base}_read1_frequencies.tsv", "read1")
    read2_analyzer.write_frequencies(f"{freq_output_base}_read2_frequencies.tsv", "read2")
    logger.info(f"Sequence base frequency results written to:\n\t{freq_output_base}_read1_frequencies.tsv\n\t{freq_output_base}_read2_frequencies.tsv")
    
    # Find and write conserved runs
    with open(f"{freq_output_base}_conserved.tsv", 'w') as f:
        f.write("read\tstart\tend\tlength\tsequence\tmedian_frequency\n")
        
        # Analyze read1
        runs = read1_analyzer.find_conserved_runs(linker_base_frequency, linker_min_length)
        for run in runs:
            f.write(f"read1\t{run['start']}\t{run['end']}\t{run['length']}\t{run['sequence']}\t{run['median_freq']:.4f}\n")
            logger.info(f"Possible linker sequence:\n\t'read1:{run['start']}-{run['end']}'\t'{run['sequence']}'")
        
        # Analyze read2
        runs = read2_analyzer.find_conserved_runs(linker_base_frequency, linker_min_length)
        for run in runs:
            f.write(f"read2\t{run['start']}\t{run['end']}\t{run['length']}\t{run['sequence']}\t{run['median_freq']:.4f}\n")
            logger.info(f"Possible linker sequence:\n\t'read2:{run['start']}-{run['end']}'\t'{run['sequence']}'")
        
    logger.info(f"Conserved sequence analysis results written to:\n\t{freq_output_base}_conserved.tsv")


def write_batch_results(read_pair: Dict, output_handler) -> None:
    """
    Write batch results to output file efficiently
    """
    for read_key in ['read1', 'read2']:
        read_info = read_pair[read_key]
        header = read_info['header']
        
        for region in read_info['regions']:
            output_handler.write(
                f"{read_key}\t{header}\t{region['whitelist']}\t{region['barcode']}\t"
                f"{region['orientation']}\t{region['start']}\t{region['end']}\t{region['mismatches']}\n"
            )