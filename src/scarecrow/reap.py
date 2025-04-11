# -*- coding: utf-8 -*-
"""
#!/usr/bin/env python3
@author: David Wragg
"""

import csv
import gc
import gzip as gz
import itertools
import json
import logging
import multiprocessing as mp
import os
import pandas as pd
import pysam
import re
import shutil
from argparse import RawTextHelpFormatter
from collections import Counter, defaultdict
from functools import lru_cache
from multiprocessing import Process, Queue, Value
from scarecrow import __version__
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import (
    generate_random_string,
    reverse_complement,
    get_process_memory_usage,
)
from scarecrow.encode import BarcodeMatcherAhoCorasick
from typing import List, Tuple, Optional, Dict


class BarcodeMatcherOptimized:
    __slots__ = [
        "matchers",
        "mismatches",
        "base_quality_threshold",
        "verbose",
        "trie_matcher",
        "logger",
    ]

    def __init__(
        self,
        barcode_files: Dict[str, str],
        mismatches: int = 0,
        base_quality_threshold: int = None,
        verbose: bool = False,
    ):
        """
        Initialize matcher with text files containing barcodes (one per line)
        """
        self.matchers = {}
        self.mismatches = mismatches
        self.base_quality_threshold = base_quality_threshold
        self.verbose = verbose
        self.trie_matcher = None
        self.logger = setup_worker_logger()
        self.logger.info("Initializing BarcodeMatcher")

        # Load and process barcode data from text files
        for whitelist, file in barcode_files.items():
            if file.endswith(".txt"):
                self.logger.info(f"Preparing whitelist {whitelist} from '{file}'")
                # Read barcodes and generate mismatch dictionary
                with open(file, "r") as f:
                    barcodes = [line.strip() for line in f if line.strip()]

                # Create exact matches set
                exact_matches = set(barcodes)

                # Generate mismatch lookup
                mismatch_lookup = {}
                if self.mismatches > 0:
                    for barcode in barcodes:
                        variants = generate_mismatches(barcode, self.mismatches)
                        mismatch_lookup[barcode] = {}
                        for n in range(1, self.mismatches + 1):
                            if n in variants:
                                mismatch_lookup[barcode][n] = set(
                                    variants[n].keys()
                                )  # Store variant sequences

                barcode_len = len(next(iter(exact_matches)))
                self.matchers[whitelist] = {
                    "exact": exact_matches,
                    "mismatches": mismatch_lookup,
                    "length": barcode_len,
                }

                if self.verbose:
                    self.logger.info(
                        f"Loaded {len(exact_matches)} exact matches and {len(mismatch_lookup)} variants for whitelist '{whitelist}'"
                    )
                    # self.logger.info(f"'{whitelist}' exact matches:\n{exact_matches}\n")
                    # self.logger.info(f"'{whitelist}' mismatches:\n{mismatch_lookup}\n")
                    # self.logger.info(f"Matchers: {self.matchers.keys()}")

            elif file.endswith(".pkl.gz"):
                # Initialize the trie matcher only once
                if self.trie_matcher is None:
                    if os.path.exists(file):
                        self.trie_matcher = BarcodeMatcherAhoCorasick(
                            barcode_sequences={},
                            pickle_file=file,
                            mismatches=mismatches,
                        )
                        self.logger.info("Imported pickle for worker")
                    else:
                        raise FileNotFoundError(
                            f"Barcode whitelist file not found: {file}"
                        )
                else:
                    # Load additional trie files into the existing matcher
                    if os.path.exists(file):
                        self.trie_matcher._load_trie(file, mismatches=mismatches)
                    else:
                        raise FileNotFoundError(
                            f"Barcode whitelist file not found: {file}"
                        )

    def find_match(
        self,
        sequence: str,
        quality_scores: str,
        whitelist: str,
        orientation: str,
        original_start: int,
        original_end: int,
        jitter: int,
        base_quality: int,
    ) -> Tuple[str, int, int]:
        """
        Find best matching barcode using precomputed exact and mismatch data
        Returns: Tuple of (barcode, mismatch_count, position)
        Always returns a tuple, using 'null' for no match
        """

        # Apply quality filtering if needed
        if base_quality is not None:
            sequence = filter_low_quality_bases(sequence, quality_scores, base_quality)[
                0
            ]

        # Handle reverse orientation
        if orientation == "reverse":
            sequence = reverse_complement(sequence)

        # Get sub sequences with jitter
        sub_sequence = self._get_sequence_with_jitter(
            sequence, original_start, original_end, jitter
        )

        if self.trie_matcher:
            # self.logger.info(f"-> start:{original_start} end:{original_end} jitter:{jitter}")

            # Use the Aho-Corasick trie for matching
            matches = self.trie_matcher.find_matches(
                sub_sequence, whitelist, orientation, original_start
            )

            # Filter matches based on jitter and select best match
            if matches:
                # self.logger.info(f"\n> reap matches: {matches}")
                filtered_matches = [
                    match
                    for match in matches
                    if abs(match.start - original_start) <= jitter
                ]
                if filtered_matches:
                    # Sort by number of mismatches and distance from expected start
                    filtered_matches.sort(key=lambda x: (x.mismatches, abs(x.distance)))
                    if len(filtered_matches) == 1 or (
                        filtered_matches[0].mismatches < filtered_matches[1].mismatches
                        or filtered_matches[0].distance < filtered_matches[1].distance
                    ):
                        # Only one match or the best match is unique
                        best_match = filtered_matches[0]
                        return (
                            best_match.barcode,
                            best_match.mismatches,
                            best_match.start,
                        )
                    else:
                        # Multiple matches with the same number of mismatches and distance
                        # self.logger.info("Multiple equidistant-error matches")
                        return "NNNNNNNN", -1, "N"
                else:
                    # No match found within the jitter range
                    # self.logger.info("No match found within jitter range")
                    return "NNNNNNNN", -1, "N"
            else:
                # No match found
                # self.logger.info("No match")
                return "NNNNNNNN", -1, "N"

        else:
            # Default to set-based method
            # NULL barcode
            NA_barcode = 'N' * len(sub_sequence[0][0])
            
            # Check whitelist available
            if whitelist not in self.matchers:
                self.logger.warning(
                    f"Whitelist '{whitelist}' not found in available whitelists: {list(self.matchers.keys())}"
                )
                return NA_barcode, -1, "N"

            # self.logger.info(f"-> start:{original_start} end:{original_end} jitter:{jitter}")
            # self.logger.info(f"{sub_sequence}")

            # First pass: Look for exact matches with best position
            exact_matches = []
            for seq, pos in sub_sequence:
                if seq in self.matchers[whitelist]["exact"]:
                    pos_distance = abs(pos - original_start)
                    exact_matches.append((seq, pos_distance, pos))

            # If we found exact matches, return the one closest to the expected position
            if exact_matches:
                exact_matches.sort(
                    key=lambda x: x[1]
                )  # Sort by distance from expected start
                # self.logger.info(f"\nExact matches: {exact_matches}")
                if len(exact_matches) == 1 or (
                    exact_matches[0][1] < exact_matches[1][1]
                ):
                    # Only one exact match
                    match = exact_matches[0]
                    return match[0], 0, match[2]
                else:
                    # Multiple exact matches with the same distance
                    # self.logger.info(f"Multiple matches")
                    return NA_barcode, -1, "N"

            # If no exact match was found, check mismatch lookup
            # self.logger.info(f"subs_sequence: {sub_sequence}")
            if self.mismatches > 0:
                mismatch_matches = []
                for seq, pos in sub_sequence:
                    # Query mismatch_lookup for all barcodes that match this sequence
                    matching_barcodes = []
                    for barcode, variants in self.matchers[whitelist][
                        "mismatches"
                    ].items():
                        for n in range(1, self.mismatches + 1):
                            if n in variants and seq in variants[n]:
                                matching_barcodes.append((barcode, n))
                                break  # No point checking n+1 if there are results form n

                    if matching_barcodes:
                        # Calculate distance from expected start
                        pos_distance = abs(pos - original_start)
                        # Add all matching barcodes to mismatch_matches
                        for barcode, n in matching_barcodes:
                            mismatch_matches.append((barcode, n, pos, pos_distance))

                if mismatch_matches:
                    # self.logger.info(f"Mismatch matches: {mismatch_matches}")
                    # Group matches by the number of mismatches and distance
                    match_groups = defaultdict(list)
                    for match in mismatch_matches:
                        key = (match[1], match[3])  # (mismatch_count, distance)
                        match_groups[key].append(match)

                    # Find the best group (fewest mismatches, smallest absolute distance)
                    best_key = min(match_groups.keys(), key=lambda x: (x[0], abs(x[1])))
                    best_matches = match_groups[best_key]

                    if len(best_matches) == 1:
                        # Only one match in the best group
                        match = best_matches[0]
                        # self.logger.info(f"Selected best mismatch match: {match[0]} at position {match[2]}")
                        return match[0], match[1], match[2]
                    else:
                        # Multiple matches in the best group
                        # self.logger.info("Multiple mismatch matches with the same distance, returning no match")
                        return NA_barcode, -1, "N"

            # No match found
            # self.logger.info("No match found")
            #null_match = "N" * len(next(iter(self.matchers[whitelist]["exact"])))
            return NA_barcode, -1, "N"

    def _get_sequence_with_jitter(
        self, full_sequence: str, start: int, end: int, jitter: int
    ) -> List[Tuple[str, int]]:
        """Generate possible sequences with position jitter"""
        barcode_length = end - start + 1
        sequences = []

        min_start = max(0, start - jitter - 1)
        max_start = min(len(full_sequence) - barcode_length + 1, start + jitter)

        for adj_start in range(min_start, max_start):
            seq = full_sequence[adj_start : adj_start + barcode_length]
            if len(seq) == barcode_length:
                sequences.append((seq, adj_start + 1))

        # Handle cases where the start position minus jitter is less than 1
        if start - jitter < 1:
            for clip_count in range(1, jitter + 1):
                # Extract the original sequence clipped at the end (start is - 1 for Python 0-based position)
                original_seq = full_sequence[
                    (start - 1) : (start - 1 + barcode_length - clip_count)
                ]

                # Insert 'N's at the start
                clipped_seq = "N" * clip_count + original_seq

                # if len(clipped_seq) == barcode_length:
                sequences.append(
                    (clipped_seq, start - clip_count - 1)
                )  # -1 skips position 0

        return sequences


# Add a new class to track statistics
class BarcodeStats:
    __slots__ = ["mismatch_counts", "position_counts"]

    def __init__(self):
        self.mismatch_counts = Counter()
        self.position_counts = Counter()

    def update(self, mismatches: List[str], positions: List[str]):
        # Convert mismatches to integers, handling 'NA' values as -1
        mismatch_values = [int(m) if m != "NA" else -1 for m in mismatches]

        # If any mismatch is negative, sum only the negative values
        if any(m < 0 for m in mismatch_values):
            total_mismatches = sum(m for m in mismatch_values if m < 0)
        else:
            total_mismatches = sum(mismatch_values)

        self.mismatch_counts[total_mismatches] += 1

        # Track each position
        for pos in positions:
            # if not pos == 'N' and pos > '0':  # Exclude invalid positions
            self.position_counts[pos] += 1

    def write_stats(self, output_prefix: str):
        # Write mismatch counts
        with open(f"{output_prefix}_mismatch_stats.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["mismatches", "count"])
            for mismatches, count in sorted(self.mismatch_counts.items()):
                writer.writerow([mismatches, count])

        # Write position counts
        with open(f"{output_prefix}_position_stats.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["position", "count"])
            for position, count in sorted(self.position_counts.items()):
                writer.writerow([position, count])


def generate_mismatches(
    sequence: str, max_mismatches: int
) -> Dict[int, Dict[str, str]]:
    """Generate all possible sequences with up to max_mismatches mismatches"""
    bases = {"A", "C", "G", "T", "N"}
    mismatch_dict = {i: {} for i in range(1, max_mismatches + 1)}
    seq_length = len(sequence)

    for num_mismatches in range(1, max_mismatches + 1):
        for positions in itertools.combinations(range(seq_length), num_mismatches):
            for replacements in itertools.product(bases, repeat=num_mismatches):
                # Convert to list for mutation
                mutated_seq = list(sequence)
                valid_mutation = False

                # Apply mutations
                for pos, new_base in zip(positions, replacements):
                    if mutated_seq[pos] != new_base:  # Ensure we introduce a mismatch
                        mutated_seq[pos] = new_base
                        valid_mutation = True

                # Only add if we actually introduced a mismatch
                if valid_mutation:
                    mutated_str = "".join(mutated_seq)
                    if mutated_str not in mismatch_dict[num_mismatches]:
                        mismatch_dict[num_mismatches][mutated_str] = sequence

    return mismatch_dict


def parser_reap(parser):
    subparser = parser.add_parser(
        "reap",
        description="""
Extract target sequence range from FASTQ file index, and use the barcode profiles generated by
scarecrow harvest to identify and correct barcodes. Output sequence and corrected barcodes to 
either SAM or interleaved FASTQ file.

Example:

scarecrow reap --threads16\n\t--fastqs R1.fastq.gz R2.fastq.gz\n\t--barcode_positions barcode_positions.csv\n\t--barcodes\tBC1:v1_whitelist:bc1_whitelist.txt\n\t\t\tBC2:v2_whitelist:bc2_whitelist.txt\n\t\t\tBC3:v1_whitelist:bc3_whitelist.txt\n\t--read1 0-64\n\t--out extracted_sequences\n\t--out_sam
---
""",
        epilog="The --barcodes <name> must match the barcode_whitelist values in the --barcode_positions file.",
        help="Extract target sequence, perform barcode correction, write to SAM or FASTQ",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument(
        "--fastqs", nargs="+", required=True, help="FASTQ files"
    )
    subparser.add_argument(
        "-o",
        "--out",
        metavar="<file>",
        help=("Prefix for output"),
        type=str,
        default="extracted",
    )
    # Add a mutually exclusive group for output format
    out_format = subparser.add_mutually_exclusive_group(required=False)
    out_format.add_argument(
        "--out_sam", action="store_true", help="Write output to SAM format [default]"
    )
    out_format.add_argument(
        "--out_fastq", action="store_true", help="Write output to FASTQ format"
    )
    subparser.add_argument(
        "-p",
        "--barcode_positions",
        metavar="<file>",
        help=("File containing barcode positions, output by scarecrow harvest"),
        type=str,
        required=True,
        default=[],
    )
    subparser.add_argument(
        "-j",
        "--jitter",
        metavar="<int>",
        type=int,
        default=2,
        help="Barcode position jitter [2]",
    )
    subparser.add_argument(
        "-m",
        "--mismatches",
        metavar="<int>",
        type=int,
        default=1,
        help="Number of allowed mismatches in barcode [1]",
    )
    subparser.add_argument(
        "-q",
        "--base_quality",
        metavar="<int>",
        type=int,
        default=None,
        help="Minimum base quality filter [None]",
    )
    subparser.add_argument(
        "-x",
        "--extract",
        metavar="<range>",
        help=("Sequence range to extract <read>:<range> (e.g. 1:1-64)"),
        type=str,
        required=True,
        default=None,
    )
    subparser.add_argument(
        "-u",
        "--umi",
        metavar="<range>",
        help=("Sequence range to extract for UMI <read>:<range> (e.g. 2:1-10)"),
        type=str,
        default=None,
    )
    subparser.add_argument(
        "-c",
        "--barcodes",
        metavar="<string>",
        nargs="+",
        required=True,
        help="Barcode whitelist files in format <barcode_name>:<whitelist_name>:<whitelist_file>\n\t(e.g. BC1:v1:barcodes1.txt BC2:v2:barcodes2.txt ...)",
    )
    subparser.add_argument(
        "-b",
        "--batch_size",
        metavar="<int>",
        help=("Number of read pairs per batch to process at a time [10000]"),
        type=int,
        default=10000,
    )
    subparser.add_argument(
        "-@",
        "--threads",
        metavar="<int>",
        help=("Number of processing threads [1]"),
        type=int,
        default=1,
    )
    subparser.add_argument(
        "-r",
        "--barcode_reverse_order",
        action="store_true",
        help="Reverse retrieval order of barcodes in barcode positions file [false]",
    )
    subparser.add_argument(
        "-z", "--gzip", 
        action="store_true", 
        help="Compress (gzip) fastq output [false]"
    )
    subparser.add_argument(
        "-v", "--verbose", 
        action="store_true", 
        help="Enable verbose output [false]"
    )
    subparser.add_argument(
        "-s", "--sift", 
        action="store_true", 
        help="Sift output to skip writing reads with any invalid barcodes [false]"
    )
    # Add a mutually exclusive group for library defaults
    lib_format = subparser.add_mutually_exclusive_group(required=False)
    lib_format.add_argument(
        "--scale_rna",
        action="store_true",
        help="Use defaults for Scale RNA data [--jitter 1] and [--mismatches 2]",
    )
    lib_format.add_argument(
        "--evercode_wtv2",
        action="store_true",
        help="Use defaults for Parse Evercode WTv2 data [--jitter 2] and [--mismatches 2]",
    )
    return subparser


def validate_reap_args(parser, args) -> None:
    """
    Validate arguments
    """
    # Global logger setup
    logfile = "{}_{}.{}".format("./scarecrow_reap", generate_random_string(), "log")
    logger = setup_logger(logfile)
    logger.info(f"scarecrow version {__version__}")
    logger.info(f"logfile: '{logfile}'")

    # Check for defaults
    if args.evercode_wtv2:
        logger.info(f"Using Parse Evercode WTv2 defaults: --jitter 2 --mismatches 2")
        args.jitter = 2
        args.mismatches = 2
    if args.scale_rna:
        logger.info(f"Using Scale RNA defaults: --jitter 1 --mismatches 2")
        args.jitter = 1
        args.mismatches = 2

    logger.info(f"{args}\n")

    # Check input files exist
    missing_files = []
    missing_files.extend(f for f in args.fastqs if not os.path.exists(f))
    if not os.path.exists(args.barcode_positions):
        missing_files.append(args.barcode_positions)
    for barcode in args.barcodes:
        key, label, file = barcode.split(":")
        if not os.path.exists(file):
            missing_files.append(file)

    # Log any missing files and raise error
    if missing_files:
        logger.error(f"The following files were not found:\n{'\n'.join(missing_files)}")
        raise FileNotFoundError

    # Check output path exists
    outpath = os.path.dirname(args.out)
    if not os.path.exists(outpath) and outpath != "":
        logger.error(f"Output directory {outpath} does not exist")
        raise FileNotFoundError
    

    run_reap(
        fastqs=[f for f in args.fastqs],
        barcode_positions=args.barcode_positions,
        barcode_reverse_order=args.barcode_reverse_order,
        output=args.out,
        extract=args.extract,
        umi=args.umi,
        barcodes=args.barcodes,
        jitter=args.jitter,
        mismatches=args.mismatches,
        base_quality=args.base_quality,
        batches=args.batch_size,
        threads=args.threads,
        FASTQ=args.out_fastq,
        SAM=args.out_sam,
        verbose=args.verbose,
        gzip=args.gzip,
        sift=args.sift,
        args_string=" ".join(
            f"--{k} {v}" for k, v in vars(args).items() if v is not None
        ),
    )


@log_errors
def run_reap(
    fastqs: List[str],
    barcode_positions: str = None,
    barcode_reverse_order: bool = False,
    output: str = "extracted.fastq",
    extract: str = None,
    umi: Optional[str] = None,
    barcodes: List[str] = None,
    jitter: int = 2,
    mismatches: int = 1,
    base_quality: int = None,
    batches: int = 10000,
    threads: int = 1,
    FASTQ: bool = False,
    SAM: bool = True,
    verbose: bool = False,
    gzip: bool = False,
    sift: bool = False,
    args_string: str = None,
) -> None:
    """
    Main function to extract sequences with barcode headers
    """
    logger = logging.getLogger("scarecrow")
    total_mb, total_gb = get_process_memory_usage()
    logger.info(
        f"Current memory usage (main process + children): {total_mb:.2f} MB ({total_gb:.2f} GB)"
    )

    # Extract barcodes and convert whitelist to set
    barcode_files = parse_seed_arguments(barcodes)
    if verbose:
        for whitelist, filename in barcode_files.items():
            logger.info(f"{whitelist}: {filename}")

    total_mb, total_gb = get_process_memory_usage()
    logger.info(
        f"Current memory usage (main process + children): {total_mb:.2f} MB ({total_gb:.2f} GB)"
    )

    # Default output to SAM format if not specified
    if FASTQ is False and SAM is False:
        logger.info("Defaulting to SAM file output")
        SAM = True

    # Append suffix to output prefix
    if FASTQ:
        outfile = f"{output}.fastq"
    if SAM:
        outfile = f"{output}.sam"
        with open(outfile, "w") as file:
            file.write("@HD\tVN:1.6\n")
            file.write(
                f"@PG\tID:reap\tPN:scarecrow\tVN:{__version__}\tDS:{args_string}\n"
            )

    logger.info(f"Results will be written to '{outfile}'")

    # Extract sequences
    extract_sequences(
        fastq_files=[f for f in fastqs],
        barcode_positions_file=barcode_positions,
        barcode_reverse_order=barcode_reverse_order,
        barcode_files=barcode_files,
        output=outfile,
        extract=extract,
        umi=umi,
        jitter=jitter,
        mismatches=mismatches,
        base_quality=base_quality,
        batch_size=batches,
        threads=threads,
        FASTQ=FASTQ,
        SAM=SAM,
        verbose=verbose,
        sift=sift
    )

    # Generate JSON file if FASTQ output is enabled
    if FASTQ:
        # Parse extract and UMI ranges
        extract_index, extract_range = extract.split(":")
        extract_range = parse_range(extract_range)
        if umi is not None:
            umi_index, umi_range = umi.split(":")
            umi_range = parse_range(umi_range)
        else:
            umi_range = None

        # Prepare barcode configurations
        barcode_positions_df = pd.read_csv(barcode_positions)
        if barcode_reverse_order:
            barcode_positions_df = barcode_positions_df[::-1].reset_index(drop=True)
        barcode_configs = prepare_barcode_configs(barcode_positions_df, jitter)

        # Generate JSON
        json_file = generate_fastq_json(
            barcode_configs = barcode_configs,
            barcode_files = barcode_files,
            umi_range = umi_range,
            output_prefix = output,
        )
        logger.info(f"Generated JSON file: {json_file}")

    # gzip
    if FASTQ and gzip:
        logger.info(f"Compressing '{outfile}'")
        with open(outfile, "rb") as f_in, gz.open(outfile + ".gz", "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(outfile)

    logger.info("Finished!")


@log_errors
def process_read_batch(
    read_batch: List[Tuple],
    barcode_configs: List[Dict],
    matcher: BarcodeMatcherOptimized,
    read_range: Tuple[int, int],
    read_index: int,
    umi_index: int,
    umi_range: Tuple[int, int],
    stats: BarcodeStats,
    base_quality: int = None,
    jitter: int = 3,
    FASTQ: bool = False,
    SAM: bool = True,
    verbose: bool = False,
    sift: bool = False
) -> List[str]:
    """
    Modified process_read_batch to handle any number of FASTQ files.
    """
    logger = setup_worker_logger()
    if verbose:
        logger.info(f"Processing batch of {len(read_batch)} reads")
    output_entries = []
    read_count = len(read_batch)

    for i, reads in enumerate(read_batch):
        if i % 100000 == 0:
            logger.debug(f"Processing read {i}/{read_count}")

        original_barcodes = []  # CR
        barcode_qualities = []  # CY
        matched_barcodes = []  # CB
        positions = []  # XP
        mismatches = []  # XM
        matched_qualities = [] # XQ

        for config in barcode_configs:
            whitelist = config["whitelist"]
            file_index = config["file_index"]
            seq = reads[file_index].sequence
            qual = reads[file_index].quality
            start, end = config["start"], config["end"]

            # Extract original barcode and its quality scores
            original_barcode = seq[start - 1 : end]
            barcode_quality = qual[start - 1 : end]
            original_barcodes.append(original_barcode)
            barcode_qualities.append(barcode_quality)            

            if verbose:
                logger.info(
                    f"Read: {reads[file_index].name} {reads[file_index].comment}"
                )
                logger.info(f"Sequence: {seq}")
                logger.info(
                    f"Looking for {original_barcode} + {barcode_quality} in range {start}-{end} with jitter {jitter}"
                )
            # logger.info(f"{reads[config['file_index']].name}")
            if matcher.trie_matcher or whitelist in matcher.matchers.keys():
                matched_barcode, mismatch_count, adj_position = matcher.find_match(
                    seq,
                    reads[file_index].quality,
                    whitelist,
                    config["orientation"],
                    start,
                    end,
                    jitter,
                    base_quality,
                )

                if verbose:
                    logger.info(
                        f"Matched barcode: {matched_barcode} with {mismatch_count} mismatches at position {adj_position}"
                    )

                matched_barcodes.append(matched_barcode)
                # Get qualities of matched barcode
                qualities = 'N' * len(matched_barcode)
                if isinstance(adj_position, int):
                    if adj_position >= 0:
                        qualities = qual[ adj_position - 1 : adj_position - 1 + len(matched_barcode) ]
                    else:
                        qualities = ('').join(['N' * abs(adj_position), qual[0:len(matched_barcode) - abs(adj_position)]])
                matched_qualities.append(qualities)                
                # Update stats
                positions.append(str(adj_position))
                mismatches.append(str(mismatch_count))
            else:
                logger.warning(
                    f"Whitelist {whitelist} not found in matcher {matcher.matchers.keys()}"
                )
                matched_barcodes.append("null")
                matched_qualities.append("null")
                positions.append(str(start))
                mismatches.append("NA")

        # Update statistics
        stats.update(mismatches, positions)

        # Check if we need to sift (skip writing reads with any barcode with an N)
        # if this test is true, then the loop continues with the next iteration, skipping FASTQ and SAM output
        if sift and any('N' in barcode for barcode in matched_barcodes):
            continue

        # Create output entry
        source_entry = reads[read_index]

        # Apply base quality filtering to extracted sequence
        if base_quality is not None:
            filtered_seq, filtered_qual = filter_low_quality_bases(
                source_entry.sequence[read_range[0] : read_range[1]],
                source_entry.quality[read_range[0] : read_range[1]],
                base_quality,
            )
        else:
            filtered_seq = source_entry.sequence[read_range[0] : read_range[1]]
            filtered_qual = source_entry.quality[read_range[0] : read_range[1]]

        #  Output interleaved FASTQ
        if FASTQ:
            # Check end of read comment has the paired read identifier (i.e. /1 or /2)            
            if not re.search(r'/[0-9]$', source_entry.comment):
                # If it doesn't end with /number, add /
                source_entry.comment += "/"
            else:
                # If it ends with /number, remove number
                source_entry.comment = re.sub(r'/[0-9]$', '/', source_entry.comment)

            # R1
            r1_header = f"@{source_entry.name} {source_entry.comment}1"
            r1_barcodes = f"{('').join(matched_barcodes)}"
            r1_quality = f"{('').join(matched_qualities)}"
            if umi_index is not None:
                umi_sequence = reads[umi_index].sequence[umi_range[0] : umi_range[1]]
                umi_quality = reads[umi_index].quality[umi_range[0] : umi_range[1]]
                r1_barcodes += f"{umi_sequence}"
                r1_quality += f"{umi_quality}"
                
            # R2
            r2_header = f"@{source_entry.name} {source_entry.comment}2"
            output_entries.append(f"{r1_header}\n{r1_barcodes}\n+\n{r1_quality}\n{r2_header}\n{filtered_seq}\n+\n{filtered_qual}\n")


        # Output SAM
        if SAM:
            qname = source_entry.name  # QNAME
            flag = "4"  # FLAG (4 means unaligned)
            rname = "*"  # RNAME (unaligned, so *)
            pos = "0"  # POS (unaligned, so 0)
            mapq = "255"  # MAPQ (255 means unavailable)
            cigar = "*"  # CIGAR (unaligned, so *)
            rnext = "*"  # RNEXT (unaligned, so *)
            pnext = "0"  # PNEXT (unaligned, so 0)
            tlen = "0"  # TLEN (unaligned, so 0)
            seq = filtered_seq  # SEQ
            qual = filtered_qual  # QUAL

            # Optional tags
            tags = []
            tags.append(f"CR:Z:{'_'.join(original_barcodes)}")
            tags.append(f"CY:Z:{'_'.join(barcode_qualities)}")
            tags.append(f"CB:Z:{'_'.join(matched_barcodes)}")
            tags.append(f"XQ:Z:{'_'.join(matched_qualities)}")
            tags.append(f"XP:Z:{'_'.join(positions)}")
            tags.append(f"XM:Z:{'_'.join(mismatches)}")

            """
            Need to check if UMI is on a read with barcodes
                if so, is it downstream of barcodes
                if barcode is jittered then UMI needs to be jittered
            """

            # Add UMI information if specified
            if umi_index is not None:
                umi_seq = reads[umi_index].sequence[umi_range[0] : umi_range[1]]
                umi_qual = reads[umi_index].quality[umi_range[0] : umi_range[1]]
                tags.append(f"UR:Z:{umi_seq}")
                tags.append(f"UY:Z:{umi_qual}")

            # Combine all fields into a SAM line
            sam_line = (
                "\t".join(
                    [
                        qname,
                        flag,
                        rname,
                        pos,
                        mapq,
                        cigar,
                        rnext,
                        pnext,
                        tlen,
                        seq,
                        qual,
                    ]
                    + tags
                )
                + "\n"
            )
            output_entries.append(sam_line)

    if verbose:
        logger.info(f"Completed processing batch of {read_count} reads")
    return output_entries, read_count, stats


@log_errors
def extract_sequences(
    fastq_files: List[str] = None,
    barcode_positions_file: str = None,
    barcode_reverse_order: bool = False,
    barcode_files: Dict[str, str] = None,
    output: str = "extracted.fastq.gz",
    extract: str = None,
    umi: Optional[str] = None,
    jitter: int = 3,
    mismatches: int = 1,
    base_quality: int = None,
    batch_size: int = 100000,
    threads: Optional[int] = None,
    FASTQ: bool = False,
    SAM: bool = True,
    verbose: bool = False,
    sift: bool = False
) -> None:
    """
    Modified to track stats
    """
    logger = setup_worker_logger()

    # Initialize configurations
    barcode_positions = pd.read_csv(barcode_positions_file)
    if barcode_reverse_order:
        barcode_positions = barcode_positions[::-1].reset_index(drop=True)
    barcode_configs = prepare_barcode_configs(barcode_positions, jitter)
    logger.info(f"Barcode configs\n{barcode_configs}\n")

    # Extract range
    extract_index, extract_range = extract.split(":")
    extract_index = int(extract_index) - 1
    extract_range = parse_range(extract_range)
    logger.info(f"FASTQ sequence range to extract: '{extract}'")

    # UMI range
    if umi is not None:
        umi_index, umi_range = umi.split(":")
        umi_index = int(umi_index) - 1
        umi_range = parse_range(umi_range)
        logger.info(f"UMI sequence range to extract: '{umi}'")
    else:
        umi_index = None
        umi_range = None

    # Set number of threads to use
    if threads is None:
        threads = min(mp.cpu_count() - 1, 8)
    else:
        threads = min(threads, mp.cpu_count())
    logger.info(f"Using {threads} threads")

    # Report memory usage before starting read processing
    total_mb, total_gb = get_process_memory_usage()
    logger.info(
        f"Current memory usage (main process + children): {total_mb:.2f} MB ({total_gb:.2f} GB)"
    )

    # Constant arguments
    constant_args = (
        barcode_configs,
        extract_range,
        extract_index,
        umi_index,
        umi_range,
        base_quality,
        jitter,
        mismatches,
        FASTQ,
        SAM,
        verbose,
        sift
    )

    # Create a queue for batch communication
    queue = Queue(
        maxsize=threads * 4
    )  # Limit queue size to avoid excessive memory usage
    stats_queue = Queue()

    # Initialize the shared counter for total reads
    total_reads_counter = init_shared_counter()

    # Start the producer process
    producer = Process(target=batch_producer, args=(fastq_files, batch_size, queue))
    producer.start()

    # Start worker processes
    workers = []
    for i in range(threads):
        worker_output = f"{output}_worker_{i}.sam"  # Each worker writes to its own file
        worker = Process(
            target = worker_process,
            args = (queue, barcode_files, worker_output, constant_args, stats_queue, total_reads_counter),
        )
        worker.start()
        workers.append(worker)
        logger.info(f"Started worker {i + 1}")

    # Wait for the producer to finish
    producer.join()

    # Signal workers to stop by adding sentinel values
    for _ in range(threads):
        queue.put(None)

    # Wait for all workers to finish
    for worker in workers:
        worker.join()

    # Collect stats from workers
    master_stats = BarcodeStats()
    for _ in range(threads):
        worker_stats = stats_queue.get()  # Get stats from the queue
        master_stats.mismatch_counts.update(worker_stats.mismatch_counts)
        master_stats.position_counts.update(worker_stats.position_counts)

    # Write final statistics
    master_stats.write_stats(output)

    # Combine worker output files into a single file
    combine_worker_outputs(output, threads)
    
    # Log the final total number of reads processed
    logger.info(f"Total reads processed: {total_reads_counter.value}")


def parse_range(range_str: str) -> Tuple[int, int]:
    """
    Parse range string
    """
    start, end = map(int, range_str.split("-"))
    start = max(0, start - 1)
    return (start, end)


def parse_seed_arguments(barcodes: List[str]) -> Dict[str, str]:
    """
    Format: <barcode_name>:<whitelist_name>:<whitelist_file.json>
    Returns: Dict mapping whitelist names to JSON file paths
    """
    logger = logging.getLogger("scarecrow")
    barcode_files = {}

    if barcodes:
        for barcode in barcodes:
            key, label, file = barcode.split(":")
            logger.info(f"Processing barcode definition: '{barcode}'")
            logger.info(f"  Name: '{key}'")
            logger.info(f"  Whitelist: '{label}'")
            logger.info(f"  File: '{file}'")

            if not os.path.exists(file):
                raise FileNotFoundError(f"Barcode whitelist file not found: {file}")
            barcode_files[f"{key}:{label}"] = file

    logger.info(f"Parsed barcode whitelists: {barcode_files}")
    return barcode_files


@lru_cache(maxsize=1024)
def filter_low_quality_bases(
    sequence: str, quality: str, threshold: int
) -> Tuple[str, str]:
    """
    Efficiently filter low-quality bases with minimal overhead
    """
    # Pre-compute quality scores to avoid repeated calculation
    qual_scores = [ord(q) - 33 for q in quality]
    # Use list comprehension for efficient filtering
    filtered_seq = "".join(
        base if score >= threshold else "N"
        for base, score in zip(sequence, qual_scores)
    )
    return filtered_seq, quality


def prepare_barcode_configs(positions: pd.DataFrame, jitter: int) -> List[Dict]:
    """
    Prepare barcode configurations
    """
    return [
        {
            "index": idx,
            "file_index": row["file_index"],
            "start": row["start"],
            "end": row["end"],
            "orientation": row["orientation"],
            "whitelist": row["barcode_whitelist"],
        }
        for idx, row in positions.iterrows()
    ]


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


def batch_producer(fastq_files, batch_size, queue):
    """
    Producer function that reads multiple FASTQ files and generates batches of read tuples.
    """
    # Open all FASTQ files
    fastq_handles = [pysam.FastxFile(f) for f in fastq_files]

    # Read all files in parallel
    read_tuples = zip(*fastq_handles)
    batch = []
    for read_tuple in read_tuples:
        batch.append(read_tuple)
        if len(batch) >= batch_size:
            queue.put(batch)
            batch = []  # Reset for the next batch
    if batch:  # Yield any remaining reads
        queue.put(batch)

    # Close all FASTQ files
    for handle in fastq_handles:
        handle.close()

    queue.put(None)  # Sentinel value to signal the end of production


def init_shared_counter():
    """
    Initialize a shared counter for tracking the total number of reads processed.
    """
    return Value("i", 0)  # Shared integer value initialized to 0


def worker_process(queue, barcode_files, output_file, constant_args, stats_queue, total_reads_counter):
    """
    Worker function that processes batches from the queue.
    """
    try:
        (
            barcode_configs,
            extract_range,
            extract_index,
            umi_index,
            umi_range,
            base_quality,
            jitter,
            mismatches,
            FASTQ,
            SAM,
            verbose,
            sift
        ) = constant_args
        matcher = BarcodeMatcherOptimized(
            barcode_files=barcode_files,
            mismatches=mismatches,
            base_quality_threshold=base_quality,
            verbose=verbose,
        )

        stats = BarcodeStats()  # Create a stats object for this worker
        logger = setup_worker_logger()  # Ensure the logger is set up
        batch_count = 0  # Track the number of batches processed

        with open(output_file, "a") as outfile:
            while True:
                batch = queue.get()
                if batch is None:  # Sentinel value to stop the worker
                    break

                # Process batch
                entries, count, batch_stats = process_read_batch(
                    batch,
                    barcode_configs,
                    matcher,
                    extract_range,
                    extract_index,
                    umi_index,
                    umi_range,
                    stats,
                    base_quality,
                    jitter,
                    FASTQ,
                    SAM,
                    verbose,
                    sift
                )
                outfile.writelines(entries)

                # Update the shared counter
                with total_reads_counter.get_lock():
                    total_reads_counter.value += count
                
                # Log the total number of reads processed periodically
                if total_reads_counter.value % 1000000 == 0:
                    logger.info(f"Total reads processed: {total_reads_counter.value}")

                batch_count += 1

                # Log memory usage after processing each batch
                if batch_count % 100 == 0:  # Log every 100 batches
                    total_mb, total_gb = get_process_memory_usage()
                    logger.info(
                        f"Worker memory usage after {batch_count} batches: {total_mb:.2f} MB ({total_gb:.2f} GB)"
                    )
                    gc.collect()

        # Log final memory usage
        total_mb, total_gb = get_process_memory_usage()
        logger.info(f"Final worker memory usage: {total_mb:.2f} MB ({total_gb:.2f} GB)")

        # Put the stats object into the stats_queue
        stats_queue.put(stats)
    except Exception as e:
        logger.error(f"Worker failed with error: {str(e)}", exc_info=True)
        raise


@log_errors
def combine_worker_outputs(output, num_workers):
    """
    Combine output files from all workers into a single file.
    """
    logger = logging.getLogger("scarecrow")
    total_lines = 0
    with open(output, "w") as outfile:
        for i in range(num_workers):
            worker_output = f"{output}_worker_{i}.sam"
            try:
                with open(worker_output, "r") as infile:
                    shutil.copyfileobj(infile, outfile)
                    logger.info(f"Combined data from {worker_output}")
                os.remove(worker_output)
            except FileNotFoundError:
                logger.warning(
                    f"Intermediate file {worker_output} not found. Skipping."
                )
            except Exception as e:
                logger.error(f"Error combining {worker_output}: {str(e)}")

    logger.info(f"Total data combined into final output: {output}")


def generate_fastq_json(barcode_configs: List[Dict], barcode_files: Dict[str, str], umi_range: Optional[Tuple[int, int]], output_prefix: str) -> None:
    """
    Generate a JSON file describing the FASTQ structure for universal compatibility.
    
    Args:
        barcode_configs: List of barcode configurations from `prepare_barcode_configs`.
        umi_range: Tuple of (start, end) for the UMI sequence, or None if no UMI.
        extract_range: Tuple of (start, end) for the cDNA sequence.
        output_prefix: Prefix for the output JSON file.
    """
    # Initialize JSON to file
    json_data = {
        "description": "scarecrow",
        "barcodes": [],
        "umi:": [],
        "kallisto-bustools": []
        }
    
    # Barcode information
    barcodes = []
    kb_x = None
    star_x = None
    current_position = 0
    whitelists = []
    for config in barcode_configs:
        barcode_length = config["end"] - config["start"]
        end_position = current_position + barcode_length + 1
        #whitelists.append = barcode_files[config["whitelist"]]
        json_data["barcodes"].append({
            "range": f"1:{current_position + 1}-{end_position}",
            "whitelist": f"{barcode_files[config["whitelist"]]}"
            })
        if kb_x is None:
            kb_x = f"0,{current_position},{end_position}"
            star_x = f"0_{current_position}_0_{end_position}"
        else:
            kb_x = f"{kb_x},0,{current_position},{end_position}"
            star_x = f"{star_x} 0_{current_position}_0_{end_position}"
        current_position = end_position
    
    # UMI information if present
    star_umi = None
    if umi_range is not None:
        umi_length = umi_range[1] - umi_range[0]
        json_data["umi:"].append({
            "range": f"1:{current_position + 1}-{current_position + umi_length}"
        })
        kb_x = f"{kb_x}:0,{current_position},{current_position + umi_length}"
        star_umi = f"0_{current_position},0,{current_position + umi_length}"

    
    json_data["kallisto-bustools"].append({
        "kb count": f"-i </path/to/transcriptome.idx> -g </path/to/transcripts_to_genes> -x {kb_x}:1,0,0 -w NONE --h5ad --inleaved -o <outdir> {output_prefix}.fastq"
    })


    json_file = f"{output_prefix}.json"
    with open(json_file, "w") as f:
        json.dump(json_data, f, indent=4)

    return json_file