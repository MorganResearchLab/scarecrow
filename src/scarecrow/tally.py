# -*- coding: utf-8 -*-
"""
#!/usr/bin/env python3
@author: David Wragg
"""

import logging
import pandas as pd
import pysam
import re
from argparse import RawTextHelpFormatter
from collections import defaultdict
from scarecrow import __version__
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string
from typing import List, Tuple


def parser_tally(parser):
    subparser = parser.add_parser(
        "tally",
        description="""
Tally barcode and barcode combination counts in fastq output of scarecrow reap.

Example:

scarecrow tally --fastq cdna.fastq.gz
---
""",
        help="Tally barcode and barcode combination counts in fastq output of scarecrow reap.",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument(
        "-i",
        "--input",
        metavar="<file>",
        help=("In FASTQ or SAM file as output by scarecrow reap"),
        type=str,
        required=True,
        default=[],
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
        "-v", "--verbose", action="store_true", help="Enable verbose output [false]"
    )
    return subparser


def validate_tally_args(parser, args) -> None:
    """
    Validate arguments
    """
    # Global logger setup
    logfile = "{}_{}.{}".format("./scarecrow_tally", generate_random_string(), "log")
    logger = setup_logger(logfile)
    logger.info(f"scarecrow version {__version__}")
    logger.info(f"logfile: '{logfile}'")

    run_tally(input_file=args.input, mismatches=args.mismatches, verbose=args.verbose)


@log_errors
def run_tally(
    input_file: str = None, mismatches: int = 1, verbose: bool = False
) -> None:
    """
    Main function to extract sequences with barcode headers
    """
    logger = logging.getLogger("scarecrow")

    # Determine file type
    is_fastq = input_file.lower().endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz"))

    # Process file headers
    logger.info("Processing sequence headers, this may take a while")
    barcode_counts, cell_barcodes = (
        process_fastq_headers(input_file)
        if is_fastq
        else process_sam_headers(input_file)
    )

    # Report barcode metrics summary
    mismatch_counts, valid_barcodes, sequence_counts, perfect_matches = (
        count_fastq_mismatch_sums(input_file)
        if is_fastq
        else count_sam_mismatch_sums(input_file)
    )

    if verbose:
        for mismatch_sum, count in sorted(mismatch_counts.items()):
            logger.info(f"Mismatch sum: {mismatch_sum}\tCount: {count}")
    logger.info(f"Total sequences processed: {sequence_counts}")
    logger.info(f"Sequences with non-null barcodes: {valid_barcodes}")
    logger.info(
        f"Barcodes with <= {mismatches} mismatches: {sum(mismatch_counts.values())}"
    )
    logger.info(f"Sequences with perfect matching barcodes: {perfect_matches}")

    # Save mismatch counts to CSV
    mismatch_df = pd.DataFrame(
        list(mismatch_counts.items()), columns=["BarcodeMismatches", "Count"]
    ).sort_values(by="BarcodeMismatches")
    mismatch_df.to_csv(f"{input_file}.mismatches.csv", index=False)

    # Count position occurrences
    position_counts = (
        count_fastq_position_occurrences(input_file)
        if is_fastq
        else count_sam_position_occurrences(input_file)
    )

    if verbose:
        for i, counts in enumerate(position_counts):
            for position, count in sorted(counts.items()):
                logger.info(
                    f"Position index: {i + 1}\tPosition: {position}\tCount: {count}"
                )

    # Save position counts to CSV
    for i, counts in enumerate(position_counts):
        positions_df = pd.DataFrame(
            list(counts.items()), columns=["Position", "Count"]
        ).sort_values(by="Position")
        positions_df.insert(0, "Index", i + 1)
        if i == 0:
            positions_df.to_csv(f"{input_file}.positions.csv", index=False)
        else:
            positions_df.to_csv(
                f"{input_file}.positions.csv", index=False, mode="a", header=False
            )

    # Original barcode counting logic
    for i, counts in enumerate(barcode_counts):
        if verbose:
            for barcode, count in counts.items():
                logger.info(
                    f"Barcode index: {i + 1}\tBarcode: {barcode}\tCount: {count}"
                )
        barcodes = pd.DataFrame(
            list(barcode_counts[i].items()), columns=["Barcode", "Count"]
        ).sort_values(by="Count")
        barcodes.insert(0, "Index", i + 1)
        if i == 0:
            barcodes.to_csv(f"{input_file}.barcodes.csv", index=False)
        else:
            barcodes.to_csv(
                f"{input_file}.barcodes.csv", index=False, mode="a", header=False
            )

    # Log the combined barcode counts
    if verbose:
        for cell, count in cell_barcodes.items():
            logger.info(f"Barcode combination: {cell}\tCount: {count}")
    barcodes = pd.DataFrame(
        list(cell_barcodes.items()), columns=["BarcodeCombination", "Count"]
    ).sort_values(by="Count")
    barcodes.to_csv(f"{input_file}.barcode.combinations.csv", index=False)

    logger.info("Finished!")


def process_fastq_headers(
    file_path: str,
) -> Tuple[List[defaultdict[str, int]], defaultdict[str, int]]:
    """
    Process fastq header to report on barcode counts
    """
    barcode_counts = []
    cell_barcodes = defaultdict(int)

    with pysam.FastxFile(file_path) as fastq_file:
        for entry in fastq_file:
            header = entry.comment
            if "CB=" in header:
                barcodes_str = re.search(r"CB=([\w_]+)", header).group(1)
                cell_barcodes[barcodes_str] += 1

                barcodes = barcodes_str.split("_")

                while len(barcode_counts) < len(barcodes):
                    barcode_counts.append(defaultdict(int))

                for i, barcode in enumerate(barcodes):
                    barcode_counts[i][barcode] += 1

    return barcode_counts, cell_barcodes


def process_sam_headers(
    file_path: str,
) -> Tuple[List[defaultdict[str, int]], defaultdict[str, int]]:
    """
    Process SAM header to report on barcode counts
    """
    barcode_counts = []
    cell_barcodes = defaultdict(int)

    with pysam.AlignmentFile(file_path, "rb", check_sq=False) as sam_file:
        # Force reading even if header is missing
        for read in sam_file.fetch(until_eof=True):
            # Check for cell barcode tag
            if read.has_tag("CB"):
                barcodes_str = read.get_tag("CB")
                cell_barcodes[barcodes_str] += 1

                barcodes = barcodes_str.split("_")

                while len(barcode_counts) < len(barcodes):
                    barcode_counts.append(defaultdict(int))

                for i, barcode in enumerate(barcodes):
                    barcode_counts[i][barcode] += 1

    return barcode_counts, cell_barcodes


def count_fastq_mismatch_sums(file_path: str) -> Tuple[dict, int, int, int]:
    """
    Count sequences for each sum of mismatch values in FASTQ
    """
    mismatch_counts = defaultdict(int)
    sequences = 0
    valid_barcodes = 0
    perfect_matches = 0

    with pysam.FastxFile(file_path) as fastq_file:
        for entry in fastq_file:
            header = entry.comment
            if header:
                sequences += 1
                if "null" not in header:
                    valid_barcodes += 1
                    if "XM=" in header:
                        mismatches_str = re.search(r"XM=([\d_]+)", header).group(1)
                        mismatch_sum = 0
                        for x in mismatches_str.split("_"):
                            mismatch_sum += int(x)
                            mismatch_counts[int(x)] += 1
                        if mismatch_sum == 0:
                            perfect_matches += 1

    return dict(mismatch_counts), valid_barcodes, sequences, perfect_matches


def count_sam_mismatch_sums(file_path: str) -> Tuple[dict, int, int, int]:
    """
    Count sequences for each sum of mismatch values in SAM
    """
    mismatch_counts = defaultdict(int)
    sequences = 0
    valid_barcodes = 0
    perfect_matches = 0

    with pysam.AlignmentFile(file_path, "r", check_sq=False) as sam_file:
        # Force reading even if header is missing
        for read in sam_file.fetch(until_eof=True):
            sequences += 1

            # Check for cell barcode
            if read.has_tag("CB"):
                valid_barcodes += 1

                # Check for mismatch tag
                if read.has_tag("XM"):
                    mismatches_str = read.get_tag("XM")
                    mismatch_sum = 0
                    for x in mismatches_str.split("_"):
                        mismatch_val = int(x)
                        mismatch_sum += mismatch_val
                        mismatch_counts[mismatch_val] += 1

                    if mismatch_sum == 0:
                        perfect_matches += 1

    return dict(mismatch_counts), valid_barcodes, sequences, perfect_matches


def count_fastq_position_occurrences(file_path: str) -> List[dict]:
    """
    Count sequences sharing each reported position in FASTQ
    """
    position_counts = []

    with pysam.FastxFile(file_path) as fastq_file:
        for entry in fastq_file:
            header = entry.comment
            if "XP=" in header:
                positions_str = re.search(r"XP=([\d_]+)", header).group(1)
                positions = positions_str.split("_")

                while len(position_counts) < len(positions):
                    position_counts.append(defaultdict(int))

                for i, pos in enumerate(positions):
                    position_counts[i][int(pos)] += 1

    return [dict(counts) for counts in position_counts]


def count_sam_position_occurrences(file_path: str) -> List[dict]:
    """
    Count sequences sharing each reported position in SAM
    """
    position_counts = []

    with pysam.AlignmentFile(file_path, "r", check_sq=False) as sam_file:
        # Force reading even if header is missing
        for read in sam_file.fetch(until_eof=True):
            if read.has_tag("XP"):
                positions_str = read.get_tag("XP")
                positions = positions_str.split("_")

                while len(position_counts) < len(positions):
                    position_counts.append(defaultdict(int))

                for i, pos in enumerate(positions):
                    position_counts[i][int(pos)] += 1

    return [dict(counts) for counts in position_counts]
