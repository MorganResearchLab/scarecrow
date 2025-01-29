# -*- coding: utf-8 -*-
"""
#!/usr/bin/env python3
@author: David Wragg
"""

import pysam
import re
import pandas as pd
from argparse import RawTextHelpFormatter
from typing import List, Tuple
from collections import defaultdict
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string

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
        "-f", "--fastq",
        metavar="<file>",
        help=("Fastq file as output by scarecrow reap"),
        type=str,
        default=[],
    )    
    subparser.add_argument(
        "-v", "--verbose",
        action='store_true',
        help='Enable verbose output [false]'
    )
    return subparser

def validate_tally_args(parser, args) -> None:
    """ 
    Validate arguments 
    """
    run_tally(fastq = args.fastq,
              verbose = args.verbose)

@log_errors
def run_tally(fastq: str = None,
              verbose: bool = False) -> None:
    """
    Main function to extract sequences with barcode headers
    """    
    # Global logger setup
    logfile = '{}_{}.{}'.format('./scarecrow_tally', generate_random_string(), 'log')
    logger = setup_logger(logfile)
    logger.info(f"logfile: '{logfile}'")

    # Process fastq header
    barcode_counts, cell_barcodes = process_fastq_headers(fastq)

    # Count complete barcodes
    # Count complete barcodes and total sequences
    total_sequences, complete_barcode_count = count_complete_barcodes(fastq)
    logger.info(f"Total sequences processed: {total_sequences}")
    logger.info(f"Sequences with complete barcodes: {complete_barcode_count} ({(complete_barcode_count/total_sequences*100):.2f}%)")

    # Count mismatch sums
    mismatch_counts = count_mismatch_sums(fastq)
    if verbose:
        for mismatch_sum, count in sorted(mismatch_counts.items()):
            logger.info(f"Mismatch sum: {mismatch_sum}\tCount: {count}")
    
    # Save mismatch counts to CSV
    mismatch_df = pd.DataFrame(list(mismatch_counts.items()), 
                              columns=["MismatchSum", "Count"]).sort_values(by="MismatchSum")
    mismatch_df.to_csv(f"{fastq}.mismatches.csv", index=False)
    
    # Count position occurrences
    position_counts = count_position_occurrences(fastq)
    if verbose:
        for i, counts in enumerate(position_counts):
            for position, count in sorted(counts.items()):
                logger.info(f"Position index: {i + 1}\tPosition: {position}\tCount: {count}")
    
    # Save position counts to CSV
    for i, counts in enumerate(position_counts):
        positions_df = pd.DataFrame(list(counts.items()), 
                                  columns=["Position", "Count"]).sort_values(by="Position")
        positions_df.insert(0, "Index", i + 1)
        if i == 0:
            positions_df.to_csv(f"{fastq}.positions.csv", index=False)
        else:
            positions_df.to_csv(f"{fastq}.positions.csv", index=False, mode="a", header=False)

    # Original barcode counting logic
    for i, counts in enumerate(barcode_counts):
        if verbose:
            for barcode, count in counts.items():
                logger.info(f"Barcode index: {i + 1}\tBarcode: {barcode}\tCount: {count}")
        barcodes = pd.DataFrame(list(barcode_counts[i].items()), 
                              columns=["Barcode", "Count"]).sort_values(by="Count")
        barcodes.insert(0, "Index", i + 1)
        if i == 0:
            barcodes.to_csv(f"{fastq}.barcodes.csv", index=False)
        else:
            barcodes.to_csv(f"{fastq}.barcodes.csv", index=False, mode="a", header=False)

    # Log the combined barcode counts
    if verbose:
        for cell, count in cell_barcodes.items():
            logger.info(f"Barcode combination: {cell}\tCount: {count}")
    barcodes = pd.DataFrame(list(cell_barcodes.items()), 
                           columns=["BarcodeCombination", "Count"]).sort_values(by="Count")
    barcodes.to_csv(f"{fastq}.barcode.combinations.csv", index=False)


def process_fastq_headers(file_path: str = None) -> Tuple[List[defaultdict[str, int]], defaultdict[str, int]]:
    """
    Process fastq header to report on barcode counts

    Args:
        file_path: fastq file to operate on
    
    Returns:
        (1) list of dictionary of barcode counts, and (2) dictionary of barcode combination counts
    """

    # Create a list of dictionaries, one for each barcode position
    barcode_counts = []
    cell_barcodes = defaultdict(int)

    # Open the FASTQ file using pysam
    with pysam.FastxFile(file_path) as fastq_file:
        for entry in fastq_file:
            # Extract the header line
            header = entry.comment

            # Check if the header contains 'barcodes='
            if 'barcodes=' in header:
                # Extract the barcodes substring
                barcodes_str = re.search(r'barcodes=([\w_]+)', header).group(1)
                cell_barcodes[barcodes_str] += 1
                
                # Split the barcodes string by underscore
                barcodes = barcodes_str.split('_')
                
                # Ensure barcode_counts has enough dictionaries for all barcode positions
                while len(barcode_counts) < len(barcodes):
                    barcode_counts.append(defaultdict(int))
                
                # Update counts in the corresponding dictionaries
                for i, barcode in enumerate(barcodes):
                    barcode_counts[i][barcode] += 1

    return barcode_counts, cell_barcodes

def count_complete_barcodes(file_path: str) -> Tuple[int, int]:
    """
    Count total sequences and sequences that have no 'null' barcodes.
    
    Args:
        file_path: Path to the fastq file
    
    Returns:
        Tuple of (total sequence count, count of sequences with complete non-null barcodes)
    """
    complete_count = 0
    total_count = 0
    
    with pysam.FastxFile(file_path) as fastq_file:
        for entry in fastq_file:
            total_count += 1
            header = entry.comment
            if 'barcodes=' in header:
                barcodes_str = re.search(r'barcodes=([\w_]+)', header).group(1)
                # Check if any barcode is 'null'
                if 'null' not in barcodes_str.lower():
                    complete_count += 1

    return total_count, complete_count

def count_mismatch_sums(file_path: str) -> dict:
    """
    Count sequences for each sum of mismatch values.
    
    Args:
        file_path: Path to the fastq file
    
    Returns:
        Dictionary with mismatch sums as keys and sequence counts as values
    """
    mismatch_counts = defaultdict(int)
    
    with pysam.FastxFile(file_path) as fastq_file:
        for entry in fastq_file:
            header = entry.comment
            if 'mismatches=' in header:
                mismatches_str = re.search(r'mismatches=([\d_]+)', header).group(1)
                # Convert string of mismatches to sum
                mismatch_sum = sum(int(x) for x in mismatches_str.split('_'))
                mismatch_counts[mismatch_sum] += 1
                
    return dict(mismatch_counts)

def count_position_occurrences(file_path: str) -> List[dict]:
    """
    Count sequences sharing each reported position for each index.
    
    Args:
        file_path: Path to the fastq file
    
    Returns:
        List of dictionaries, one for each position index, containing position counts
    """
    position_counts = []
    
    with pysam.FastxFile(file_path) as fastq_file:
        for entry in fastq_file:
            header = entry.comment
            if 'positions=' in header:
                positions_str = re.search(r'positions=([\d_]+)', header).group(1)
                positions = positions_str.split('_')
                
                # Ensure we have enough dictionaries for all positions
                while len(position_counts) < len(positions):
                    position_counts.append(defaultdict(int))
                
                # Count occurrences of each position at each index
                for i, pos in enumerate(positions):
                    position_counts[i][int(pos)] += 1
                    
    return [dict(counts) for counts in position_counts]