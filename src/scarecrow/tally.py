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
        metavar="fastq",
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

    # Log the barcode counts for each position
    for i, counts in enumerate(barcode_counts):
        if verbose:
            for barcode, count in counts.items():
                logger.info(f"Barcode index: {i + 1}\tBarcode: {barcode}\tCount: {count}")
        barcodes = pd.DataFrame(list(barcode_counts[i].items()), columns=["Barcode", "Count"]).sort_values(by="Count")
        barcodes.insert(0, "Index", i + 1)
        if i == 0:
            barcodes.to_csv('{}.{}'.format(fastq, 'barcodes.csv'), index = False)
        else:
            barcodes.to_csv('{}.{}'.format(fastq, 'barcodes.csv'), index = False, mode = "a", header = False)

    # Log the combined barcode counts (i.e. cell sequence counts)
    if verbose:
        for cell, count in cell_barcodes.items():
            logger.info(f"Barcode combination: {cell}\tCount: {count}")
    barcodes = pd.DataFrame(list(cell_barcodes.items()), columns=["BarcodeCombination", "Count"]).sort_values(by="Count")
    barcodes.to_csv('{}.{}'.format(fastq, 'barcode.combinations.csv'), index = False)


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
