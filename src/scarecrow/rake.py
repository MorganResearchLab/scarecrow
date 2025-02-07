#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
"""

from argparse import RawTextHelpFormatter
import logging
import json
import itertools
from typing import List, Dict
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string


def parser_rake(parser):
    subparser = parser.add_parser(
        "rake",
        description="""
Generate barcode mismatch lists.

Example:

scarecrow rake --barcodes whitelist.txt --max_mismatch 3 --out barcode_mismatches.json
---
""",
        help="Rake barcodes to output mismatch lists",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument(
        "-b", "--barcodes",
        metavar="<file>",
        help=("Path to barcode text file"),
        type=str,
        required=True,
        default=None,
    )
    subparser.add_argument(
        "-o", "--out",
        metavar="<file>",
        help=("Path to output barcode mismatch file"),
        type=str,
        default=None
    )
    subparser.add_argument(
        "-m", "--max_mismatches",
        metavar="<int>",
        help=("Maximum number of mismatches in a barcode to characterise [1]"),
        type=int,
        default=1,
    )        
    return subparser

def validate_rake_args(parser, args):
    run_rake(barcodes = args.barcodes,
            output_file = args.out,
            max_mismatches = args.max_mismatches)
    
@log_errors
def run_rake(barcodes: str = None,
             output_file: str = None,
             max_mismatches: int = 1) -> None:
    """
    Function to rake barcodes for mismatches
    """
    # Setup logging
    logfile = f'./scarecrow_rake_{generate_random_string()}.log'
    logger = setup_logger(logfile)
    logger.info(f"logfile: '{logfile}'")
    
    if barcodes:
        logger.info(f"Processing barcodes in {barcodes}")
        if output_file is None:
            output_file = f'{barcodes}.json'
        process_barcodes(barcodes, 2, output_file)

    
def generate_mismatches(sequence: str, max_mismatches: int) -> Dict[int, Dict[str, str]]:
    bases = {'A', 'C', 'G', 'T'}
    mismatch_dict = {i: {} for i in range(1, max_mismatches + 1)}
    
    seq_length = len(sequence)
    for num_mismatches in range(1, max_mismatches + 1):
        for positions in itertools.combinations(range(seq_length), num_mismatches):
            for replacements in itertools.product(bases, repeat=num_mismatches):
                mutated_seq = list(sequence)
                for pos, new_base in zip(positions, replacements):
                    if mutated_seq[pos] != new_base:  # Ensure we introduce a mismatch
                        mutated_seq[pos] = new_base
                mutated_seq = "".join(mutated_seq)
                if mutated_seq not in mismatch_dict[num_mismatches]:
                    mismatch_dict[num_mismatches][mutated_seq] = sequence
    
    return mismatch_dict

def process_barcodes(input_file: str, max_mismatches: int, output_file: str):
    with open(input_file, 'r') as f:
        barcodes = [line.strip() for line in f if line.strip()]
    
    result = {0: {barcode: [barcode] for barcode in barcodes}}  # Mismatch=0 is the original
    
    for barcode in barcodes:
        mismatches = generate_mismatches(barcode, max_mismatches)
        for mismatch_level, seq_dict in mismatches.items():
            if mismatch_level not in result:
                result[mismatch_level] = {}
            result[mismatch_level].update(seq_dict)
    
    with open(output_file, 'w') as f:
        json.dump(result, f, indent=4)
    
    print(f"JSON file saved: {output_file}")    