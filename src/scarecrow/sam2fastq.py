# -*- coding: utf-8 -*-
"""
#!/usr/bin/env python3
@author: David Wragg
"""

import pysam
import logging
from pathlib import Path
from argparse import RawTextHelpFormatter
from scarecrow import __version__
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string

def parser_sam2fastq(parser):
    subparser = parser.add_parser(
        "sam2fastq",
        description="""
Tally barcode and barcode combination counts in fastq output of scarecrow reap.

Example:

scarecrow sam2fastq --sam cdna.sam
---
""",
        help="Extracts reads from SAM and writes to FASTQ retaining read tags in sequence header.",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument(
        "-s", "--sam",
        metavar="<file>",
        help=("SAM file to extract reads from"),
        type=str,
        required=True,
        default=[],
    )
    return subparser

def validate_sam2fastq_args(parser, args) -> None:
    """ 
    Validate arguments 
    """
    # Global logger setup
    logfile = '{}_{}.{}'.format('./scarecrow_sam2fastq', generate_random_string(), 'log')
    logger = setup_logger(logfile)
    logger.info(f"scarecrow version {__version__}")
    logger.info(f"logfile: '{logfile}'")

    run_sam2fastq(sam_file = args.sam)

@log_errors
def run_sam2fastq(sam_file: str = None) -> None:
    """
    Main function to extract sequences with barcode headers
    """    
    logger = logging.getLogger('scarecrow')

    # Validate file exists
    if sam_file:
        if Path(sam_file).exists():

            fastq_file = f'{sam_file}.fastq'
            logger.info(f"Extracting reads from '{sam_file}' to '{fastq_file}'")
        
            with pysam.AlignmentFile(sam_file, "rb", check_sq=False) as sam, open(fastq_file, "w") as fq:
                # Force reading even if header is missing
                for read in sam.fetch(until_eof=True):
                    # Extract tags and incorporate into header
                    tags = {k: str(v) for k, v in read.tags}
                    tag_string = "|".join([f"{k}:{v}" for k, v in tags.items()])
                    header = f"@{read.query_name}|{tag_string}"
                    # Convert Phred scores correctly (add 33 for ASCII)
                    quality_string = "".join(chr(q + 33) for q in read.query_qualities)
                    fq.write(f"{header}\n{read.query_sequence}\n+\n{quality_string}\n")
    
        else:
            logger.info(f"'{sam_file}' does not exist")
    else:
        logger.info("No SAM file provided")