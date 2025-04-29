# -*- coding: utf-8 -*-
"""
#!/usr/bin/env python3
@author: David Wragg
"""

import gzip
import json
import pysam
import logging
from collections import defaultdict
from pathlib import Path
from argparse import RawTextHelpFormatter
from scarecrow import __version__
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string

def parser_sift(parser):
    subparser = parser.add_parser(
        "sift",
        description="""
Sift a SAM file to remove reads with invalid barcodes.

Example:

scarecrow sift --in cdna.sam
---
""",
        help="Removes reads with invalid barcodes from SAM file.",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument(
        "-i",
        "--in",
        dest="infile",
        metavar="<file>",
        help=("Interleaved FASTQ file or SAM file to sift for valid barcodes"),
        type=str,
        required=True,
        default=None,
    )
    subparser.add_argument(
        "-j",
        "--json",
        metavar="<file>",
        help=("JSON file to accompany interleaved FASTQ file, to identify barcode ranges"),
        type=str,
        default=None,
    )    
    return subparser


def validate_sift_args(parser, args) -> None:
    """
    Validate arguments
    """
    # Global logger setup
    logfile = "{}_{}.{}".format(
        "./scarecrow_sift", generate_random_string(), "log"
    )
    logger = setup_logger(logfile)
    logger.info(f"scarecrow version {__version__}")
    logger.info(f"logfile: '{logfile}'")

    run_sift(input_file=args.infile, json_file=args.json)


@log_errors
def run_sift(input_file: str = None, json_file: str = None) -> None:
    """
    Main function to sift reads to remove those with invalid barcodes
    """
    logger = logging.getLogger("scarecrow")

    # Validate file exists
    if Path(input_file).exists():
        if input_file.endswith(".sam"):
            sift_sam(input_file)
        elif input_file.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
            if Path(json_file).exists():
                if json_file.endswith((".json")):
                    sift_fastq(input_file, json_file)
                else:
                    logger.info("Passed JSON file suffix is not .json")
                    raise IOError
            else:
                raise FileNotFoundError
        else:
            logger.info("Input file suffix is not one of : '.sam', '.fastq', '.fastq.gz', '.fq', '.fq.gz'")
            raise IOError
    else:
        raise FileNotFoundError

    logger.info("Finished!")

@log_errors
def sift_sam(sam_file: str = None):
    """
    Read SAM file and sift invalid barcodes
    """
    logger = logging.getLogger("scarecrow")
    output_sam = sam_file.replace(".sam", "_sift.sam")
    logger.info(f"Sifting '{sam_file}', results will be written to '{output_sam}")
    with pysam.AlignmentFile(sam_file, "r", check_sq=False) as sam, pysam.AlignmentFile(output_sam, "w", template=sam) as outfile:
        # Force reading even if header is missing
        for read in sam.fetch(until_eof=True):
            cb_tag = read.get_tag("CB") if read.has_tag("CB") else ""
            if "N" not in cb_tag:
                outfile.write(read)

@log_errors
def sift_fastq(fastq_file: str = None, json_file: str = None):
    """
    Read interleaved FASTQ file and sift invalid barcodes
    """
    logger = logging.getLogger("scarecrow")

    # Load JSON file
    with open(json_file) as f:
        config = json.load(f)
    
    # Extract barcode ranges
    barcode_ranges = []
    for bc in config['barcodes']:
        barcode_ranges.append(parse_range(bc['range']))
    logger.info(f"barcode ranges: {barcode_ranges}")

    # Open input and output files
    output_fastq = fastq_file.replace(".fastq", "_sift.fastq")
    output_fastq = output_fastq[:-3] if output_fastq.endswith('.gz') else output_fastq
    open_func = gzip.open if fastq_file.endswith('.gz') else open
    logger.info(f"Sifting '{fastq_file}', results will be written to '{output_fastq}' and '{json_file}'")            
    with open_func(fastq_file, 'rt') as infile, open(output_fastq, 'wt') as outfile:
        while True:
            # Read 4 lines (one record)
            lines = [infile.readline() for _ in range(8)]
            if not lines[0]:  # end of file
                break
            
            # Extract read 1 sequence (2nd line of the first 4-line block)
            seq = lines[1].strip()
            #logger.info(f"{seq}")

            # Check if any barcode range contains 'N'
            keep = True
            for read_num, start, end in barcode_ranges:
                # Get the sequence line (2nd line of the record)
                barcode = seq[start:end+1]  # end is inclusive                        
                #logger.info(f"{barcode}")
                if 'N' in barcode:
                    keep = False
                    break
            
            if keep:
                outfile.writelines(lines)
    
def parse_range(range_str):
    """Parse a range string like '1:0-7' into (read_num, start, end)."""
    read_part, coords = range_str.split(':')
    read_num = int(read_part) - 1  # convert to 0-based index
    start, end = map(int, coords.split('-'))
    return (read_num, start, end)