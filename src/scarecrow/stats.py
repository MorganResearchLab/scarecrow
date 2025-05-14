# -*- coding: utf-8 -*-
"""
#!/usr/bin/env python3
@author: David Wragg
"""

import pysam
import gzip
import logging
from collections import defaultdict
from pathlib import Path
from typing import Tuple, Dict, DefaultDict
from argparse import RawTextHelpFormatter
from scarecrow import __version__
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string

def parser_stats(parser):
    subparser = parser.add_parser(
        "stats",
        description="""
Generate barcode metrics from SAM or FASTQ read tags/header, including position and mismatch counts.

Example:

scarecrow stats --in cDNA.sam
---
""",
        help="Generate barcode metrics from SAM or FASTQ read tags/header, including position and mismatch counts",
        formatter_class=RawTextHelpFormatter,
    )
    # Add a mutually exclusive group for output format
    subparser.add_argument(
        "-i",
        "--in",
        dest="infile",
        metavar="<file>",
        help=("SAM or FASTQ file to generate stats from"),
        type=str,
        required=True,
    )
    return subparser


def validate_stats_args(parser, args) -> None:
    """
    Validate arguments
    """
    # Global logger setup
    logfile = "{}_{}.{}".format(
        "./scarecrow_stats", generate_random_string(), "log"
    )
    logger = setup_logger(logfile)
    logger.info(f"scarecrow version {__version__}")
    logger.info(f"logfile: '{logfile}'")

    run_stats(in_file=args.infile)


@log_errors
def run_stats(in_file: str = None) -> None:
    """
    Main function to extract stats from SAM or FASTQ read tags
    """
    logger = logging.getLogger("scarecrow")

    # List of accepted file types
    suffixes = ['.fastq.gz', '.fq.gz', '.fastq', '.fq', '.sam']
    suffix = next((s for s in suffixes if str(in_file).endswith(s)), None)

    # Get count data from SAM file
    names = [
        "CR_counts",
        "CB_counts",
        "XP_counts",
        "XM_counts",
        "UR_counts",
        "combined_CR_counts",
        "combined_CB_counts",
    ]
    
    # Validate file exists
    if Path(in_file).exists():
        
        if suffix == '.sam':            
            counts_list = list(parse_sam_tags(in_file))

        elif suffix in ['.fastq', '.fq', '.fastq.gz', '.fq.gz']:
            counts_list = list(parse_fastq_tags(in_file))        

        if counts_list:
            # Iterate through results and write each to file
            for name, counts in zip(names, counts_list):
                try:
                    out_file = f"{in_file}.{name}.txt"
                    if "combined" in name:
                        write_combined_counts_to_file(counts, out_file)
                    else:
                        write_counts_to_file(counts, out_file)
                    logger.info(f"{name} written to '{out_file}'")
                except ValueError:
                    logger.info(f"Encountered a problem with {name}")
        else:
            raise ValueError(f"Unsupported file type for {in_file}")


    logger.info("Finished!")


@log_errors
def parse_fastq_tags(fastq_file: str = None) -> Tuple[
    DefaultDict[int, DefaultDict[str, int]],  # CR_counts
    DefaultDict[int, DefaultDict[str, int]],  # CB_counts
    DefaultDict[int, DefaultDict[str, int]],  # XP_counts
    DefaultDict[int, DefaultDict[str, int]],  # XM_counts
    DefaultDict[int, DefaultDict[str, int]],  # UR_counts
    DefaultDict[str, int],                   # combined_CR_counts
    DefaultDict[str, int]                    # combined_CB_counts
]:
    """
    Read FASTQ file and extract read TAG data from header of read 1
    """
    logger = logging.getLogger("scarecrow")
    logger.info(f"Processing read tags in header of read 1 from '{fastq_file}'")

    # Initialize dictionaries to store counts
    CR_counts = defaultdict(lambda: defaultdict(int))
    CB_counts = defaultdict(lambda: defaultdict(int))
    XP_counts = defaultdict(lambda: defaultdict(int))
    XM_counts = defaultdict(lambda: defaultdict(int))
    UR_counts = defaultdict(lambda: defaultdict(int))
    combined_CR_counts = defaultdict(int)
    combined_CB_counts = defaultdict(int)

    # Handle gzipped or plain fastq
    open_func = gzip.open if fastq_file.endswith('.gz') else open
    
    with open_func(fastq_file, 'rt') as f:
        while True:
            # Read R1 header
            header = f.readline().strip()
            if not header:
                break
            
            # Skip R1 sequence/quality and R2 entirely
            for _ in range(7):
                f.readline()

            # Extract tag values with proper underscore handling
            def get_tag_parts(tag_name):
                tag_pos = header.find(tag_name + ':')
                if tag_pos == -1:
                    return None
                
                value_start = tag_pos + len(tag_name) + 1
                value_end = len(header)
                
                # Find next tag or end of string
                for end_marker in [':', '/']:
                    end_pos = header.find(end_marker, value_start)
                    if end_pos != -1:
                        value_end = min(value_end, end_pos)
                
                value = header[value_start:value_end]
                return value.split('_') if value else None

            # Process CR tag
            cr_parts = get_tag_parts('CR')
            if cr_parts:
                for idx, barcode in enumerate(cr_parts):
                    CR_counts[idx][barcode] += 1
                combined_CR_counts["_".join(cr_parts)] += 1  # Store with commas to match SAM

            # Process CB tag
            cb_parts = get_tag_parts('CB')
            if cb_parts:
                for idx, barcode in enumerate(cb_parts):
                    CB_counts[idx][barcode] += 1
                combined_CB_counts["_".join(cb_parts)] += 1  # Store with commas to match SAM

            # Process XP tag
            xp_parts = get_tag_parts('XP')
            if xp_parts:
                for idx, pos in enumerate(xp_parts):
                    XP_counts[idx][pos] += 1

            # Process XM tag
            xm_parts = get_tag_parts('XM')
            if xm_parts:
                for idx, mm in enumerate(xm_parts):
                    XM_counts[idx][mm] += 1

            # Process UR tag (single value)
            ur_val = get_tag_parts('UR')
            if ur_val:
                UR_counts[0][ur_val[0]] += 1

    return (
        CR_counts,
        CB_counts,
        XP_counts,
        XM_counts,
        UR_counts,
        combined_CR_counts,
        combined_CB_counts,
    )





@log_errors
def parse_sam_tags(sam_file: str = None) -> Tuple[
    DefaultDict[int, DefaultDict[str, int]],  # CR_counts
    DefaultDict[int, DefaultDict[str, int]],  # CB_counts
    DefaultDict[int, DefaultDict[str, int]],  # XP_counts
    DefaultDict[int, DefaultDict[str, int]],  # XM_counts
    DefaultDict[int, DefaultDict[str, int]],  # UR_counts
    DefaultDict[str, int],                   # combined_CR_counts
    DefaultDict[str, int]                    # combined_CB_counts
]:
    """
    Read SAM file and extract read TAG data
    """
    logger = logging.getLogger("scarecrow")
    logger.info(f"Processing read tags from '{sam_file}'")

    # Initialize dictionaries to store counts
    CR_counts = defaultdict(lambda: defaultdict(int))
    CB_counts = defaultdict(lambda: defaultdict(int))
    XP_counts = defaultdict(lambda: defaultdict(int))
    XM_counts = defaultdict(lambda: defaultdict(int))
    UR_counts = defaultdict(lambda: defaultdict(int))
    combined_CR_counts = defaultdict(int)
    combined_CB_counts = defaultdict(int)

    # Open SAM file
    with pysam.AlignmentFile(sam_file, "rb", check_sq=False) as sam:
        # Force reading even if header is missing
        for read in sam.fetch(until_eof=True):
            # Extract tags and incorporate into header
            tags = {tag: value for tag, value in read.tags}

            # Process CR tag
            if "CR" in tags:
                CR_barcodes = tags["CR"].split(",")
                for idx, barcode in enumerate(CR_barcodes):
                    CR_counts[idx][barcode] += 1
                combined_CR_counts["_".join(CR_barcodes)] += 1

            # Process CB tag
            if "CB" in tags:
                CB_barcodes = tags["CB"].split(",")
                for idx, barcode in enumerate(CB_barcodes):
                    CB_counts[idx][barcode] += 1
                combined_CB_counts["_".join(CB_barcodes)] += 1

            # Process XP tag
            if "XP" in tags:
                XP_barcodes = tags["XP"].split(",")
                for idx, barcode in enumerate(XP_barcodes):
                    XP_counts[idx][barcode] += 1

            # Process XM tag
            if "XM" in tags:
                XM_barcodes = tags["XM"].split(",")
                for idx, barcode in enumerate(XM_barcodes):
                    XM_counts[idx][barcode] += 1

            # Process UR tag
            if "UR" in tags:
                UR_barcodes = tags["UR"].split(",")
                for idx, barcode in enumerate(UR_barcodes):
                    UR_counts[idx][barcode] += 1

        return (
            CR_counts,
            CB_counts,
            XP_counts,
            XM_counts,
            UR_counts,
            combined_CR_counts,
            combined_CB_counts,
        )


def write_counts_to_file(counts, filename):
    with open(filename, "w") as file:
        for idx, barcode_dict in sorted(counts.items()):
            for barcode, count in sorted(barcode_dict.items()):
                file.write(f"Index {idx}\t{barcode}\t{count}\n")


def write_combined_counts_to_file(counts, filename):
    with open(filename, "w") as file:
        for barcode, count in sorted(counts.items()):
            file.write(f"{barcode}\t{count}\n")
