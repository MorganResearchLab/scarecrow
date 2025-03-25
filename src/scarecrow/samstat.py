# -*- coding: utf-8 -*-
"""
#!/usr/bin/env python3
@author: David Wragg
"""

import pysam
import logging
from collections import defaultdict
from pathlib import Path
from argparse import RawTextHelpFormatter
from scarecrow import __version__
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string


def parser_samstat(parser):
    subparser = parser.add_parser(
        "samstat",
        description="""
Generate barcode metrics from SAM file including position and mismatch counts.

Example:

scarecrow samstat --sam cDNA.sam
---
""",
        help="Generate barcode metrics from SAM file including position and mismatch counts",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument(
        "-s",
        "--sam",
        metavar="<file>",
        help=("SAM file to process"),
        type=str,
        required=True,
        default=[],
    )
    return subparser


def validate_samstat_args(parser, args) -> None:
    """
    Validate arguments
    """
    # Global logger setup
    logfile = "{}_{}.{}".format(
        "./scarecrow_sam2fastq", generate_random_string(), "log"
    )
    logger = setup_logger(logfile)
    logger.info(f"scarecrow version {__version__}")
    logger.info(f"logfile: '{logfile}'")

    run_samstat(sam_file=args.sam)


@log_errors
def run_samstat(sam_file: str = None) -> None:
    """
    Main function to extract sequences with barcode headers
    """
    logger = logging.getLogger("scarecrow")

    # Validate file exists
    if sam_file:
        if Path(sam_file).exists():
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
            counts_list = list(parse_sam_tags(sam_file))

            # Iterate through results and write each to file
            for name, counts in zip(names, counts_list):
                try:
                    out_file = f"{sam_file}.{name}.txt"
                    if "combined" in name:
                        write_combined_counts_to_file(counts, out_file)
                    else:
                        write_counts_to_file(counts, out_file)
                    logger.info(f"{name} written to '{out_file}'")
                except ValueError:
                    logger.info(f"Encountered a problem with {name}")
        else:
            logger.info(f"'{sam_file}' does not exist")
    else:
        logger.info("No SAM file provided")

    logger.info("Finished!")


@log_errors
def parse_sam_tags(sam_file: str = None):
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
                CR_barcodes = tags["CR"].split("_")
                for idx, barcode in enumerate(CR_barcodes):
                    CR_counts[idx][barcode] += 1
                combined_CR_counts["_".join(CR_barcodes)] += 1

            # Process CB tag
            if "CB" in tags:
                CB_barcodes = tags["CB"].split("_")
                for idx, barcode in enumerate(CB_barcodes):
                    CB_counts[idx][barcode] += 1
                combined_CB_counts["_".join(CB_barcodes)] += 1

            # Process XP tag
            if "XP" in tags:
                XP_barcodes = tags["XP"].split("_")
                for idx, barcode in enumerate(XP_barcodes):
                    XP_counts[idx][barcode] += 1

            # Process XM tag
            if "XM" in tags:
                XM_barcodes = tags["XM"].split("_")
                for idx, barcode in enumerate(XM_barcodes):
                    XM_counts[idx][barcode] += 1

            # Process UR tag
            if "UR" in tags:
                UR_barcodes = tags["UR"].split("_")
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
