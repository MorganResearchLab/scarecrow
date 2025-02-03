#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
"""

import logging
import pysam
import os
from argparse import RawTextHelpFormatter
from typing import List, Dict
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string


def parser_samtag(parser):
    subparser = parser.add_parser(
        "samtag",
        description="""
Update SAM file with barcode and UMI tags from scarecrow reap fastq header.

Example:

scarecrow samtag --fastq in.fastq --sam in.sam
---
""",
        help="Update SAM file with barcode and UMI tags from scarecrow reap fastq header",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument(
        "-f", "--fastq",
        metavar="<file>",
        help=("Path to scarecrow reap fastq file"),
        type=str,
        default=None,
    )
    subparser.add_argument(
        "-s", "--sam",
        metavar="<file>",
        help=("Path to SAM file to update tags"),
        type=str,
        default=None,
    )
    return subparser

def validate_samtag_args(parser, args):
    run_samtag(fastq_file = args.fastq,
               sam_file = args.sam)
    
@log_errors
def run_samtag(fastq_file: str = None,
               sam_file: str = None) -> None:
    """
    Function to rake barcodes for mismatches
    """
    # Setup logging
    logfile = f'./scarecrow_samtag_{generate_random_string()}.log'
    logger = setup_logger(logfile)
    logger.info(f"logfile: '{logfile}'")
    
    # Preprocess fastq headers
    if fastq_file:
        logger.info(f"Pre-processing {fastq_file}")
        read_tags = preprocess_fastq_headers(fastq_file)    

    # Add tags to SAM file
    if sam_file:
        name, _ = os.path.splitext(fastq_file)
        output = name + ".tagged.sam"
        logger.info(f"Adding tags and writing to {output}")
        add_tags_to_sam(sam_file, output, read_tags)

def preprocess_fastq_headers(fastq_file):
    """Preprocess FASTQ headers to create a mapping of read names to tags."""
    read_tags = {}
    with open(fastq_file, "r") as fq:
        for line in fq:
            if line.startswith("@"):
                # Extract the read name and tags
                fastq_header = line.strip()
                read_name = fastq_header.split(" ")[0][1:]  # Remove "@" and split to get the read name
                attributes = {item.split("=")[0]: item.split("=")[1] for item in fastq_header.split() if "=" in item}
                
                # Store barcodes and UMI if they exist
                barcodes = attributes.get("barcodes", None)
                umi = attributes.get("UMI", None)
                read_tags[read_name] = {"CB": barcodes, "ZU": umi}

    return read_tags

@log_errors
def add_tags_to_sam(input_sam, output_sam, read_tags):
    """Add tags to SAM file based on the preprocessed FASTQ headers."""
    logger = logging.getLogger('scarecrow')
    with pysam.AlignmentFile(input_sam, "r") as infile, pysam.AlignmentFile(output_sam, "w", header=infile.header) as outfile:
        for read in infile:
            if read.query_name in read_tags:
                #logger.info(f"{read.query_name} in read_tags")
                tags = read_tags[read.query_name]
                if tags["CB"]:
                    read.set_tag("CB", tags["CB"], value_type="Z")
                if tags["UR"]:
                    read.set_tag("UR", tags["UR"], value_type="Z")
            #else:
                #logger.info(f"{read.query_name} not in in read_tags")
            outfile.write(read)