# -*- coding: utf-8 -*-
"""
#!/usr/bin/env python3
@author: David Wragg
"""

import logging
import pysam
import os
from argparse import RawTextHelpFormatter
from typing import Dict, Generator
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
    run_samtag(fastq_file=args.fastq, sam_file=args.sam)
    
@log_errors
def run_samtag(fastq_file: str = None, sam_file: str = None) -> None:
    """
    Function to process barcodes with minimal memory usage
    """
    # Setup logging
    logfile = f'./scarecrow_samtag_{generate_random_string()}.log'
    logger = setup_logger(logfile)
    logger.info(f"logfile: '{logfile}'")
    
    # Preprocess fastq headers
    if fastq_file and sam_file:
        logger.info(f"Processing {fastq_file} and {sam_file}")
        name, _ = os.path.splitext(sam_file)
        output = name + ".tagged.sam"
        
        # Use a generator to iterate through read tags
        read_tag_generator = fastq_header_tags(fastq_file)
        
        # Process SAM file with minimal memory
        process_sam_with_tags(sam_file, output, read_tag_generator)

def fastq_header_tags(fastq_file: str) -> Generator[Dict[str, str], None, None]:
    """
    Generate read tags from FASTQ headers with minimal memory usage.
    Yields a dictionary of tags for each read.
    """
    with open(fastq_file, "r") as fq:
        while True:
            header = fq.readline().strip()
            if not header:
                break
            
            # Skip the next 3 lines (sequence, '+', quality)
            fq.readline()
            fq.readline()
            fq.readline()
            
            if header.startswith("@"):
                # Extract the read name and tags
                read_name = header.split(" ")[0][1:]  # Remove "@" and split to get the read name
                attributes = {item.split("=")[0]: item.split("=")[1] 
                              for item in header.split() if "=" in item}
                
                # Prepare the tags dictionary
                tags = {
                    "read_name": read_name,
                    "CR": attributes.get("CR"),
                    "CY": attributes.get("CY"),
                    "CB": attributes.get("CB"),
                    "XP": attributes.get("XP"),
                    "XM": attributes.get("XM"),
                    "UR": attributes.get("UR"),
                    "UY": attributes.get("UY")
                }
                
                yield tags

def process_sam_with_tags(input_sam: str, output_sam: str, read_tag_generator: Generator):
    """
    Process SAM file and add tags with minimal memory usage.
    """
    logger = logging.getLogger('scarecrow')
    logger.info(f"Processing SAM file {input_sam}")
    
    # Create a dictionary to store read tags
    read_tags = {}
    
    # Pre-load read tags into memory
    for tag_dict in read_tag_generator:
        read_tags[tag_dict['read_name']] = tag_dict
    
    # Process SAM file
    with pysam.AlignmentFile(input_sam, "r") as infile, \
         pysam.AlignmentFile(output_sam, "w", header=infile.header) as outfile:
        
        for read in infile:
            # Check if read name exists in tags
            if read.query_name in read_tags:
                tags = read_tags[read.query_name]
                
                # Add tags if they exist
                tag_mapping = [
                    ("CR", "Z"), ("CY", "Z"), ("CB", "Z"),
                    ("XP", "Z"), ("XM", "Z"), 
                    ("UR", "Z"), ("UY", "Z")
                ]
                
                for tag_name, val_type in tag_mapping:
                    tag_value = tags.get(tag_name)
                    if tag_value:
                        read.set_tag(tag_name, tag_value, value_type=val_type)
            
            # Write the read to output
            outfile.write(read)
    
    logger.info(f"Completed processing. Output written to {output_sam}")