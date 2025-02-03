#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
"""

import logging
import pysam
import os
import tempfile
from argparse import RawTextHelpFormatter
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
    
    # Validate inputs
    if not fastq_file or not sam_file:
        logger.error("Both FASTQ and SAM files must be provided")
        return
    
    # Prepare output filename
    name, _ = os.path.splitext(sam_file)
    output = name + ".tagged.sam"
    
    # Create temporary file for read tags
    with tempfile.NamedTemporaryFile(mode='w+', delete=True) as tag_temp_file:
        # Step 1: Extract tags from FASTQ file
        extract_fastq_tags(fastq_file, tag_temp_file)
        tag_temp_file.flush()
        tag_temp_file.seek(0)
        
        # Step 2: Update SAM file with tags
        process_sam_with_tags(sam_file, output, tag_temp_file)
    
    logger.info(f"Processing complete. Output written to {output}")

def extract_fastq_tags(fastq_file: str, tag_temp_file):
    """
    Extract read tags from FASTQ file and write to a temporary file.
    
    Temporary file format: 
    read_name\tCR\tCY\tCB\tXP\tXM\tUR\tUY
    """
    logger = logging.getLogger('scarecrow')
    logger.info(f"Extracting tags from {fastq_file}")
    
    tag_names = ["CR", "CY", "CB", "XP", "XM", "UR", "UY"]
    
    with open(fastq_file, "r") as fq:
        while True:
            # Read header line
            header = fq.readline().strip()
            if not header:
                break
            
            # Skip the next 3 lines (sequence, '+', quality)
            fq.readline()
            fq.readline()
            fq.readline()
            
            if header.startswith("@"):
                # Extract read name
                read_name = header.split(" ")[0][1:]  # Remove "@"
                
                # Extract tag values
                tag_dict = {}
                for item in header.split():
                    if "=" in item:
                        key, value = item.split("=")
                        tag_dict[key] = value
                
                # Prepare tag line
                tag_values = [tag_dict.get(tag, "") for tag in tag_names]
                tag_line = f"{read_name}\t{'\t'.join(tag_values)}\n"
                tag_temp_file.write(tag_line)

def process_sam_with_tags(input_sam: str, output_sam: str, tag_temp_file):
    """
    Process SAM file and add tags with absolute minimal memory usage.
    """
    logger = logging.getLogger('scarecrow')
    logger.info(f"Processing SAM file {input_sam}")
    
    # Create tag lookup function
    def get_read_tags(read_name):
        tag_temp_file.seek(0)
        for line in tag_temp_file:
            parts = line.strip().split('\t')
            if parts[0] == read_name:
                return parts[1:]
        return None
    
    # Tag names in order
    tag_names = ["CR", "CY", "CB", "XP", "XM", "UR", "UY"]
    
    # Process SAM file
    with pysam.AlignmentFile(input_sam, "r") as infile, \
         pysam.AlignmentFile(output_sam, "w", header=infile.header) as outfile:
        
        for read in infile:
            # Get tags for this read
            tags = get_read_tags(read.query_name)
            
            if tags:
                # Add tags if they exist
                for tag_name, tag_value in zip(tag_names, tags):
                    if tag_value:
                        read.set_tag(tag_name, tag_value, value_type="Z")
            
            # Write the read to output
            outfile.write(read)
    
    logger.info(f"Completed processing. Output written to {output_sam}")