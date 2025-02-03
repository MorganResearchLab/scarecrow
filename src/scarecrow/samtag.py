#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
Refactored for multiprocessing batch processing
"""

import logging
import pysam
import os
import multiprocessing as mp
from functools import partial
from argparse import RawTextHelpFormatter
from typing import List, Dict, Tuple
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
    subparser.add_argument(
        "-o", "--out",
        metavar="<file>",
        help=("Path to SAM file to output"),
        type=str,
        default=None,
    )
    subparser.add_argument(
        "-m", "--max-read-buffer", 
        metavar="<int>",       
        type=int, 
        default=10000, 
        help="Maximum number of reads to buffer (default: 10000)"
    )
    
    return subparser

def validate_samtag_args(parser, args):
    run_samtag(fastq_file=args.fastq, sam_file=args.sam, out_file=args.out, max_read_buffer = args.max_read_buffer)

@log_errors
def run_samtag(fastq_file: str = None, sam_file: str = None, out_file: str = None, max_read_buffer: int = None) -> None:
    """
    Multiprocessing function to process SAM tags efficiently
    """
    # Setup logging
    logfile = f'./scarecrow_samtag_{generate_random_string()}.log'
    logger = setup_logger(logfile)
    logger.info(f"logfile: '{logfile}'")

    try:
        # Determine output path
        if not out_file:
            base, ext = os.path.splitext(sam_file)
            out_file = f"{base}.tagged{ext}"
        
        # Stream tags from FASTQ
        tag_generator = stream_fastq_tags(fastq_file)
        
        # Process SAM/BAM with tags
        process_sam_with_tags(
            sam_file, 
            out_file, 
            tag_generator,
            max_read_buffer = max_read_buffer
        )
        
        logger.info("Processing completed successfully.")
    
    except Exception as e:
        logger.error(f"Error processing files: {e}")
        sys.exit(1)



def stream_fastq_tags(fastq_path: str) -> Generator[Dict[str, Optional[str]], None, None]:
    """
    Stream tags from FASTQ file with minimal memory usage.
    
    Yields dictionaries of tags for each read, preserving memory efficiency.
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Extracting tags from FASTQ: {fastq_path}")
    
    tag_names = ["CR", "CY", "CB", "XP", "XM", "UR", "UY"]
    
    with open(fastq_path, 'r') as fq:
        while True:
            # Read header line
            header = fq.readline().strip()
            if not header:
                break
            
            # Skip sequence, '+', and quality lines
            fq.readline()
            fq.readline()
            fq.readline()
            
            # Process only FASTQ headers starting with '@'
            if header.startswith('@'):
                # Extract read name (first part of header)
                read_name = header.split()[0][1:]
                
                # Extract tags
                tags = {}
                for item in header.split():
                    if '=' in item:
                        key, value = item.split('=')
                        tags[key] = value
                
                # Prepare tag dictionary with specified tag names
                yield {
                    'read_name': read_name,
                    **{tag: tags.get(tag) for tag in tag_names}
                }

def process_sam_with_tags(
    input_path: str, 
    output_path: str, 
    tag_generator: Generator[Dict[str, Optional[str]], None, None],
    max_read_buffer: int = 10000
):
    """
    Process SAM/BAM file with extremely low memory usage.
    
    Streams reads, looks up tags, and writes to output with minimal memory overhead.
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Processing SAM/BAM: {input_path}")
    
    # Prepare tag lookup
    tag_cache: Dict[str, Dict[str, Optional[str]]] = {}
    tag_names = ["CR", "CY", "CB", "XP", "XM", "UR", "UY"]
    
    # Open input and output files
    with pysam.AlignmentFile(input_path, 'rb') as infile, \
         pysam.AlignmentFile(output_path, 'wb', header=infile.header) as outfile:
        
        # Stream tags and build a small, rotating cache
        for tag_entry in tag_generator:
            read_name = tag_entry.pop('read_name')
            tag_cache[read_name] = tag_entry
            
            # Limit cache size
            if len(tag_cache) > max_read_buffer:
                oldest_read = min(tag_cache.keys())
                del tag_cache[oldest_read]
        
        # Reset file pointer
        infile.reset()
        
        # Process reads
        for read in infile:
            # Check if read name exists in current tag cache
            if read.query_name in tag_cache:
                tags = tag_cache[read.query_name]
                
                # Add tags if they exist
                for tag_name in tag_names:
                    tag_value = tags.get(tag_name)
                    if tag_value:
                        read.set_tag(tag_name, tag_value, value_type='Z')
            
            # Write read to output
            outfile.write(read)
        
        logger.info(f"Completed processing: {output_path}")

