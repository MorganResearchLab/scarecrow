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
        "-@", "--threads",
        metavar="<int>",
        help=("Number of processing threads [1]"),
        type=int,
        default=1,
    )
    return subparser

def validate_samtag_args(parser, args):
    run_samtag(fastq_file=args.fastq, sam_file=args.sam, threads=args.threads)

@log_errors
def run_samtag(fastq_file: str = None, sam_file: str = None, threads: int = None) -> None:
    """
    Multiprocessing function to process SAM tags efficiently
    """
    # Setup logging
    logfile = f'./scarecrow_samtag_{generate_random_string()}.log'
    logger = setup_logger(logfile)
    logger.info(f"logfile: '{logfile}'")
    
    # Validate inputs
    if not fastq_file or not sam_file:
        logger.error("Both FASTQ and SAM files must be provided")
        return
    
    # Process files with minimal overhead
    if threads is None:
        threads = min(mp.cpu_count() - 1, 8)
    else:
        threads = min(threads, mp.cpu_count())
    logger.info(f"Using {threads} threads")

    # Prepare output filename
    name, _ = os.path.splitext(sam_file)
    output = name + ".tagged.sam"
    
    # Extract read tags from FASTQ
    read_tags = extract_fastq_tags(fastq_file)
    logger.info(f"Extracted {len(read_tags)} read tags from FASTQ")
    
    # Process SAM file in batches
    process_sam_with_tags_parallel(sam_file, output, read_tags, threads)
    
    logger.info(f"Processing complete. Output written to {output}")

def extract_fastq_tags(fastq_file: str) -> Dict[str, Dict[str, str]]:
    """
    Extract all read tags from FASTQ file into a dictionary
    
    Returns a dictionary with read names as keys and tag dictionaries as values
    """
    read_tags = {}
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
                read_name = header.split(" ")[0][1:]  # Remove "@"
                
                # Parse tags from header
                tags = {}
                for item in header.split():
                    if "=" in item:
                        key, value = item.split("=")
                        tags[key] = value
                
                read_tags[read_name] = {
                    "CR": tags.get("CR"),
                    "CY": tags.get("CY"),
                    "CB": tags.get("CB"),
                    "XP": tags.get("XP"),
                    "XM": tags.get("XM"),
                    "UR": tags.get("UR"),
                    "UY": tags.get("UY")
                }
    
    return read_tags

def process_batch(batch: List[Tuple], read_tags: Dict[str, Dict[str, str]]) -> List[Tuple]:
    """
    Process a batch of reads and add tags
    
    Args:
    batch: List of reads (as pysam AlignedSegment objects or their equivalent)
    read_tags: Dictionary of read tags
    
    Returns:
    List of processed reads with tags added
    """
    processed_batch = []
    for read in batch:
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
        
        processed_batch.append(read)
    
    return processed_batch

def process_sam_with_tags_parallel(input_sam: str, output_sam: str, 
                                    read_tags: Dict[str, Dict[str, str]], 
                                    threads: int):
    """
    Process SAM file in parallel using multiprocessing
    """
    logger = logging.getLogger('scarecrow')
    logger.info(f"Processing SAM file {input_sam} in parallel")
    
    # Determine batch size (adjust based on your system's capabilities)
    batch_size = 10000
    
    # Open input and output SAM files
    with pysam.AlignmentFile(input_sam, "r") as infile, \
         pysam.AlignmentFile(output_sam, "w", header=infile.header) as outfile:
        
        # Create a pool of worker processes
        with mp.Pool(processes=threads) as pool:
            # Partial function with read tags
            process_batch_func = partial(process_batch, read_tags=read_tags)
            
            # Batch reads and process
            batch = []
            for read in infile:
                batch.append(read)
                
                # Process batch when it reaches the specified size
                if len(batch) >= batch_size:
                    # Use map to process batches in parallel
                    processed_batches = pool.map(process_batch_func, [batch])
                    
                    # Write processed reads
                    for processed_batch in processed_batches:
                        for processed_read in processed_batch:
                            outfile.write(processed_read)
                    
                    # Reset batch
                    batch = []
            
            # Process any remaining reads
            if batch:
                processed_batches = pool.map(process_batch_func, [batch])
                for processed_batch in processed_batches:
                    for processed_read in processed_batch:
                        outfile.write(processed_read)
    
    logger.info(f"Parallel processing complete")