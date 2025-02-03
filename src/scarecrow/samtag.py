#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
Refactored for multiprocessing batch processing
"""

import os
import sys
import time
import logging
import pysam
import multiprocessing as mp
from typing import Dict, Generator, Optional, List, Tuple
from functools import partial
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
    subparser.add_argument(
        "-o", "--out",
        metavar="<file>",
        help=("Path to SAM file to output"),
        type=str,
        default=None,
    )
    subparser.add_argument(
        "-b", "--batch_size",
        metavar="<int>",
        help=("Number of read pairs per batch to process at a time [10000]"),
        type=int,
        default=10000,
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
    run_samtag(fastq_file=args.fastq, sam_file=args.sam, out_file=args.out, batch_size = args.batch_size)

@log_errors
def run_samtag(fastq_file: str = None, sam_file: str = None, out_file: str = None, batch_size: int = None) -> None:
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
        
        # Extract tags from FASTQ (in a memory-efficient manner)
        read_tags = stream_fastq_tags(fastq_file)
        
        # Process files with minimal overhead
        if threads is None:
            threads = min(mp.cpu_count() - 1, 8)
        else:
            threads = min(threads, mp.cpu_count())
        logger.info(f"Using {threads} threads")

        # Process SAM/BAM with multiprocessing
        process_sam_multiprocessing(
            sam_file, 
            out_file, 
            read_tags,
            num_processes = threads,
            batch_size = batch_size
        )
        
        logger.info("Processing completed successfully.")
    
    except Exception as e:
        logger.error(f"Error processing files: {e}")
        sys.exit(1)



def get_total_reads(input_path: str) -> int:
    """
    Count total number of reads in the SAM/BAM file.
    
    Args:
    input_path: Path to input SAM/BAM file
    
    Returns:
    Total number of reads in the file
    """
    with pysam.AlignmentFile(input_path, 'rb') as infile:
        return infile.count()

def stream_fastq_tags(fastq_path: str) -> Dict[str, Dict[str, Optional[str]]]:
    """
    Extract all tags from FASTQ file with minimal memory usage.
    
    Returns a dictionary of read tags.
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Extracting tags from FASTQ: {fastq_path}")
    
    tag_names = ["CR", "CY", "CB", "XP", "XM", "UR", "UY"]
    read_tags = {}
    
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
                
                # Store tags for the read
                read_tags[read_name] = {
                    tag: tags.get(tag) for tag in tag_names
                }
    
    logger.info(f"Extracted tags for {len(read_tags)} reads")
    return read_tags

def process_read_batch(
    batch: List[Tuple], 
    read_tags: Dict[str, Dict[str, Optional[str]]]
) -> List[Tuple]:
    """
    Process a batch of reads, adding tags from the read_tags dictionary.
    
    Args:
    batch: List of reads to process
    read_tags: Dictionary of read tags
    
    Returns:
    List of processed reads with tags added
    """
    tag_names = ["CR", "CY", "CB", "XP", "XM", "UR", "UY"]
    processed_batch = []
    
    for read in batch:
        # Check if read name exists in tags
        if read.query_name in read_tags:
            tags = read_tags[read.query_name]
            
            # Add tags if they exist
            for tag_name in tag_names:
                tag_value = tags.get(tag_name)
                if tag_value:
                    read.set_tag(tag_name, tag_value, value_type='Z')
        
        processed_batch.append(read)
    
    return processed_batch

def process_sam_multiprocessing(
    input_path: str, 
    output_path: str, 
    read_tags: Dict[str, Dict[str, Optional[str]]],
    num_processes: Optional[int] = None,
    batch_size: int = 10000
):
    """
    Process SAM/BAM file using multiprocessing with low memory overhead.
    
    Args:
    input_path: Path to input SAM/BAM file
    output_path: Path to output SAM/BAM file
    read_tags: Dictionary of read tags
    num_processes: Number of processes to use (defaults to CPU count)
    batch_size: Number of reads to process in each batch
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Processing SAM/BAM: {input_path}")
    
    # Determine number of processes
    if num_processes is None:
        num_processes = mp.cpu_count()
    logger.info(f"Using {num_processes} processes")
    
    # Get total number of reads for progress tracking
    total_reads = get_total_reads(input_path)
    logger.info(f"Total reads to process: {total_reads}")
    
    # Start time for performance tracking
    start_time = time.time()
    
    # Open input and output files
    with pysam.AlignmentFile(input_path, 'rb') as infile, \
         pysam.AlignmentFile(output_path, 'wb', header=infile.header) as outfile:
        
        # Create process pool
        with mp.Pool(processes=num_processes) as pool:
            # Partial function for batch processing
            process_batch_func = partial(process_read_batch, read_tags=read_tags)
            
            # Process reads in batches
            batch = []
            processed_reads = 0
            for read in infile:
                batch.append(read)
                
                # Process batch when it reaches the specified size
                if len(batch) >= batch_size:
                    # Process batch in parallel
                    processed_batches = pool.map(process_batch_func, [batch])
                    
                    # Write processed reads
                    for processed_batch in processed_batches:
                        for processed_read in processed_batch:
                            outfile.write(processed_read)
                    
                    # Update progress
                    processed_reads += len(batch)
                    elapsed_time = time.time() - start_time
                    reads_per_sec = processed_reads / elapsed_time if elapsed_time > 0 else 0
                    estimated_remaining = (total_reads - processed_reads) / reads_per_sec if reads_per_sec > 0 else 0
                    
                    logger.info(f"Processed {processed_reads}/{total_reads} reads "
                                f"({processed_reads/total_reads*100:.2f}%) "
                                f"| Speed: {reads_per_sec:.2f} reads/sec "
                                f"| Est. remaining: {estimated_remaining:.2f} sec")
                    
                    # Reset batch
                    batch = []
            
            # Process any remaining reads
            if batch:
                processed_batches = pool.map(process_batch_func, [batch])
                for processed_batch in processed_batches:
                    for processed_read in processed_batch:
                        outfile.write(processed_read)
                
                # Update final progress
                processed_reads += len(batch)
                elapsed_time = time.time() - start_time
                logger.info(f"Processed {processed_reads}/{total_reads} reads "
                            f"({processed_reads/total_reads*100:.2f}%) "
                            f"| Total time: {elapsed_time:.2f} sec")
    
    logger.info("Multiprocessing complete")
