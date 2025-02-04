#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
Refactored for on-demand FASTQ tag extraction and multiprocessing
"""

import os
import sys
import time
import logging
import pysam
import mmap
import multiprocessing as mp
from typing import Dict, Optional, List, Tuple, Union
from functools import partial
from pyfastx import Fastq
from argparse import RawTextHelpFormatter
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string

def parser_samtag(parser):
    """
    Add samtag subparser with command-line arguments
    """
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
        required=True,
    )
    subparser.add_argument(
        "-s", "--sam",
        metavar="<file>",
        help=("Path to SAM file to update tags"),
        type=str,
        required=True,
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
        help=("Number of read pairs per batch to process at a time [20000]"),
        type=int,
        default=20000,
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
    """
    Validate and run samtag processing
    """
    run_samtag(
        fastq_file=args.fastq, 
        sam_file=args.sam, 
        out_file=args.out, 
        batch_size=args.batch_size,
        threads=args.threads
    )

@log_errors
def run_samtag(
    fastq_file: str = None, 
    sam_file: str = None, 
    out_file: str = None, 
    batch_size: int = 200000,
    threads: Optional[int] = None
) -> None:
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
        
        # Determine number of threads
        if threads is None:
            threads = min(mp.cpu_count() - 1, 8)
        else:
            threads = min(threads, mp.cpu_count())
        logger.info(f"Using {threads} threads")
        logger.info(f"Processing reads in batches of {batch_size}")

        # Create temporary directory for intermediate files
        temp_dir = f"samtag_temp_{generate_random_string()}"
        os.makedirs(temp_dir, exist_ok=True)
        logger.info(f"Using temporary directory: {temp_dir}")

        # Process SAM/BAM with multiprocessing
        process_sam_multiprocessing(
            sam_file, 
            out_file, 
            fastq_file,
            temp_dir,
            num_processes = threads,
            chunk_size = batch_size
        )
        
        # Clean up temporary directory
        import shutil
        shutil.rmtree(temp_dir)
        
        logger.info("Processing completed successfully.")
    
    except Exception as e:
        logger.error(f"Error processing files: {e}")
        sys.exit(1)

def extract_tags_from_fastq(fastq_path: str, read_names: List[str]) -> Dict[str, Dict[str, str]]:
    """
    Extract tags for a specific set of read names from the FASTQ file.
    """
    tag_names = ["CR", "CY", "CB", "XP", "XM", "UR", "UY"]
    read_tags = {}
    read_names_set = set(read_names)

    with pysam.FastxFile(fastq_path) as fastq:
        for entry in fastq:
            read_name = entry.name.split()[0]  # Extract read name from header
            if read_name in read_names_set:
                tags = {}
                for item in entry.comment.split():  # Extract tags from the comment field
                    if '=' in item:
                        key, value = item.split('=')
                        tags[key] = value
                read_tags[read_name] = {tag: tags.get(tag) for tag in tag_names}
                read_names_set.remove(read_name)
                if not read_names_set:
                    break  # Stop once all required reads are found

    return

@log_errors
def process_chunk(args: Tuple) -> Tuple[bool, str]:
    """
    Process a chunk of the BAM file and write tagged reads to a temporary file.
    """
    #logger = logging.getLogger('scarecrow')
    logger = setup_worker_logger()
    input_path, fastq_path, header_dict, chunk_reads, temp_file = args

    try:
        # Extract tags for the reads in this chunk
        read_names = [read.query_name for read in chunk_reads]
        read_tags = extract_tags_from_fastq(fastq_path, read_names)

        # Write tagged reads to the temporary file
        with pysam.AlignmentFile(temp_file, 'wb', header=header_dict) as outfile:
            for read in chunk_reads:
                read_name = read.query_name
                if read_name in read_tags:
                    tags = read_tags[read_name]
                    for tag_name, tag_value in tags.items():
                        if tag_value:
                            read.set_tag(tag_name, tag_value, value_type='Z')
                outfile.write(read)

        return True, temp_file

    except Exception as e:
        logger.warning(f"Error processing chunk {temp_file}: {e}")
        return False, str(e)

@log_errors
def process_sam_multiprocessing(
    input_path: str,
    output_path: str,
    fastq_path: str,
    temp_dir: str,
    num_processes: int = 4,
    chunk_size: int = 100000
):
    """
    Process the BAM file in parallel using multiprocessing.
    """
    logger = setup_worker_logger()
    logger.info("Starting sequence extraction")

    # Get BAM file header
    with pysam.AlignmentFile(input_path, 'rb') as infile:
        header_dict = infile.header.to_dict()

    # Create temporary directory
    os.makedirs(temp_dir, exist_ok=True)

    # Read the entire BAM file and split into chunks
    logger.info("Reading BAM file and splitting into chunks...")
    chunks = []
    current_chunk = []
    with pysam.AlignmentFile(input_path, 'rb') as infile:
        for read in infile:
            current_chunk.append(read)
            if len(current_chunk) >= chunk_size:
                chunks.append(current_chunk)
                current_chunk = []
        if current_chunk:  # Add the last chunk
            chunks.append(current_chunk)

    logger.info(f"Split into {len(chunks)} chunks")

    # Prepare arguments for each chunk
    process_args = [
        (input_path, fastq_path, header_dict, chunk, os.path.join(temp_dir, f"temp_{i}.bam"))
        for i, chunk in enumerate(chunks)
    ]

    # Process chunks in parallel
    logger.info(f"Processing {len(chunks)} chunks...")
    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(process_chunk, process_args)

    # Merge temporary files
    logger.info("Merging temporary files...")
    with pysam.AlignmentFile(output_path, 'wb', header=header_dict) as outfile:
        for success, temp_file in results:
            if success:
                with pysam.AlignmentFile(temp_file, 'rb') as infile:
                    for read in infile:
                        outfile.write(read)
                os.remove(temp_file)  # Clean up temporary file
            else:
                logger.error(f"Failed to process temporary file: {temp_file}")

    logger.info(f"Processing complete. Output written to {output_path}")

def setup_worker_logger(log_file: str = None):
    """Configure logger for worker processes with file output"""
    logger = logging.getLogger('scarecrow')
    if not logger.handlers:  # Only add handlers if none exist
        # Create formatters
        formatter = logging.Formatter('%(asctime)s - %(processName)s - %(name)s - %(levelname)s - %(message)s')
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
        
        # File handler if log_file provided
        if log_file:
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
        
        logger.setLevel(logging.INFO)
    return logger