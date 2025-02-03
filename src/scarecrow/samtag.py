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
import multiprocessing as mp
from typing import Dict, Optional, List, Tuple, Union
from functools import partial
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

def get_total_reads(input_path: str) -> int:
    """
    Count total number of reads in the SAM/BAM file.
    
    Handles both indexed and non-indexed files.
    """
    try:
        # Try to use indexed counting first
        with pysam.AlignmentFile(input_path, 'rb') as infile:
            return infile.count()
    except ValueError:
        # If indexing fails, count manually
        total_reads = 0
        with pysam.AlignmentFile(input_path, 'rb') as infile:
            for _ in infile:
                total_reads += 1
        return total_reads

@log_errors
def run_samtag(
    fastq_file: str = None, 
    sam_file: str = None, 
    out_file: str = None, 
    batch_size: int = 10000,
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
            batch_size = batch_size
        )
        
        # Clean up temporary directory
        import shutil
        shutil.rmtree(temp_dir)
        
        logger.info("Processing completed successfully.")
    
    except Exception as e:
        logger.error(f"Error processing files: {e}")
        sys.exit(1)

def extract_batch_tags(
    batch_read_names: List[str], 
    fastq_path: str
) -> Dict[str, Dict[str, Optional[str]]]:
    """
    Extract tags for a specific batch of reads from FASTQ file.
    
    Args:
    batch_read_names: List of read names to extract tags for
    fastq_path: Path to FASTQ file
    
    Returns:
    Dictionary of read tags for the requested reads
    """
    tag_names = ["CR", "CY", "CB", "XP", "XM", "UR", "UY"]
    read_tags = {}
    read_names_set = set(batch_read_names)
    
    with open(fastq_path, 'r') as fq:
        while True:
            # Read header line
            header = fq.readline().strip()
            if not header:
                break
            
            # Skip sequence, '+', and quality lines
            seq_line = fq.readline()
            plus_line = fq.readline()
            qual_line = fq.readline()
            
            # Process only FASTQ headers starting with '@'
            if header.startswith('@'):
                # Extract read name (first part of header)
                read_name = header.split()[0][1:]
                
                # Check if this read is in the current batch
                if read_name in read_names_set:
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
                    
                    # Remove from set to reduce iteration time
                    read_names_set.remove(read_name)
                    
                    # Exit early if all batch reads processed
                    if not read_names_set:
                        break
    
    return read_tags

def serialize_read(read: Union[pysam.AlignedSegment, None]) -> Tuple[Optional[str], Optional[List]]:
    """
    Convert a pysam AlignedSegment to a picklable representation.
    
    Args:
    read: A pysam AlignedSegment or None
    
    Returns:
    A tuple of (read name, [serialized read data])
    """
    if read is None:
        return None, None
    
    # Extract core read information
    return read.query_name, [
        read.query_sequence,
        read.qual,
        read.flag,
        read.reference_name,
        read.reference_start,
        read.mapping_quality,
        list(read.get_tags())
    ]

def deserialize_read(
    header: pysam.AlignmentHeader, 
    serialized_read: Tuple[Optional[str], Optional[List]]
) -> Optional[pysam.AlignedSegment]:
    """
    Reconstruct a pysam AlignedSegment from its serialized representation.
    
    Args:
    header: The AlignmentHeader from the original file
    serialized_read: Tuple of (read name, [serialized read data])
    
    Returns:
    A reconstructed pysam AlignedSegment
    """
    if serialized_read is None or serialized_read[0] is None:
        return None
    
    read_name, read_data = serialized_read
    
    # Recreate the AlignedSegment
    read = pysam.AlignedSegment(header)
    read.query_name = read_name
    read.query_sequence = read_data[0]
    read.qual = read_data[1]
    read.flag = read_data[2]
    read.reference_name = read_data[3]
    read.reference_start = read_data[4]
    read.mapping_quality = read_data[5]
    
    # Restore tags
    for tag in read_data[6]:
        read.set_tag(tag[0], tag[1], tag[2])
    
    return read

def process_read_batch(
    serialized_batch: List[Tuple[Optional[str], Optional[List]]], 
    fastq_path: str,
    header_dict: Dict[str, Union[str, int]]
) -> List[Tuple[Optional[str], Optional[List]]]:
    """
    Process a batch of reads by extracting and adding tags from FASTQ.
    
    Args:
    serialized_batch: List of serialized reads
    fastq_path: Path to FASTQ file
    header_dict: Dictionary representation of the BAM header
    
    Returns:
    List of processed serialized reads with tags added
    """
    # Reconstruct header
    header = pysam.AlignmentHeader.from_dict(header_dict)
    print(f"\nheader:\n{header}")

    # Deserialize batch
    batch = [
        deserialize_read(header, read_data) 
        for read_data in serialized_batch
    ]
    
    # Extract read names
    batch_read_names = [read.query_name for read in batch if read is not None]
    print(f"\nbatch_read_names:\n{batch_read_names}")
    
    # Get tags for this batch
    print(f"\nbatch_read_names:\n{batch_read_names}")
    read_tags = extract_batch_tags(batch_read_names, fastq_path)
    print(f"\nread_tags:\n{read_tags}")

    tag_names = ["CR", "CY", "CB", "XP", "XM", "UR", "UY"]
    processed_batch = []
    
    for read in batch:
        if read is not None:
            # Check if read name exists in tags
            if read.query_name in read_tags:
                tags = read_tags[read.query_name]
                
                # Add tags if they exist
                for tag_name in tag_names:
                    tag_value = tags.get(tag_name)
                    if tag_value:
                        read.set_tag(tag_name, tag_value, value_type='Z')
            
            # Serialize processed read
            processed_batch.append(serialize_read(read))
    
    return processed_batch

def process_sam_multiprocessing(
    input_path: str, 
    output_path: str, 
    fastq_path: str,
    temp_dir: str,
    num_processes: Optional[int] = None,
    batch_size: int = 10000
):
    """
    Process SAM/BAM file using multiprocessing with on-demand FASTQ tag extraction.
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
    
    # Prepare multiprocessing pool and intermediate files
    temp_files = [os.path.join(temp_dir, f"batch_{i}.sam") for i in range(num_processes)]
    
    # Open input file and get header
    with pysam.AlignmentFile(input_path, 'rb') as infile:
        # Convert header to dictionary for serialization
        header_dict = infile.header.to_dict()
        
        # Create process pool
        with mp.Pool(processes=num_processes) as pool:
            # Create temporary output files for each process
            temp_outfiles = [
                pysam.AlignmentFile(temp_file, 'w', header=infile.header) 
                for temp_file in temp_files
            ]
            
            try:
                # Process reads in batches
                batch = []
                processed_reads = 0
                
                for i, read in enumerate(infile):
                    # Serialize read to make it picklable
                    serialized_read = serialize_read(read)
                    batch.append(serialized_read)
                    
                    # Process batch when it reaches the specified size
                    if len(batch) >= batch_size:

                        # Distribute batches across processes
                        process_batch_func = partial(
                            process_read_batch, 
                            fastq_path = fastq_path,
                            header_dict = header_dict
                        )
                        
                        # Distribute batch processing
                        processed_batches = pool.map(process_batch_func, [batch])
                        
                        # Write processed reads to thread-specific files
                        for process_idx, processed_batch in enumerate(processed_batches):
                            for processed_read_data in processed_batch:
                                print(f"{process_idx} : {processed_batch}")
                                processed_read = deserialize_read(infile.header, processed_read_data)
                                if processed_read:
                                    temp_outfiles[process_idx].write(processed_read)
                        
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
                    process_batch_func = partial(
                        process_read_batch, 
                        fastq_path=fastq_path,
                        header_dict=header_dict
                    )
                    processed_batches = pool.map(process_batch_func, [batch])
                    
                    # Write processed reads to thread-specific files
                    for process_idx, processed_batch in enumerate(processed_batches):
                        for processed_read_data in processed_batch:
                            processed_read = deserialize_read(infile.header, processed_read_data)
                            if processed_read:
                                temp_outfiles[process_idx].write(processed_read)
                    
                    processed_reads += len(batch)
                    elapsed_time = time.time() - start_time
                    reads_per_sec = processed_reads / elapsed_time if elapsed_time > 0 else 0

                # Close temporary output files
                for outfile in temp_outfiles:
                    outfile.close()
                    
                # Concatenate temporary files
                with pysam.AlignmentFile(output_path, 'w', header=infile.header) as final_outfile:
                    for temp_file in temp_files:
                        with pysam.AlignmentFile(temp_file, 'r') as temp_infile:
                            for temp_read in temp_infile:
                                final_outfile.write(temp_read)
                    
                # Log final progress
                elapsed_time = time.time() - start_time
                logger.info(f"Processed {processed_reads}/{total_reads} reads "
                            f"({processed_reads/total_reads*100:.2f}%) "
                            f"| Total time: {elapsed_time:.2f} sec")
            
            finally:
                # Ensure files are closed
                for outfile in temp_outfiles:
                    if not outfile.is_closed:
                        outfile.close()
    
    logger.info("Multiprocessing complete")