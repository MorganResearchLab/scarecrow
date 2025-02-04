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
        help=("Number of read pairs per batch to process at a time [200000]"),
        type=int,
        default=200000,
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
            batch_size = batch_size
        )
        
        # Clean up temporary directory
        import shutil
        shutil.rmtree(temp_dir)
        
        logger.info("Processing completed successfully.")
    
    except Exception as e:
        logger.error(f"Error processing files: {e}")
        sys.exit(1)

def extract_batch_tags(batch_read_names: List[str], fastq_path: str) -> Dict[str, Dict[str, str]]:
    """
    Optimized tag extraction using memory mapping and set lookups.
    """
    tag_names = ["CR", "CY", "CB", "XP", "XM", "UR", "UY"]
    read_tags = {}
    read_names_set = set(batch_read_names)
    
    # Memory map the FASTQ file
    with open(fastq_path, 'rb') as f:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        
        current_pos = 0
        while current_pos < mm.size():
            # Read header line
            header_end = mm.find(b'\n', current_pos)
            if header_end == -1:
                break
                
            header = mm[current_pos:header_end].decode('utf-8')
            
            # Skip to next record
            for _ in range(3):  # Skip sequence, '+', and quality lines
                current_pos = mm.find(b'\n', current_pos) + 1
                if current_pos == 0:  # Not found
                    break
            
            if header.startswith('@'):
                read_name = header.split()[0][1:]
                if read_name in read_names_set:
                    # Extract tags
                    tags = {}
                    for item in header.split():
                        if '=' in item:
                            key, value = item.split('=')
                            tags[key] = value
                    
                    read_tags[read_name] = {
                        tag: tags.get(tag) for tag in tag_names
                    }
                    
                    read_names_set.remove(read_name)
                    if not read_names_set:
                        break
            
            current_pos = header_end + 1
        
        mm.close()
    
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
        list(read.get_tags(with_value_type=True))
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
    try:
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
    
    except IndexError:
        print(f"Index error for {read_data}")
        print("Returning read as None")
        read = None
    
    return read

def chunk_iterator(file_path: str, chunk_size: int, num_chunks: int):
    """
    Memory-efficient iterator that yields chunks of file positions with improved boundary handling.
    """
    try:
        file_size = os.path.getsize(file_path)
        chunk_positions = []
        current_pos = 0
        last_query_name = None
        reads_in_chunk = 0
        
        with pysam.AlignmentFile(file_path, 'rb', check_sq=False) as infile:
            try:
                for read in infile:
                    reads_in_chunk += 1
                    
                    # Start new chunk if:
                    # 1. We've reached chunk size AND
                    # 2. Current read has different name from last (don't split paired reads)
                    if reads_in_chunk >= chunk_size and (
                        not last_query_name or read.query_name != last_query_name
                    ):
                        next_pos = infile.tell()
                        if next_pos > current_pos:
                            chunk_positions.append((current_pos, next_pos))
                            current_pos = next_pos
                            reads_in_chunk = 0
                            
                            if len(chunk_positions) >= num_chunks:
                                yield chunk_positions
                                chunk_positions = []
                    
                    last_query_name = read.query_name
                    
            except OSError as e:
                logging.error(f"Error reading file during chunking: {e}")
                raise
            
            # Handle remaining portion
            if current_pos < file_size:
                chunk_positions.append((current_pos, file_size))
                if chunk_positions:
                    yield chunk_positions
                    
    except OSError as e:
        logging.error(f"Error accessing file {file_path}: {e}")
        raise

def process_chunk(args: tuple) -> tuple:
    """
    Process a chunk of SAM/BAM file with improved chunk boundary handling and error reporting.
    Returns tuple of (success, result/error_message)
    """
    input_path, fastq_path, header_dict, chunk_start, chunk_end, temp_file, chunk_size = args
    
    try:
        if not os.path.exists(input_path) or not os.path.exists(fastq_path):
            return (False, f"Input files not found: {input_path} or {fastq_path}")
            
        print(f"Processing chunk: start={chunk_start}, end={chunk_end}", file=sys.stderr)
        
        with pysam.AlignmentFile(input_path, 'rb', check_sq=False, require_index=False) as infile:
            header = pysam.AlignmentHeader.from_dict(header_dict)
            with pysam.AlignmentFile(temp_file, 'wb', header=header) as outfile:
                if chunk_start >= chunk_end:
                    return (False, f"Invalid chunk positions: start={chunk_start}, end={chunk_end}")
                
                try:
                    infile.seek(chunk_start)
                    batch_reads = []
                    batch_positions = []
                    current_pos = chunk_start
                    reads_processed = 0
                    last_read = None
                    
                    while True:
                        try:
                            read = next(infile)
                            current_pos = infile.tell()
                            
                            # Break if we've exceeded chunk boundary and read is different from last
                            if current_pos > chunk_end and (not last_read or read.query_name != last_read.query_name):
                                break
                                
                            reads_processed += 1
                            last_read = read
                            
                            if reads_processed % 1000 == 0:
                                print(f"Chunk progress: {reads_processed} reads, pos={current_pos}/{chunk_end}", 
                                      file=sys.stderr)
                            
                            batch_reads.append(read.query_name)
                            batch_positions.append((read, current_pos))
                            
                            # Process batch when full or at chunk boundary
                            if len(batch_reads) >= chunk_size or current_pos >= chunk_end:
                                if batch_reads:
                                    read_tags = extract_batch_tags(batch_reads, fastq_path)
                                    for read, _ in batch_positions:
                                        if read.query_name in read_tags:
                                            tags = read_tags[read.query_name]
                                            for tag_name, tag_value in tags.items():
                                                if tag_value:
                                                    read.set_tag(tag_name, tag_value, value_type='Z')
                                        outfile.write(read)
                                
                                batch_reads = []
                                batch_positions = []
                                
                                # Break if we're past chunk boundary
                                if current_pos >= chunk_end:
                                    break
                        
                        except StopIteration:
                            # Process any remaining reads in the last batch
                            if batch_reads:
                                read_tags = extract_batch_tags(batch_reads, fastq_path)
                                for read, _ in batch_positions:
                                    if read.query_name in read_tags:
                                        tags = read_tags[read.query_name]
                                        for tag_name, tag_value in tags.items():
                                            if tag_value:
                                                read.set_tag(tag_name, tag_value, value_type='Z')
                                    outfile.write(read)
                            break
                
                except Exception as e:
                    return (False, f"Error in processing loop: {str(e)}")
        
        if os.path.getsize(temp_file) == 0:
            return (False, f"Generated empty output file: {temp_file}")
            
        return (True, temp_file)
    
    except Exception as e:
        error_msg = f"Fatal error in process_chunk: {str(e)}"
        print(error_msg, file=sys.stderr)
        if os.path.exists(temp_file):
            try:
                os.remove(temp_file)
            except:
                pass
        return (False, error_msg)


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
    #print(f"\nheader:\n{header}")

    # Deserialize batch
    batch = [
        deserialize_read(header, read_data) 
        for read_data in serialized_batch
    ]
    
    # Extract read names
    batch_read_names = [read.query_name for read in batch if read is not None]
    #print(f"\nbatch_read_names:\n{batch_read_names}")
    
    # Get tags for this batch
    #print(f"\nbatch_read_names:\n{batch_read_names}")
    read_tags = extract_batch_tags(batch_read_names, fastq_path)
    #print(f"\nread_tags:\n{read_tags}")

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
    batch_size: int = 200000
):
    """
    High-performance SAM/BAM processing with improved error handling and progress tracking.
    """
    print(f"Processing SAM/BAM: {input_path}", file=sys.stderr)
    
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
    if not os.path.exists(fastq_path):
        raise FileNotFoundError(f"FASTQ file not found: {fastq_path}")
    
    os.makedirs(temp_dir, exist_ok=True)
    
    if num_processes is None:
        num_processes = max(1, mp.cpu_count() - 1)
    print(f"Using {num_processes} processes", file=sys.stderr)
    
    try:
        with pysam.AlignmentFile(input_path, 'rb', check_sq=False, require_index=False) as infile:
            print("Counting total reads...", file=sys.stderr)
            total_reads = sum(1 for _ in infile)
            header_dict = infile.header.to_dict()
        
        print(f"Total reads to process: {total_reads}", file=sys.stderr)
        processed_reads = 0
        start_time = time.time()
        successful_temp_files = []
        failed_chunks = []
        
        with mp.Pool(processes=num_processes) as pool:
            try:
                chunk_groups = list(chunk_iterator(input_path, batch_size, num_processes))
                
                for group_idx, chunk_group in enumerate(chunk_groups):
                    process_args = [
                        (
                            input_path,
                            fastq_path,
                            header_dict,
                            chunk_start,
                            chunk_end,
                            os.path.join(temp_dir, f"temp_{group_idx}_{i}.bam"),
                            batch_size
                        )
                        for i, (chunk_start, chunk_end) in enumerate(chunk_group)
                    ]
                    
                    # Process chunks and collect results
                    results = []
                    for arg in process_args:
                        result = process_chunk(arg)
                        results.append(result)
                        
                        # Handle result immediately
                        success, result_data = result
                        if success:
                            successful_temp_files.append(result_data)
                        else:
                            failed_chunks.append((arg, result_data))
                            print(f"Warning: Chunk processing failed: {result_data}", file=sys.stderr)
                    
                    # Update progress
                    processed_reads += batch_size * len(chunk_group)
                    elapsed_time = time.time() - start_time
                    reads_per_sec = processed_reads / elapsed_time if elapsed_time > 0 else 0
                    
                    print(f"Processed approximately {processed_reads}/{total_reads} reads "
                          f"({processed_reads/total_reads*100:.2f}%) "
                          f"| Speed: {reads_per_sec:.2f} reads/sec", file=sys.stderr)
                
                if failed_chunks:
                    raise RuntimeError(f"Some chunks failed processing: {failed_chunks}")
                
                # Merge successful files
                print("Merging output files...", file=sys.stderr)
                with pysam.AlignmentFile(output_path, 'wb', header=header_dict) as final_out:
                    for temp_file in successful_temp_files:
                        if os.path.exists(temp_file) and os.path.getsize(temp_file) > 0:
                            with pysam.AlignmentFile(temp_file, 'rb') as temp_in:
                                for read in temp_in:
                                    final_out.write(read)
                        else:
                            print(f"Warning: Empty or missing temp file: {temp_file}", file=sys.stderr)
                
            except Exception as e:
                print(f"Error in main processing loop: {str(e)}", file=sys.stderr)
                raise
            
            finally:
                # Clean up temp files
                for temp_file in successful_temp_files:
                    if os.path.exists(temp_file):
                        try:
                            os.remove(temp_file)
                        except Exception as e:
                            print(f"Warning: Failed to remove temp file {temp_file}: {e}", 
                                  file=sys.stderr)
        
        elapsed_time = time.time() - start_time
        print(f"Processing complete. Total time: {elapsed_time:.2f} sec "
              f"| Average speed: {total_reads/elapsed_time:.2f} reads/sec", file=sys.stderr)
    
    except Exception as e:
        print(f"Fatal error in processing: {str(e)}", file=sys.stderr)
        raise