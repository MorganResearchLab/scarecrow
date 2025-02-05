#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
Refactored for on-demand FASTQ tag extraction and multiprocessing
"""

import pysam
import multiprocessing as mp
import os
import heapq
import bisect
import sqlite3
from typing import Optional
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
        fastq_file = args.fastq, 
        bam_file = args.sam, 
        out_file = args.out, 
        threads = args.threads
    )

def split_bam_file(bam_file, num_chunks):
    """Split the BAM file into chunks, distributing reads evenly."""
    bam_chunks = [f"chunk_{i}.bam" for i in range(num_chunks)]
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Count the total number of reads
        total_reads = sum(1 for _ in bam)
        reads_per_chunk = total_reads // num_chunks

        # Reset the BAM file iterator
        bam.reset()

        # Open all chunk files for writing
        chunk_files = [pysam.AlignmentFile(chunk, "wb", template=bam) for chunk in bam_chunks]

        # Distribute reads across chunks
        for i, read in enumerate(bam):
            chunk_idx = min(i // reads_per_chunk, num_chunks - 1)
            chunk_files[chunk_idx].write(read)

        # Close all chunk files
        for chunk_file in chunk_files:
            chunk_file.close()

    return bam_chunks

@log_errors
def run_samtag(
    fastq_file: str = None, 
    bam_file: str = None, 
    out_file: str = None, 
    threads: Optional[int] = None
) -> None:
    """
    Multiprocessing function to process SAM tags efficiently
    """
    # Setup logging
    logfile = f'./scarecrow_samtag_{generate_random_string()}.log'
    logger = setup_logger(logfile)
    logger.info(f"logfile: '{logfile}'")

    # Create or load the FASTQ index
    index_db = f'{fastq_file}.db'
    if not os.path.exists(index_db):
        logger.info("Creating FASTQ index...")
        create_fastq_index(fastq_file, index_db)
        logger.info("FASTQ index created.")
    else:
        logger.info(f"Using existing FASTQ index db: '{index_db}'")

    # Connect to the SQLite database and output first record
    conn = sqlite3.connect(index_db)
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM fastq_index LIMIT 1")
    result = cursor.fetchone()
    if result:
        logger.info(f"First record in db: ${result}")
        logger.info(f"result[0]: ${result[0]}")
        tags = get_tags_from_fastq(fastq_file, result[0])
        for tag, value in tags.items():
            logger.info(f"{tag} : {value}")
    else:
        logger.info("No records found.")
    conn.close()

    # Split the BAM file into chunks
    print("Splitting BAM file into chunks...")
    bam_chunks = split_bam_file(bam_file, threads)

    # Process each chunk in parallel
    pool = mp.Pool(processes = threads)
    results = []
    for i, chunk in enumerate(bam_chunks):
        output_sam = f"chunk_{i}.sam"
        results.append(pool.apply_async(process_chunk, args=(chunk, fastq_file, index_db, output_sam)))
    pool.close()
    pool.join()

    # Combine the processed SAM files into a single BAM file
    with pysam.AlignmentFile(out_file, "wb", template=pysam.AlignmentFile(bam_file, "rb")) as out_bam:
        for i in range(threads):
            output_sam = f"chunk_{i}.sam"
            with open(output_sam, "r") as in_sam:
                for line in in_sam:
                    out_bam.write(pysam.AlignedSegment.fromstring(line, out_bam.header))
    
    # Clean up temporary files
    for i in range(threads):
        os.remove(f"chunk_{i}.bam")
        os.remove(f"chunk_{i}.sam")
        

def parse_fastq_header(header):
    """Parse the FASTQ header to extract tags."""
    tags = {}
    for tag in header.split()[1:]:
        key, value = tag.split('=', 1)
        tags[key] = value
    return tags

def create_fastq_index(fastq_file, index_db):
    """Create an SQLite database mapping read names to file offsets."""
    # Connect to the SQLite database
    conn = sqlite3.connect(index_db)
    cursor = conn.cursor()
    # Create the index table
    cursor.execute("CREATE TABLE IF NOT EXISTS fastq_index (read_name TEXT PRIMARY KEY, offset INTEGER)")
    conn.commit()

    # Read the FASTQ file and populate the database
    offset = 0
    with open(fastq_file, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                break
            read_name = header.split()[0][1:]  # Remove '@' and take the first part
            # Insert the read name and offset into the database
            cursor.execute("INSERT INTO fastq_index (read_name, offset) VALUES (?, ?)", (read_name, offset))
            offset = f.tell()  # Save the current file offset
            # Skip the next 3 lines (sequence, '+', quality)
            f.readline()
            f.readline()
            f.readline()
    conn.commit()
    conn.close()

def load_fastq_index(index_file):
    """Load the FASTQ index into memory."""
    index = []
    with open(index_file, "r") as idx:
        for line in idx:
            read_name, offset = line.strip().split("\t")
            index.append((read_name, int(offset)))
    return index

def get_tags_from_fastq(fastq_file, offset):
    """Retrieve tags from the FASTQ file at a specific offset."""
    with open(fastq_file, "r") as f:
        f.seek(offset)
        header = f.readline().strip()
        return parse_fastq_header(header)

def process_chunk(bam_chunk, fastq_file, index_db, output_sam):
    """Process a chunk of the BAM file and add tags from the FASTQ file."""
    # Connect to the SQLite database
    conn = sqlite3.connect(index_db)
    cursor = conn.cursor()

    with pysam.AlignmentFile(bam_chunk, "rb") as bam, open(output_sam, "w") as out_sam:
        for read in bam:
            read_name = read.query_name
            # Query the database for the read's offset
            cursor.execute("SELECT offset FROM fastq_index WHERE read_name = ?", (read_name,))
            result = cursor.fetchone()
            if result:
                offset = result[0]
                tags = get_tags_from_fastq(fastq_file, offset)
                for tag, value in tags.items():
                    read.set_tag(tag, value)
            out_sam.write(read.to_string() + "\n")

    conn.close()