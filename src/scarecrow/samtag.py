#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
Refactored for on-demand FASTQ tag extraction and multiprocessing
"""

import pysam
import multiprocessing as mp
import os
import bisect
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
        fastq_file = args.fastq, 
        bam_file = args.sam, 
        out_file = args.out, 
        batch_size = args.batch_size,
        threads = args.threads
    )

@log_errors
def run_samtag(
    fastq_file: str = None, 
    bam_file: str = None, 
    out_file: str = None, 
    batch_size: int = 200000,
    threads: Optional[int] = None
) -> None:
    """
    Multiprocessing function to process SAM tags efficiently
    """
    # Setup logging
    #logfile = f'./scarecrow_samtag_{generate_random_string()}.log'
    #logger = setup_logger(logfile)
    #logger.info(f"logfile: '{logfile}'")

    # Create or load the FASTQ index
    index_file = f'{fastq_file}.idx'
    if not os.path.exists(index_file):
        print("Creating FASTQ index...")
        create_fastq_index(fastq_file, index_file)
    print("Loading FASTQ index...")
    fastq_index = load_fastq_index(index_file)

    # Split the BAM file into chunks
    chunk_size = os.path.getsize(bam_file) // threads
    bam_chunks = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for i in range(threads):
            chunk_file = f"chunk_{i}.bam"
            with pysam.AlignmentFile(chunk_file, "wb", template=bam) as chunk:
                for j, read in enumerate(bam):
                    if j >= chunk_size:
                        break
                    chunk.write(read)
            bam_chunks.append(chunk_file)

    # Process each chunk in parallel
    pool = mp.Pool(processes=threads)
    results = []
    for i, chunk in enumerate(bam_chunks):
        output_sam = f"chunk_{i}.sam"
        results.append(pool.apply_async(process_chunk, args=(chunk, fastq_file, fastq_index, output_sam)))
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

def create_fastq_index(fastq_file, index_file):
    """Create an index file mapping read names to file offsets."""
    index = []
    with open(fastq_file, "r") as f:
        while True:
            offset = f.tell()
            header = f.readline()
            if not header:
                break
            read_name = header.split()[0][1:]  # Remove '@' and take the first part
            index.append((read_name, offset))
            # Skip the next 3 lines (sequence, '+', quality)
            f.readline()
            f.readline()
            f.readline()
    # Sort the index by read name for efficient lookup
    index.sort(key=lambda x: x[0])
    # Write the index to a file
    with open(index_file, "w") as idx:
        for read_name, offset in index:
            idx.write(f"{read_name}\t{offset}\n")

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

def process_chunk(bam_chunk, fastq_file, fastq_index, output_sam):
    """Process a chunk of the BAM file and add tags from the FASTQ file."""
    with pysam.AlignmentFile(bam_chunk, "rb") as bam, open(output_sam, "w") as out_sam:
        for read in bam:
            read_name = read.query_name
            # Perform a binary search on the FASTQ index
            idx = bisect.bisect_left(fastq_index, (read_name, 0))
            if idx < len(fastq_index) and fastq_index[idx][0] == read_name:
                offset = fastq_index[idx][1]
                tags = get_tags_from_fastq(fastq_file, offset)
                for tag, value in tags.items():
                    read.set_tag(tag, value)
            out_sam.write(read.to_string() + "\n")

