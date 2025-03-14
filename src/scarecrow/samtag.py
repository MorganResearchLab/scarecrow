#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
Refactored for on-demand FASTQ tag extraction and multiprocessing
"""

import pysam
import multiprocessing as mp
import os
import sqlite3
import logging
from typing import Optional
from pathlib import Path
from argparse import RawTextHelpFormatter
from scarecrow import __version__
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
        "-f",
        "--fastq",
        metavar="<file>",
        help=("Path to scarecrow reap fastq file"),
        type=str,
        required=True,
    )
    subparser.add_argument(
        "-s",
        "--sam",
        metavar="<file>",
        help=("Path to SAM file to update tags"),
        type=str,
        required=True,
    )
    subparser.add_argument(
        "-o",
        "--out",
        metavar="<file>",
        help=("Path to SAM file to output"),
        type=str,
        default="samtag.sam",
    )
    subparser.add_argument(
        "-@",
        "--threads",
        metavar="<int>",
        help=("Number of processing threads [1]"),
        type=int,
        default=1,
    )

    return subparser


def validate_samtag_args(parser, args):
    """
    Validate arguments
    """

    # Global logger setup
    logfile = f"./scarecrow_seed_{generate_random_string()}.log"
    logger = setup_logger(logfile)
    logger.info(f"scarecrow version {__version__}")
    logger.info(f"logfile: '{logfile}'")

    run_samtag(
        fastq_file=args.fastq,
        bam_file=args.sam,
        out_file=args.out,
        threads=args.threads,
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
        chunk_files = [
            pysam.AlignmentFile(chunk, "wb", template=bam) for chunk in bam_chunks
        ]

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
    threads: Optional[int] = None,
) -> None:
    """
    Multiprocessing function to process SAM tags efficiently
    """
    logger = logging.getLogger("scarecrow")

    # Create or load the FASTQ index
    index_db = f"{fastq_file}.db"
    if not os.path.exists(index_db):
        logger.info("Creating FASTQ index...")
        create_fastq_index(fastq_file, index_db)
        logger.info("FASTQ index created.")
    else:
        logger.info(f"Using existing FASTQ index db: '{index_db}'")

    # Split the BAM file into chunks
    logger.info("Splitting BAM file into chunks...")
    bam_chunks = split_bam_file(bam_file, threads)

    # Process each chunk in parallel
    logger.info("Processing BAM chunks...")
    pool = mp.Pool(processes=threads)
    results = []
    for i, chunk in enumerate(bam_chunks):
        output_sam = f"chunk_{i}.sam"
        results.append(
            pool.apply_async(
                process_chunk, args=(chunk, fastq_file, index_db, output_sam)
            )
        )
    pool.close()
    pool.join()

    # Combine the processed SAM files into a single BAM file
    logger.info(f"Writing SAM chunks to {out_file}...")
    with pysam.AlignmentFile(
        out_file, "wb", template=pysam.AlignmentFile(bam_file, "rb")
    ) as out_bam:
        for i in range(threads):
            output_sam = f"chunk_{i}.sam"
            with open(output_sam, "r") as in_sam:
                for line in in_sam:
                    out_bam.write(pysam.AlignedSegment.fromstring(line, out_bam.header))

    # Report disk space used
    bam_size = sum(f.stat().st_size for f in Path(".").glob("chunk*.bam"))
    sam_size = sum(f.stat().st_size for f in Path(".").glob("chunk*.sam"))
    logger.info(f"Total BAM chunk disk space used: {bam_size}")
    logger.info(f"Total SAM chunk disk space used: {sam_size}")

    # Clean up temporary files
    logger.info("Removing chunks...")
    for i in range(threads):
        os.remove(f"chunk_{i}.bam")
        os.remove(f"chunk_{i}.sam")

    logger.info("Finished!")


def parse_fastq_header(header):
    """Parse the FASTQ header to extract tags."""
    tags = {}
    for tag in header.split()[1:]:
        key, value = tag.split("=", 1)
        tags[key] = value
    return tags


def create_fastq_index(fastq_file, index_db):
    """Create an SQLite database mapping read names to file offsets."""
    # Connect to the SQLite database
    conn = sqlite3.connect(index_db)
    cursor = conn.cursor()
    # Create the index table
    cursor.execute(
        "CREATE TABLE IF NOT EXISTS fastq_index (read_name TEXT PRIMARY KEY, offset INTEGER)"
    )
    conn.commit()

    # Read the FASTQ file in binary mode and populate the database
    with open(fastq_file, "rb") as f:
        while True:
            # Record the offset **before** reading the header line
            offset = f.tell()
            header = f.readline()
            if not header:
                break  # End of file
            # Decode the header line from bytes to string
            header = header.decode("utf-8").strip()
            read_name = header.split()[0][1:]  # Remove '@' and take the first part
            # Insert the read name and offset into the database
            cursor.execute(
                "INSERT INTO fastq_index (read_name, offset) VALUES (?, ?)",
                (read_name, offset),
            )
            # Skip the next 3 lines (sequence, '+', quality)
            f.readline()
            f.readline()
            f.readline()
            # logger.info(f"{read_name}, {offset}")
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
    """Retrieve specific tags from the FASTQ file at a given offset."""
    allowed_keys = {"CR", "CY", "CB", "XP", "XM", "UR", "UY"}

    with open(fastq_file, "rb") as f:
        f.seek(offset)
        header = f.readline().decode("utf-8").strip()
        # The header starts with '@' followed by the read name and then the tags
        if header.startswith("@"):
            parts = header.split()
            tags = {
                key: value
                for part in parts[1:]
                if "=" in part and (key := part.split("=", 1)[0]) in allowed_keys
                for value in [part.split("=", 1)[1]]
            }
            return tags

    return {}


def process_chunk(bam_chunk, fastq_file, index_db, output_sam):
    """Process a chunk of the BAM file and add tags from the FASTQ file."""
    # Connect to the SQLite database
    conn = sqlite3.connect(index_db)
    cursor = conn.cursor()

    with pysam.AlignmentFile(bam_chunk, "rb") as bam, open(output_sam, "w") as out_sam:
        for read in bam:
            read_name = read.query_name
            # Query the database for the read's offset
            cursor.execute(
                "SELECT offset FROM fastq_index WHERE read_name = ?", (read_name,)
            )
            result = cursor.fetchone()
            if result:
                offset = result[0]
                tags = get_tags_from_fastq(fastq_file, offset)
                for tag, value in tags.items():
                    read.set_tag(tag, value)
            out_sam.write(read.to_string() + "\n")

    conn.close()
