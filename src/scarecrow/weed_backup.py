#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
Refactored for on-demand FASTQ tag extraction and multiprocessing
"""

import pysam
import multiprocessing as mp
import os
import gzip
import shutil
import sqlite3
import logging
from typing import Optional
from pathlib import Path
from argparse import RawTextHelpFormatter
from scarecrow import __version__
from scarecrow.logger import log_errors, setup_logger
from scarecrow.reap import setup_worker_logger, parse_seed_arguments, BarcodeMatcherOptimized
from scarecrow.tools import generate_random_string

def parser_weed(parser):
    """
    Add samtag subparser with command-line arguments
    """
    subparser = parser.add_parser(
        "weed",
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
        help=("Path to fastq file to weed barcode from sequence header"),
        type=str,
        required=True,
    )
    subparser.add_argument(
        "-i", "--barcode_index",
        metavar="<file>",
        help=("Barcode index to extract (i.e. first sequence in header is -i 1) [1]"),
        type=int,
        default=1
    )
    subparser.add_argument(
        "-c", "--barcodes",
        metavar="<string>",
        type=str,
        required=True,
        help='A single barcode whitelist file in format <barcode_name>:<whitelist_name>:<whitelist_file>\n\t(e.g. BC1:v1:barcodes1.txt)',
    )
    subparser.add_argument(
        "-s", "--sam",
        metavar="<file>",
        help=("Path to SAM file to update barcode tags"),
        type=str,
        required=True,
    )
    subparser.add_argument(
        "-o", "--out",
        metavar="<file>",
        help=("Path to SAM file to output"),
        type=str,
        default="samtag.sam",
    )
    subparser.add_argument(
        "-m", "--mismatches",
        metavar="<int>",
        type=int,
        default=1,
        help='Number of allowed mismatches in barcode [1]',
    )
    subparser.add_argument(
        "-q", "--base_quality",
        metavar="<int>",
        type=int,
        default=None,
        help='Minimum base quality filter [None]',
    )
    subparser.add_argument(
        "-v", "--verbose",
        action='store_true',
        help='Enable verbose output [false]'
    )
    subparser.add_argument(
        "-@", "--threads",
        metavar="<int>",
        help=("Number of processing threads [1]"),
        type=int,
        default=1,
    )
    
    return subparser

def validate_weed_args(parser, args):
    """
    Validate and run samtag processing
    """
        # Get outfile dirname for writing temp chunk files to
    outpath = os.path.dirname(args.out)
    if not outpath:
        outpath = "./"
    else:
        outpath = f"{outpath}/"
    

    # Setup logging
    rnd_string = generate_random_string()
    logfile = f'{outpath}scarecrow_weed_{rnd_string}.log'
    logger = setup_logger(logfile)
    logger.info(f"scarecrow version {__version__}")
    logger.info(f"logfile: '{logfile}'")
    logger.info(f"outpath: '{outpath}'")

    run_weed(
        fastq_file = args.fastq, 
        barcode_index = args.barcode_index,
        barcodes = args.barcodes,
        bam_file = args.sam, 
        out_file = args.out, 
        outpath = outpath,
        rnd_string = rnd_string,
        mismatches = args.mismatches,
        base_quality = args.base_quality,
        verbose = args.verbose,
        threads = args.threads
    )

def count_reads_no_header(bam_file):
    """Count the number of reads in a BAM file without a header."""
    read_count = 0
    with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam:
        for _ in bam.fetch(until_eof=True):  
            read_count += 1
    return read_count

def split_bam_file(bam_file, num_chunks, out):
    """Split the BAM file into chunks, distributing reads evenly."""
    bam_chunks = [f"{out}_chunk_{i}.bam" for i in range(num_chunks)]
    with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam:
        # Count the total number of reads in the SAM file
        with open(bam_file, "r") as infile:
            total_reads = sum(1 for line in infile if not line.startswith("@"))  # Skip header lines if any

        reads_per_chunk = total_reads // num_chunks

        # Reset the BAM file iterator
        bam.reset()

        # Open all chunk files for writing
        chunk_files = [pysam.AlignmentFile(chunk, "wb", template=bam) for chunk in bam_chunks]

        # Distribute reads across chunks
        for i, read in enumerate(bam.fetch(until_eof=True)):
            chunk_idx = min(i // reads_per_chunk, num_chunks - 1)
            chunk_files[chunk_idx].write(read)

        # Close all chunk files
        for chunk_file in chunk_files:
            chunk_file.close()

    return bam_chunks


@log_errors
def run_weed(
    fastq_file: str = None, 
    barcode_index: int = 1,
    barcodes: str = None,
    bam_file: str = None, 
    out_file: str = None, 
    outpath: str = None,
    rnd_string: str = None,
    mismatches: int = 1,
    base_quality: int = 10,
    verbose: bool = False,
    threads: Optional[int] = None
) -> None:
    """
    Multiprocessing function to process SAM tags efficiently
    """
    logger = logging.getLogger('scarecrow')

    # Check if FASTQ data is compressed, and if so uncompress it
    remove_fq = False
    if is_gz_file(fastq_file):
        logger.info(f"Decompressing '{fastq_file}' for processing.")
        fastq_file = decompress_gz_file(fastq_file)
        remove_fq = True

    # Create or load the FASTQ index
    index_db = f'{fastq_file}.db'
    if not os.path.exists(index_db):
        logger.info("Creating FASTQ index")
        create_fastq_index(fastq_file, index_db)
        logger.info("FASTQ index created.")
    else:
        logger.info(f"Using existing FASTQ index db: '{index_db}'")

    # Split the BAM file into chunks
    logger.info("Splitting BAM file into chunks")
    bam_chunks = split_bam_file(bam_file, threads, f'{outpath}{rnd_string}')
    logger.info(f"{bam_chunks}")

    # Extract barcodes and convert whitelist to set
    barcode_files = parse_seed_arguments([barcodes])
    for whitelist, filename in barcode_files.items():
        pass
    
    # Create matcher
    matcher = BarcodeMatcherOptimized(
        barcode_files = barcode_files,
        mismatches = mismatches,
        base_quality_threshold = base_quality,
        verbose = verbose
    )

    # Process each chunk in parallel
    logger.info("Processing BAM chunks")
    pool = mp.Pool(processes = threads)
    results = []
    for i, chunk in enumerate(bam_chunks):
        output_sam = f"{outpath}{rnd_string}_chunk_{i}.sam"
        results.append(pool.apply_async(process_chunk, args=(chunk, fastq_file, index_db, output_sam, barcode_index, matcher, whitelist)))
    pool.close()
    pool.join()

    # Combine the processed SAM files into a single BAM file
    logger.info(f"Writing SAM chunks to: '{out_file}'")
    with pysam.AlignmentFile(out_file, "w", template=pysam.AlignmentFile(bam_file, "rb", check_sq=False)) as out_bam:
        for i in range(threads):
            output_sam = f"{outpath}{rnd_string}_chunk_{i}.sam"
            with open(output_sam, "r") as in_sam:
                for read in in_sam:
                    out_bam.write(pysam.AlignedSegment.fromstring(read.strip(), out_bam.header))
    
    # Report disk space used
    bam_size = sum(f.stat().st_size for f in Path(".").glob(f"{outpath}/{rnd_string}_chunk*.bam"))
    sam_size = sum(f.stat().st_size for f in Path(".").glob(f"{outpath}/{rnd_string}_chunk*.sam"))
    logger.info(f"Total BAM chunk disk space used: {bam_size}")
    logger.info(f"Total SAM chunk disk space used: {sam_size}")
    
    # Clean up temporary files
    for i in range(threads):
        if os.path.exists(f"{outpath}{rnd_string}_chunk_{i}.bam"):
            os.remove(f"{outpath}{rnd_string}_chunk_{i}.bam")
        if os.path.exists(f"{outpath}{rnd_string}_chunk_{i}.sam"):
            os.remove(f"{outpath}{rnd_string}_chunk_{i}.sam")
    
    if remove_fq is True:
        logger.info(f"Removing decompressed FASTQ file that was generated: '{fastq_file}'")
        os.remove(fastq_file)
    logger.info(f"Removing FASTQ index that was generated: '{index_db}'")
    os.remove(index_db)

    logger.info("Finished!")
        
def is_gz_file(filepath):
    """Check if a file is a gzip-compressed file."""
    try:
        with gzip.open(filepath, 'rb') as f:
            f.read(1)  # Try to read a single byte
        return True
    except (OSError, gzip.BadGzipFile):
        return False

def decompress_gz_file(filepath):
    """Decompress a .gz file without removing the original."""
    if not is_gz_file(filepath):
        raise ValueError(f"{filepath} is not a valid gzip file.")
    
    decompressed_filepath = filepath.rstrip('.gz')
    with gzip.open(filepath, 'rb') as f_in:
        with open(decompressed_filepath, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return decompressed_filepath

def parse_fastq_header(header):
    """Parse the FASTQ header to extract tags."""
    tags = {}
    for tag in header.split()[1:]:
        key, value = tag.split('=', 1)
        tags[key] = value
    return tags

@log_errors
def create_fastq_index(fastq_file, index_db):
    """Create an SQLite database mapping read names to file offsets."""
    logger = logging.getLogger('scarecrow')

    # Connect to the SQLite database
    conn = sqlite3.connect(index_db)
    cursor = conn.cursor()
    # Create the index table
    cursor.execute("CREATE TABLE IF NOT EXISTS fastq_index (read_name TEXT PRIMARY KEY, offset INTEGER)")
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
            cursor.execute("INSERT INTO fastq_index (read_name, offset) VALUES (?, ?)", (read_name, offset))
            # Skip the next 3 lines (sequence, '+', quality)
            f.readline()
            f.readline()
            f.readline()
            #logger.info(f"{read_name}, {offset}")
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

def get_barcode_from_fastq(fastq_file: str, offset: int, barcode_index: int) -> Optional[str]:
    """
    Retrieve the barcode sequence from the FASTQ file at a given offset.
    The barcode is extracted from the read description section, which contains colon-delimited fields.
    The barcode sequences are identified by the presence of a '+' symbol in the description.
    """
    with open(fastq_file, "r") as fastq:
        # Seek to the specific read using the offset
        fastq.seek(offset)
        header = fastq.readline().strip()
        sequence = fastq.readline().strip()
        plus_line = fastq.readline().strip()
        quality = fastq.readline().strip()

        # Extract the description section
        description = header.split(" ", 1)[1] if " " in header else ""

        # Split the description into colon-delimited fields
        fields = description.split(':')

        # Look for the field containing the '+' symbol (barcode section)
        barcode_field = None
        for field in fields:
            if '+' in field:
                barcode_field = field
                break

        if barcode_field:
            # Split the barcode field on '+' to get individual barcodes
            barcodes = barcode_field.split('+')
            if len(barcodes) > barcode_index:
                return barcodes[barcode_index - 1] # index is 0-based, user input is 1-based

    return None

@log_errors        
def process_chunk(bam_chunk, fastq_file, index_db, output_sam, barcode_index, matcher, whitelist):
    """Process a chunk of the BAM file and add tags from the FASTQ file."""
    logger = setup_worker_logger()

    # Connect to the SQLite database
    conn = sqlite3.connect(index_db)
    cursor = conn.cursor()

    with pysam.AlignmentFile(bam_chunk, "rb", check_sq=False) as bam, open(output_sam, "w") as out_sam:
        for read in bam.fetch(until_eof=True):
            #logger.info(f"{read}")
            read_name = read.query_name
            if not read_name:
                logger.warning(f"Skipping read with no name in {bam_chunk}")
                continue

            # Query the database for the read's offset
            cursor.execute("SELECT offset FROM fastq_index WHERE read_name = ?", (read_name,))
            result = cursor.fetchone()
            if result:
                offset = result[0]
                barcode = get_barcode_from_fastq(fastq_file, offset, barcode_index)

                # Search for barcode in forward orientation
                matched_barcode, mismatch_count, adj_position = matcher.find_match(
                    barcode, None, whitelist, 'forward', 1, len(barcode), 0, None)
                # If none found, search in reverse orientation
                if 'NN' in matched_barcode:
                    matched_barcode, mismatch_count, adj_position = matcher.find_match(
                        barcode, None, whitelist, 'reverse', 1, len(barcode), 0, None)

                # Update read tags
                if barcode:
                    # Update CR and CY tags
                    cr_tag = read.get_tag('CR') if read.has_tag('CR') else ''
                    cy_tag = read.get_tag('CY') if read.has_tag('CY') else ''
                    cb_tag = read.get_tag('CB') if read.has_tag('CB') else ''
                    xm_tag = read.get_tag('XM') if read.has_tag('XM') else ''
                    xp_tag = read.get_tag('XP') if read.has_tag('XP') else ''
                    read.set_tag('CR', f"{cr_tag}_{barcode}")
                    read.set_tag('CY', f"{cy_tag}_{'!' * len(barcode)}")
                    if matched_barcode:
                        read.set_tag('CB', f"{cb_tag}_{matched_barcode}")
                        read.set_tag('XM', f"{xm_tag}_{mismatch_count}")
                        read.set_tag('XP', f"{xp_tag}_0")
            
                out_sam.write(read.to_string() + "\n")

    conn.close()

