# -*- coding: utf-8 -*-
"""
#!/usr/bin/env python3
@author: David Wragg
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
from scarecrow.reap import parse_seed_arguments, BarcodeMatcherOptimized
from scarecrow.tools import generate_random_string


def parser_weed(parser):
    """
    Add samtag subparser with command-line arguments
    """
    subparser = parser.add_parser(
        "weed",
        description="""
Weeds out a sequence index from a FASTQ header and appends to barcode tags (CR, CY, CB) in corresponding SAM file.

Example:
scarecrow weed --fastq in.fastq --in in.sam --barcode_index 1 --barcodes BC3:P7:BC3.txt --mismatches 1
---
""",
        help="Update SAM file with barcode tag from a FASTQ sequence header",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument(
        "-f",
        "--fastq",
        metavar="<file>",
        help=("Path to fastq file to weed barcode from sequence header"),
        type=str,
        required=True,
    )
    subparser.add_argument(
        "-x",
        "--barcode_index",
        metavar="<file>",
        help=("Barcode index to extract (i.e. first sequence in header is -i 1) [1]"),
        type=int,
        default=1,
    )
    subparser.add_argument(
        "-c",
        "--barcodes",
        metavar="<string>",
        type=str,
        required=True,
        help="A single barcode whitelist file in format <barcode_name>:<whitelist_name>:<whitelist_file>\n\t(e.g. BC1:v1:barcodes1.txt)",
    )
    subparser.add_argument(
        "-i",
        "--in",
        dest="infile",
        metavar="<file>",
        help=("Path to SAM file to update barcode tags"),
        type=str,
        required=True,
    )
    subparser.add_argument(
        "-o",
        "--out",
        metavar="<file>",
        help=("Path to SAM file to output"),
        type=str,
    )
    subparser.add_argument(
        "-m",
        "--mismatches",
        metavar="<int>",
        type=int,
        default=1,
        help="Number of allowed mismatches in barcode [1]",
    )
    subparser.add_argument(
        "-q",
        "--base_quality",
        metavar="<int>",
        type=int,
        default=None,
        help="Minimum base quality filter [None]",
    )
    subparser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose output [false]"
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
    logfile = f"{outpath}scarecrow_weed_{rnd_string}.log"
    logger = setup_logger(logfile)
    logger.info(f"scarecrow version {__version__}")
    logger.info(f"logfile: '{logfile}'")
    logger.info(f"outpath: '{outpath}'")

    run_weed(
        fastq_file=args.fastq,
        barcode_index=args.barcode_index,
        barcodes=args.barcodes,
        infile=args.infile,
        out_file=args.out,
        outpath=outpath,
        rnd_string=rnd_string,
        mismatches=args.mismatches,
        base_quality=args.base_quality,
        verbose=args.verbose,
        threads=args.threads,
        args_string=" ".join(
            f"--{k} {v}" for k, v in vars(args).items() if v is not None
        ),
    )


def count_reads_no_header(bam_file):
    """Count the number of reads in a BAM file without a header."""
    read_count = 0
    with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam:
        for _ in bam.fetch(until_eof=True):
            read_count += 1
    return read_count


def split_bam_file(bam_file, num_chunks, out, args_string):
    """Split the BAM file into chunks, distributing reads evenly."""
    bam_chunks = [f"{out}_chunk_{i}.bam" for i in range(num_chunks)]

    with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam:
        # Extract headers
        headers = bam.header.to_dict()
        # Update PG part of header
        if "PG" not in headers:
            headers["PG"] = []
        headers["PG"].append(
            {"ID": "weed", "PN": "scarecrow", "VN": __version__, "DS": args_string}
        )

        # Open all chunk files for writing
        chunk_files = []
        for chunk in bam_chunks:
            chunk_file = pysam.AlignmentFile(chunk, "wb", header=headers)
            chunk_files.append(chunk_file)

        # Distribute reads across chunks
        for i, read in enumerate(bam.fetch(until_eof=True)):
            chunk_idx = i % num_chunks  # Distribute reads evenly across chunks
            chunk_files[chunk_idx].write(read)

        # Close all chunk files
        for chunk_file in chunk_files:
            chunk_file.close()

    return bam_chunks, headers


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
    threads: Optional[int] = None,
    args_string: str = None,
) -> None:
    """
    Multiprocessing function to process SAM tags efficiently
    """
    logger = logging.getLogger("scarecrow")

    # Check if FASTQ data is compressed, and if so uncompress it
    remove_fq = False
    if is_gz_file(fastq_file):
        logger.info(f"Decompressing '{fastq_file}' for processing.")
        fastq_file = decompress_gz_file(fastq_file)
        remove_fq = True

    # Create or load the FASTQ index
    index_db = f"{fastq_file}.db"
    if not os.path.exists(index_db):
        logger.info("Creating FASTQ index")
        create_fastq_index(fastq_file, index_db)
        logger.info("FASTQ index created.")
    else:
        logger.info(f"Using existing FASTQ index db: '{index_db}'")

    # Split the BAM file into chunks
    logger.info("Splitting BAM file into chunks")
    bam_chunks, headers = split_bam_file(
        bam_file, threads, f"{outpath}{rnd_string}", args_string
    )
    logger.info(f"{bam_chunks}")

    # Extract barcodes and convert whitelist to set
    barcode_files = parse_seed_arguments([barcodes])
    for whitelist, filename in barcode_files.items():
        pass

    # Create matcher
    matcher = BarcodeMatcherOptimized(
        barcode_files=barcode_files,
        mismatches=mismatches,
        base_quality_threshold=base_quality,
        verbose=verbose,
    )

    # Process each chunk in parallel
    logger.info("Processing BAM chunks")
    with mp.Pool(processes=threads) as pool:
        # Prepare arguments for worker tasks
        args = [
            (
                chunk,
                fastq_file,
                index_db,
                f"{outpath}{rnd_string}_chunk_{i}.sam",
                barcode_index,
                matcher,
                whitelist,
            )
            for i, chunk in enumerate(bam_chunks)
        ]

        # Use imap_unordered to process chunks
        temp_files = pool.imap_unordered(process_chunk, args, chunksize=1)

        # Combine the temporary SAM files into the final output
        with pysam.AlignmentFile(out_file, "w", header=headers) as out_sam:
            for temp_file in temp_files:
                with pysam.AlignmentFile(temp_file, "r", check_sq=False) as in_sam:
                    for read in in_sam.fetch(until_eof=True):
                        out_sam.write(read)
                os.remove(temp_file)  # Clean up the temporary file

    # Report disk space used
    bam_size = sum(
        f.stat().st_size for f in Path(".").glob(f"{outpath}/{rnd_string}_chunk*.bam")
    )
    sam_size = sum(
        f.stat().st_size for f in Path(".").glob(f"{outpath}/{rnd_string}_chunk*.sam")
    )
    logger.info(f"Total BAM chunk disk space used: {bam_size}")
    logger.info(f"Total SAM chunk disk space used: {sam_size}")

    # Clean up temporary files
    for i in range(threads):
        if os.path.exists(f"{outpath}{rnd_string}_chunk_{i}.bam"):
            os.remove(f"{outpath}{rnd_string}_chunk_{i}.bam")
        if os.path.exists(f"{outpath}{rnd_string}_chunk_{i}.sam"):
            os.remove(f"{outpath}{rnd_string}_chunk_{i}.sam")

    if remove_fq is True:
        logger.info(
            f"Removing decompressed FASTQ file that was generated: '{fastq_file}'"
        )
        os.remove(fastq_file)
    logger.info(f"Removing FASTQ index that was generated: '{index_db}'")
    os.remove(index_db)

    logger.info(f"Updated SAM file written to '{out_file}'")
    logger.info("Finished!")


def is_gz_file(filepath):
    """Check if a file is a gzip-compressed file."""
    try:
        with gzip.open(filepath, "rb") as f:
            f.read(1)  # Try to read a single byte
        return True
    except (OSError, gzip.BadGzipFile):
        return False


def decompress_gz_file(filepath):
    """Decompress a .gz file without removing the original."""
    if not is_gz_file(filepath):
        raise ValueError(f"{filepath} is not a valid gzip file.")

    decompressed_filepath = filepath.rstrip(".gz")
    with gzip.open(filepath, "rb") as f_in:
        with open(decompressed_filepath, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    return decompressed_filepath


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


def get_barcode_from_fastq(
    fastq_file: str, offset: int, barcode_index: int
) -> Optional[str]:
    """
    Retrieve the barcode sequence from the FASTQ file at a given offset.
    The barcode is extracted from the read description section, which contains colon-delimited fields.
    The barcode sequences are identified by the presence of a '+' symbol in the description.
    """
    with open(fastq_file, "r") as fastq:
        # Seek to the specific read using the offset
        fastq.seek(offset)
        header = fastq.readline().strip()
        # sequence = fastq.readline().strip()
        # plus_line = fastq.readline().strip()
        # quality = fastq.readline().strip()

        # Extract the description section
        description = header.split(" ", 1)[1] if " " in header else ""

        # Split the description into colon-delimited fields
        fields = description.split(":")

        # Look for the field containing the '+' symbol (barcode section)
        barcode_field = None
        for field in fields:
            if "+" in field:
                barcode_field = field
                break

        if barcode_field:
            # Split the barcode field on '+' to get individual barcodes
            barcodes = barcode_field.split("+")
            if len(barcodes) > barcode_index:
                return barcodes[
                    barcode_index - 1
                ]  # index is 0-based, user input is 1-based

    return None


def process_chunk(args):
    """Process a chunk of the BAM file and add tags from the FASTQ file."""
    bam_chunk, fastq_file, index_db, output_sam, barcode_index, matcher, whitelist = (
        args
    )

    # Connect to the SQLite database
    conn = sqlite3.connect(index_db)
    cursor = conn.cursor()

    # Open the output SAM file for writing
    with pysam.AlignmentFile(
        output_sam, "w", template=pysam.AlignmentFile(bam_chunk, "rb", check_sq=False)
    ) as out_sam:
        with pysam.AlignmentFile(bam_chunk, "rb", check_sq=False) as bam:
            for read in bam.fetch(until_eof=True):
                read_name = read.query_name
                if not read_name:
                    continue

                # Query the database for the read's offset
                cursor.execute(
                    "SELECT offset FROM fastq_index WHERE read_name = ?", (read_name,)
                )
                result = cursor.fetchone()
                if result:
                    offset = result[0]
                    barcode = get_barcode_from_fastq(fastq_file, offset, barcode_index)

                    # Search for barcode in forward orientation
                    matched_barcode, mismatch_count, adj_position = matcher.find_match(
                        barcode, None, whitelist, "forward", 1, len(barcode), 0, None
                    )
                    # If none found, search in reverse orientation
                    if "NN" in matched_barcode:
                        matched_barcode, mismatch_count, adj_position = (
                            matcher.find_match(
                                barcode,
                                None,
                                whitelist,
                                "reverse",
                                1,
                                len(barcode),
                                0,
                                None,
                            )
                        )

                    # Update read tags
                    if barcode:
                        # Update CR and CY tags
                        cr_tag = read.get_tag("CR") if read.has_tag("CR") else ""
                        cy_tag = read.get_tag("CY") if read.has_tag("CY") else ""
                        cb_tag = read.get_tag("CB") if read.has_tag("CB") else ""
                        xm_tag = read.get_tag("XM") if read.has_tag("XM") else ""
                        xp_tag = read.get_tag("XP") if read.has_tag("XP") else ""
                        read.set_tag("CR", f"{cr_tag},{barcode}")
                        read.set_tag("CY", f"{cy_tag},{'!' * len(barcode)}")
                        if matched_barcode:
                            read.set_tag("CB", f"{cb_tag},{matched_barcode}")
                            read.set_tag("XM", f"{xm_tag},{mismatch_count}")
                            read.set_tag("XP", f"{xp_tag},0")

                out_sam.write(read)  # Write the modified read to the output SAM file

    conn.close()
    return output_sam  # Return the path to the temporary SAM file
