# -*- coding: utf-8 -*-
"""
#!/usr/bin/env python3
@author: David Wragg
"""

import gzip
import sys
from collections import defaultdict
from argparse import RawTextHelpFormatter
import logging
from typing import Dict
from scarecrow import __version__
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string

def parser_rake(parser):
    subparser = parser.add_parser(
        "rake",
        description="""
Rake valid barcodes from scarecrow SAM file and inject into matching reads in source FASTQ file.

Example:

scarecrow rake --barcodes whitelist.txt --max_mismatch 3 --out barcode_mismatches.json
---
""",
        help="Rake valid barcodes from scarecrow SAM file and inject into matching reads in target FASTQ file",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument(
        "-s",
        "--sam",
        metavar="<file>",
        help=("Input scarecrow SAM file"),
        type=str,
        required=True,
        default=None,
    )
    subparser.add_argument(
        "-f",
        "--fastq",
        metavar="<file>",
        help=("Input source FASTQ file"),
        type=str,
        required=True,
        default=None,
    )
    subparser.add_argument(
        "-o",
        "--out",
        metavar="<file>",
        help=("Output corrected FASTQ file"),
        type=str,
        required=True,
        default=None,
    )
    subparser.add_argument(
        "-p",
        "--position",
        metavar="<int>",
        help=("Position for barcode sequence injection [1]"),
        type=int,
        default=1,
    )
    subparser.add_argument(
        "-l",
        "--length",
        metavar="<int>",
        help=("Expected barcode length [16]"),
        type=int,
        default=16,
    )
    return subparser


def validate_rake_args(parser, args):
    """
    Validate arguments
    """
    # Global logger setup
    logfile = f"./scarecrow_rake_{generate_random_string()}.log"
    logger = setup_logger(logfile)
    logger.info(f"scarecrow version {__version__}")
    logger.info(f"logfile: '{logfile}'")

    if args.position < 1:
        sys.exit("ERROR: --insert-pos must be >= 1")

    run_rake(
        sam_file=args.sam, fastq_file=args.fastq, output_file=args.out, position=args.position, length=args.length,
    )


@log_errors
def run_rake(
    sam_file: str = None, fastq_file: str = None, output_file: str = None, position: int = 1, length: int = 16
) -> None:
    """
    Function to rake barcodes from SAM file and inject into FASTQ file
    """
    logger = logging.getLogger("scarecrow")

    sys.stderr.write("[INFO] Extracting valid CB tags from SAM...\n")
    cb_map = extract_valid_cb_tags(sam_file, length)

    if not cb_map:
        sys.stderr.write("[WARNING] No valid CB tags found\n")

    sys.stderr.write("[INFO] Modifying FASTQ...\n")
    modify_fastq(
        fastq_gz_in=fastq_file,
        fastq_gz_out=output_file,
        cb_map=cb_map,
        insert_pos=position,
    )


    logger.info("Finished!")


def extract_valid_cb_tags(sam_path, cb_length):
    """
    Returns dict: read_id -> CB sequence
    """
    cb_map = {}
    total = 0
    kept = 0

    with open(sam_path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue

            total += 1
            fields = line.rstrip().split("\t")
            qname = fields[0]
            tags = fields[11:]

            cb = None
            for tag in tags:
                if tag.startswith("CB:Z:"):
                    cb = tag[5:]
                    break

            if cb is None:
                continue

            if len(cb) != cb_length:
                continue

            if set(cb) == {"N"}:
                continue

            cb_map[qname] = cb
            kept += 1

    sys.stderr.write(
        f"[INFO] SAM reads processed: {total}\n"
        f"[INFO] Valid CB reads kept: {kept}\n"
    )

    return cb_map

def modify_fastq(
    fastq_gz_in,
    fastq_gz_out,
    cb_map,
    insert_pos,
    progress_interval = 1000000,
):
    insert_idx = insert_pos - 1

    total = 0
    modified = 0

    with gzip.open(fastq_gz_in, "rt") as fin, \
         gzip.open(fastq_gz_out, "wt") as fout:

        while True:
            header = fin.readline()
            if not header:
                break

            seq = fin.readline().rstrip()
            plus = fin.readline()
            qual = fin.readline().rstrip()

            total += 1

            read_id = header.split()[0][1:]

            if read_id in cb_map:
                cb = cb_map[read_id]
                k = len(cb)

                if insert_idx + k <= len(seq):
                    seq = seq[:insert_idx] + cb + seq[insert_idx + k:]
                    modified += 1

            fout.write(header)
            fout.write(seq + "\n")
            fout.write(plus)
            fout.write(qual + "\n")

            if total % progress_interval == 0:
                sys.stderr.write(
                    f"[INFO] FASTQ reads processed: {total:,} | "
                    f"modified: {modified:,}\n"
                )

    sys.stderr.write(
        f"[DONE] FASTQ reads processed: {total:,}\n"
        f"[DONE] Reads modified: {modified:,}\n"
    )
