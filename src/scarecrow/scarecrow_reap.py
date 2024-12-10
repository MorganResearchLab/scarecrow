#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
"""

from argparse import RawTextHelpFormatter

def parser_reap(parser):
    subparser = parser.add_parser(
        "reap",
        description="""
Extract sequence range from fastq files

Example:
---
""",
        help="Extract sequence range from fastq files",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument("fastqs", nargs="+", help="List of FASTQ files")
    subparser.add_argument(
        "-o", "--out",
        metavar="out.fastq",
        help=("Path to output fastq file"),
        type=str,
        default="./cDNA.fq",
    )
    # To do - need to pass something here for header regions 
    # either need to directly pass barcode positions, or parse a config file (i.e. seqpsec yaml)
    subparser.add_argument(
        "-r", "--header_regions",
        metavar="region_id",
        help=("List of elements to include in sequence header (e.g. )"),
        nargs="*",
        type=str,
        default=[],
    )
    group = subparser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-1", "--read1",
        metavar="read1_range",
        help=("Sequence range to extract from read1 (e.g. (0-100))"),
        type=str,
        default=None
    )
    group.add_argument(
        "-2", "--read2",
        metavar="read2_range",
        help=("Sequence range to extract from read2 (e.g. 0-100)"),
        type=str,
        default=None
    )
    subparser.add_argument(
        "-b", "--batch_size",
        metavar="batch_size",
        help=("Number of read pairs per batch to process at a time [10000]"),
        type=int,
        default=10000,
    )
    subparser.add_argument(
        "-x", "--max_batches",
        metavar="max_batches",
        help=("Maximum number of read batches to process"),
        type=int,
        default=None,
    )
    subparser.add_argument(
        "-@", "--threads",
        metavar="threads",
        help=("Number of processing threads [4]"),
        type=int,
        default=4,
    )
    subparser.add_argument(
        "-l", "--logfile",
        metavar="logfile",
        help=("File to write log to"),
        type=str,
        default="./scarecrow.log",
    )
    return subparser

def validate_reap_args(parser, args):
    run_extract(logfile = args.logfile)
