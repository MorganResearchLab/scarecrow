#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 09:37:12 2024

@author: s14dw4
"""

from seqspec.utils import load_spec
from argparse import RawTextHelpFormatter


def parser_extract(parser):
    subparser = parser.add_parser(
        "extract",
        description="""
Extract cDNA sequence from fastq files

Examples:
scarecrow extract spec.yaml -o ~/path/to/output
---
""",
        help="Extract cDNA from fastqs",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument("yaml", help="Sequencing specification yaml file")
    subparser.add_argument(
        "-o",
        metavar="OUT",
        help=("Path to output cDNA fastq files"),
        type=str,
        default=None,
    )
    return subparser

def validate_extract_args(parser, args):
    run_extract(yaml = args.yaml, outdir = args.o)

def run_extract(yaml, outdir):
    spec = load_spec(yaml)
    format(spec)
    if outdir:
        spec.to_YAML(outdir)
    else:
        print(spec.to_YAML())