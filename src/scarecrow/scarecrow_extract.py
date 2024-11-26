#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 09:37:12 2024

@author: David Wragg
"""

from seqspec import Assay
from seqspec.utils import load_spec
from seqspec.seqspec_print import run_seqspec_print
from seqspec.seqspec_index import get_index_by_primer
from argparse import RawTextHelpFormatter
from scarecrow.read_fastqs import process_paired_fastq_batches

def parser_extract(parser):
    subparser = parser.add_parser(
        "extract",
        description="""
Extract cDNA sequence from fastq files

Examples:
scarecrow extract spec.yaml fastqs -v -o ~/path/to/output
---
""",
        help="Extract cDNA from fastqs",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument("yaml", help="Sequencing specification yaml file")
    subparser.add_argument("fastqs", nargs="+", help="List of FASTQ files")
    subparser.add_argument(
        "-o",
        metavar="OUT",
        help=("Path to output cDNA fastq files"),
        type=str,
        default=None,
    )
    subparser.add_argument(
        "-p",
        help=("Paired-end sequencing"),
        action='store_true'
    )
    subparser.add_argument(
        "-v",
        help=("Verbose output"),
        action='store_true'
    )
    return subparser


def validate_extract_args(parser, args):
    run_extract(yaml = args.yaml, fastqs = [f for f in args.fastqs], 
                outdir = args.o, paired = args.p, verbose = args.v)


def run_extract(yaml, fastqs, outdir, paired, verbose):
    """
    Employs seqspec functions to (1) output library spec and (2) identify elements contained in sequencing reads.
    The identified elements are then extracted from paired-end fastq files in batches and written to file (jsonl).
    """
    
    # Run seqspec print to get format
    print(f"\033[32m\nseqspec print {yaml}:\033[0m\n")
    run_seqspec_print(yaml, fmt="library-ascii", o = None)

    # import seqspec.yaml
    spec = load_spec(yaml)

    # Run seqspec get_index_by_primer to identify library elements contained in reads
    elements = region_indices(spec, fastqs)

    # Extract elements from sequencing reads
    results = process_paired_fastq_batches(elements, batch_size=1000, max_batches=5, output_file='fastq_regions.jsonl')
    

    if outdir:
        # This needs to be edited still
        print("")
    else:
        print("")
    return spec



def region_indices(spec: Assay, fastqs):
    """
    Identify library elements contained in sequencing reads
    """
    print(f"\033[32m\nLibrary elements identified by seqspec.get_index_by_primer\033[0m")
    indices = []
    for fastq in fastqs:
        index = get_index_by_primer(spec, "rna", fastq)
        indices.append(index)
    
    for index in indices:
        for file, regions in index.items():
            if file in fastqs:
                print(f"\033[34m\n{file}\033[0m")
                for region in regions:
                    print(f"\t{region.region_id}: {region.start}-{region.stop}")        

    return indices
