#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
"""

from seqspec import Assay
from seqspec.utils import load_spec
from seqspec.seqspec_print import run_seqspec_print
from seqspec.seqspec_index import get_index_by_primer, format_kallisto_bus
from argparse import RawTextHelpFormatter
from scarecrow.tools import process_paired_fastq_batches
from scarecrow.fastq_logging import logger, log_errors, setup_logger

def parser_extract(parser):
    subparser = parser.add_parser(
        "extract",
        description="""
Extract cDNA sequence from fastq files

Example extracting sequence elements using regions from spec.yaml:
scarecrow extract spec.yaml R1.fastq.gz R2.fastq.gz -o ~/path/to/output -r UMI Round_1_BC Round_2_BC Round_3_BC

Example identifying barcode elements using whitelists to help with debugging (results recorded to log file):
scarecrow extract spec.yaml R1.fastq.gz R2.fastq.gz --barcodes  BC1:/Users/s14dw4/Documents/Repos/scarecrow/specs/evercode/BC1.txt BC2:/Users/s14dw4/Documents/Repos/scarecrow/specs/evercode/BC2-3.txt BC3:/Users/s14dw4/Documents/Repos/scarecrow/specs/evercode/BC2-3.txt
---
""",
        help="Extract cDNA from fastqs",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument("yaml", help="Sequencing specification yaml file")
    subparser.add_argument("fastqs", nargs="+", help="List of FASTQ files")
    subparser.add_argument(
        "-o",
        metavar="out",
        help=("Path to output cDNA fastq file"),
        type=str,
        default="./cDNA.fq",
    )
    subparser.add_argument(
        "-r",
        metavar="region_id",
        help=("List of regions for cell barcode"),
        nargs="*",
        type=str,
        default=[],
    )
    subparser.add_argument(
        "-b",
        metavar="batches",
        help=("Number of read batches to process at a time before writing to file [1000]"),
        type=int,
        default=1000,
    )
    subparser.add_argument(
        "-m",
        metavar="max_batches",
        help=("Maximum number of read batches to process"),
        type=int,
        default=None,
    )
    subparser.add_argument(
        "-t",
        metavar="threads",
        help=("Number of processing threads [4]"),
        type=int,
        default=4,
    )
    subparser.add_argument(
        "--barcodes",
        metavar="barcodes",
        nargs='+', 
        help='Barcode files in format BC1:path/to/barcodes1.txt BC2:path/to/barcodes2.txt',
    )
    subparser.add_argument(
        "-l",
        metavar="logfile",
        help=("File to write log to"),
        type=str,
        default="./scarecrow.log",
    )
    return subparser

def validate_extract_args(parser, args):
    run_extract(yaml = args.yaml, fastqs = [f for f in args.fastqs], 
                output_file = args.o, batches = args.b, regions = args.r, logfile = args.l,
                threads = args.t, max_batches = args.m, barcodes = args.barcodes)

def run_extract(yaml, fastqs, output_file, batches, max_batches, regions, threads, barcodes):
    """
    Employs seqspec functions to (1) output library spec and (2) identify elements contained in sequencing reads.
    The identified elements are then extracted from paired-end fastq files in batches and written to file.
    """
    
    # Global logger setup
    logger = setup_logger(logfile)

    # Run seqspec print to get format
    print(f"\033[32m\nseqspec print \033[34m{yaml}\033[0m\n")
    run_seqspec_print(yaml, fmt="library-ascii", o = None)

    # import seqspec.yaml
    spec = load_spec(yaml)

    # Run seqspec get_index_by_primer to identify library elements contained in reads
    elements = region_indices(spec, fastqs)

    # Open file for writing output
    if output_file:
        f = open(f"{output_file}", 'w')

    # Extract elements from sequencing reads
    process_paired_fastq_batches(elements, batch_size = batches, max_batches = max_batches,
                                 num_workers = threads, region_ids = regions, output_handler = f,
                                 barcodes = None)

    # Return kallisto bus -x string, equivalent to:
    # seqspec index -t kb -m rna -r {R1.fastq.gz},{R2.fastq.gz} {yaml}
    x = format_kallisto_bus(elements)
    print(f"\033[32m\nkallisto bus -x \033[34m{x}\033[0m\n")

    return 

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


