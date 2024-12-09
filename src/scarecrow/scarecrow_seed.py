#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
"""

from seqspec.utils import load_spec
from scarecrow.fastq_logging import log_errors, setup_logger, logger
from argparse import RawTextHelpFormatter
from scarecrow.tools import process_paired_fastq_batches, region_indices

def parser_seed(parser):
    subparser = parser.add_parser(
        "seed",
        description="""
Search fastq reads for barcodes in whitelists

Example:

scarecrow seed spec.yaml R1.fastq.gz R2.fastq.gz\n\t--barcodes BC1:BC1.txt BC2:BC2.txt BC3:BC3.txt\n\t--out barcode_counts.csv 
---
""",
        help="Search fastq reads for barcodes",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument("yaml", help="Sequencing specification yaml file")
    subparser.add_argument("fastqs", nargs="+", help="List of FASTQ files")
    subparser.add_argument(
        "-c", "--barcodes",
        metavar="barcodes",
        nargs='+', 
        help='Barcode files in format BC1:barcodes1.txt BC2:barcodes2.txt',
    )
    subparser.add_argument(
        "-o", "--out",
        metavar="out",
        help=("CSV file to write barcode counts to"),
        type=str,
        default="./barcode_counts.csv",
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

def validate_seed_args(parser, args):
    run_seed(yaml = args.yaml, 
             fastqs = [f for f in args.fastqs],
             barcodes = args.barcodes, 
             output_file = args.out, 
             logfile = args.logfile,
             batches = args.batch_size, 
             threads = args.threads, 
             max_batches = args.max_batches)

@log_errors
def run_seed(yaml, fastqs, barcodes, output_file, batches, max_batches, threads, logfile):
    """
    Search for barcodes in fastq reads, write summary to file.
    """
    # Global logger setup
    logger = setup_logger(logfile)

    # import seqspec.yaml
    spec = load_spec(yaml)

    # Run seqspec get_index_by_primer
    # This structures fastq_info in an appropriate format for different tools
    fastq_info = region_indices(spec, fastqs)

    # Load barcodes
    expected_barcodes = parse_seed_arguments(barcodes)  
    logger.info(f"Expected barcodes")
    for key, barcode in expected_barcodes.items():
        logger.info(f"{key}: {barcode}")

    # Open file for writing output
    if output_file:
        f = open(f"{output_file}", 'w')
        f.write("read\tname\tbarcode_whitelist\tbarcode\torientation\tstart\tend\tmismatches\n")

    # Extract elements from sequencing reads
    process_paired_fastq_batches(fastq_info = fastq_info, batch_size = batches, max_batches = max_batches,
                                 num_workers = threads, region_ids = None, target = None, output_handler = f,
                                 barcodes = expected_barcodes)

    return 


@log_errors
def parse_seed_arguments(barcode_args):
    """
    Parse seed arguments from command line.
    
    Args:
        barcode_args (List[str]): List of barcode arguments in format 'KEY:FILE'
    
    Returns:
        Dict[str, List[str]]: Dictionary of barcodes with keys as region identifiers
    """

    expected_barcodes = {}
    
    for arg in barcode_args:
        try:
            # Split the argument into key and file path
            key, file_path = arg.split(':')
            
            # Read barcodes from the file
            barcodes = read_barcode_file(file_path)
            
            # Store barcodes in the dictionary
            if barcodes:
                expected_barcodes[key] = barcodes
                logger.info(f"Loaded {len(barcodes)} barcodes for {key} from {file_path}")
            else:
                logger.warning(f"No barcodes loaded for {key} from {file_path}")
        
        except ValueError:
            logger.error(f"Invalid barcode argument format: {arg}. Use 'KEY:FILE'")
    
    return expected_barcodes


@log_errors
def read_barcode_file(file_path):
    """
    Read barcode sequences from a text file.
    
    Args:
        file_path (str): Path to the barcode file
    
    Returns:
        List[str]: List of unique barcode sequences
    """

    try:
        with open(file_path, 'r') as f:
            # Read lines, strip whitespace, remove empty lines
            barcodes = [line.strip() for line in f if line.strip()]
        
        # Remove duplicates while preserving order
        unique_barcodes = list(dict.fromkeys(barcodes))
        
        if not unique_barcodes:
            logger.warning(f"No barcodes found in file: {file_path}")
        
        return unique_barcodes
    
    except FileNotFoundError:
        logger.error(f"Barcode file not found: {file_path}")
        return []
    except Exception as e:
        logger.error(f"Error reading barcode file {file_path}: {e}")
        return []