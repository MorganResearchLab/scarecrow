#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
"""

import os
import pysam
from argparse import RawTextHelpFormatter
from concurrent.futures import ProcessPoolExecutor, as_completed
from seqspec.utils import load_spec
from scarecrow.fastq_logging import log_errors, setup_logger, logger
from scarecrow.tools import FastqProcessingError, reverse_complement, count_fastq_reads, region_indices, generate_random_string
from typing import List, Dict, Optional, Set, Union
from tqdm import tqdm

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
        help='Barcode whitelist files in format <name>:<whitelist>:<file> (e.g. BC1:v1:barcodes1.txt BC2:n198:barcodes2.txt)',
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
    return subparser

def validate_seed_args(parser, args):
    run_seed(yaml = args.yaml, 
             fastqs = [f for f in args.fastqs],
             barcodes = args.barcodes, 
             output_file = args.out, 
             batches = args.batch_size, 
             threads = args.threads, 
             max_batches = args.max_batches)

@log_errors
def run_seed(yaml, fastqs, barcodes, output_file, batches, max_batches, threads):
    """
    Search for barcodes in fastq reads, write summary to file.
    """
    # Global logger setup
    logfile = '{}_{}.{}'.format('./scarecrow', generate_random_string(), 'log')
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
    process_read_pair_batches(fastq_info = fastq_info, batch_size = batches, max_batches = max_batches,
                              num_workers = threads, output_handler = f, barcodes = expected_barcodes)

    return 


@log_errors
def parse_seed_arguments(barcode_args):
    """
    Parse seed arguments from command line.
    
    Args:
        barcode_args (List[str]): List of barcode arguments in format 'KEY:WHITELIST:FILE'
    
    Returns:
        Dict[str, List[str]]: Dictionary of barcodes with keys as region identifiers
    """

    expected_barcodes = {}
    
    for arg in barcode_args:
        try:
            # Split the argument into key, whitelist and file path
            key, label, file_path = arg.split(':')
            
            # Read barcodes from the file
            barcodes = read_barcode_file(file_path)
            
            # Store barcodes in the dictionary
            if barcodes:
                expected_barcodes[key,label] = barcodes                
                logger.info(f"Loaded {len(barcodes)} barcodes for {key} from {label} at {file_path}")
            else:
                logger.warning(f"No barcodes loaded for {key} from {label} at {file_path}")
        
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


@log_errors
def process_read_pair_batches(
    fastq_info: List[Dict],
    batch_size: int = 10000,
    max_batches: Optional[int] = None,
    output_handler: Optional[str] = None,
    num_workers: int = None,
    barcodes: Dict = None
) -> Dict[str, int]:
    """
    Process paired-end FASTQ files with improved parallel and memory-efficient processing.
    
    Args:
        fastq_info (List[Dict]): FASTQ file information
        batch_size (int): Number of read pairs per batch
        max_batches (int, optional): Limit on number of batches
        output_handler (str, optional): Output handling method
        region_ids (List[str], optional): Regions to extract
        num_workers (int, optional): Number of parallel processing workers
        barcodes (Dict, optional): Barcode information
        target (str, optional): Target specification
    
    Returns:
        Dict[str, int]: Barcode count distribution
    """
    # Use all available CPU cores if not specified
    if num_workers is None:
        num_workers = multiprocessing.cpu_count()
    
    # Prepare chunk arguments
    r1_file_path = list(fastq_info[0].keys())[0]
    r1_regions = fastq_info[0][r1_file_path]
    r1_strand = fastq_info[0]['strand']
    
    r2_file_path = list(fastq_info[1].keys())[0]
    r2_regions = fastq_info[1][r2_file_path]
    r2_strand = fastq_info[1]['strand']
    
    # Chunk arguments for parallel processing
    chunk_args = (
        r1_file_path, r2_file_path, 
        batch_size, max_batches, barcodes,
        r1_regions, r2_regions,
        r1_strand, r2_strand
    )
    
    barcode_counts = {}
    total_read_pairs = 0
    
    try:
        # Parallel processing of FASTQ chunks
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            # Submit chunk processing
            future = executor.submit(process_fastq_chunk, chunk_args)
            
            # Process results
            batches = future.result()
            
            # Progress tracking
            with tqdm(total=len(batches), desc="Processing batches") as pbar:
                for batch in batches:
                    for read_pair in batch:
                        try:
                            write_barcodes_CSV(read_pair, output_handler)
                        
                        except Exception as e:
                            logger.error(f"Error processing read pair: {e}")
                    
                    total_read_pairs += len(batch)
                    pbar.update(1)
        
        return barcode_counts
    
    except Exception as e:
        logger.error(f"Processing failed: {e}")
        raise

@log_errors
def process_fastq_chunk(chunk_args):
    """
    Process a chunk of FASTQ files with improved parallel processing.
    
    Args:
        chunk_args (tuple): Contains:
            - r1_file_path (str): Path to Read 1 FASTQ file
            - r2_file_path (str): Path to Read 2 FASTQ file
            - batch_size (int): Number of read pairs per batch
            - max_batches (Optional[int]): Limit on number of batches
            - barcodes (Optional[Dict]): Barcode information
            - r1_regions (List): Regions for Read 1
            - r2_regions (List): Regions for Read 2
            - r1_strand (str): Strand information for Read 1
            - r2_strand (str): Strand information for Read 2
    
    Returns:
        List[Dict]: Processed batches of read pair information
    """
    r1_file_path, r2_file_path, batch_size, max_batches, barcodes, \
    r1_strand, r2_strand = chunk_args
    
    processed_batches = []
    
    try:
        reads = count_fastq_reads(r1_file_path)
        with tqdm(total = reads, desc = "Processing reads") as pbar:

            with pysam.FastxFile(r1_file_path) as r1_fastq, \
                pysam.FastxFile(r2_file_path) as r2_fastq:

                batch = []
                batch_count = 0
                
                while True:
                    try:
                        r1_entry = next(r1_fastq)
                        r2_entry = next(r2_fastq)
                    except StopIteration:
                        break
                    
                    # Barcode or region extraction
                    r1_region_details = extract_barcodes(r1_entry, rev = False, expected_barcodes = barcodes)                    
                    r2_region_details = extract_barcodes(r2_entry, rev = True, expected_barcodes = barcodes)
                    
                    # Combine pair information
                    read_pair_info = {
                        'read1': {
                            'file': os.path.basename(r1_file_path),
                            'strand': r1_strand,
                            'header': r1_entry.name,
                            'regions': r1_region_details
                        },
                        'read2': {
                            'file': os.path.basename(r2_file_path),
                            'strand': r2_strand,
                            'header': r2_entry.name,
                            'regions': r2_region_details
                        }
                    }
                    
                    batch.append(read_pair_info)
                    pbar.update(1)
                    
                    # Yield batch when full
                    if len(batch) >= batch_size:
                        processed_batches.append(batch)
                        batch = []
                        batch_count += 1
                        
                        # Optional batch limit
                        if max_batches and batch_count >= max_batches:
                            break
                
                # Yield final batch if not empty
                if batch:
                    processed_batches.append(batch)
        
        return processed_batches
    
    except IOError as e:
        logger.error(f"File reading error: {e}")
        raise FastqProcessingError(f"Unable to process FASTQ files: {e}")


@log_errors
def extract_barcodes(
    entry: pysam.FastxRecord,
    rev: bool = False,
    expected_barcodes: Optional[Union[Dict[str, List[str]], List[str]]] = None
) -> List[Dict]:
    """
    Extract either specific regions or find barcodes in a sequence.
    
    Args:
        entry (pysam.FastxRecord): Input sequencing read
        regions (Optional[List[RegionCoordinate]]): Regions to extract coordinates
        expected_barcodes (Optional[Union[Dict[str, List[str]], List[str]]]): 
            Barcodes to cross-reference
        rev (bool): Whether to reverse the sequence
    
    Returns:
        List[Dict]: Extracted region or barcode details
    
    Raises:
        ValueError: If neither regions nor barcodes are provided
    """
    try:
        # Validate input
        if expected_barcodes is None:
            raise ValueError("Must provide whitelist of barcodes to search")
        
        # Prepare sequence
        full_sequence = entry.sequence
                
        # Find barcodes
        # Prepare variables to track dictionary keys
        barcode_dict_keys = {}
        if isinstance(expected_barcodes, dict):
            # Create a reverse mapping of barcodes to their original dictionary keys
            barcode_dict_keys = {}
            for key, barcode_list in expected_barcodes.items():
                for barcode in barcode_list:
                    if barcode not in barcode_dict_keys:
                        barcode_dict_keys[barcode] = []
                    barcode_dict_keys[barcode].append(key)
            
            # Flatten barcodes for searching
            barcodes_to_search = list(barcode_dict_keys.keys())
        elif isinstance(expected_barcodes, list):
            barcodes_to_search = expected_barcodes
        else:
            raise ValueError(f"Unexpected barcode type: {type(expected_barcodes)}")
        
        # Find barcode positions
        barcode_matches = find_barcode_positions(full_sequence, barcodes_to_search)

        logger.info(f">{entry.name} {entry.comment}")
        logger.info(f"{full_sequence}")
        
        # Prepare barcode match details
        barcode_details = []
        for match in barcode_matches:
            barcode_detail = {
                'barcode': match['barcode'],
                'orientation': match['orientation'],
                'sequence': match['sequence'],
                'start': match['start'],
                'end': match['end'],
                'mismatches': match['mismatches']
            }
            
            # Add all matching dictionary keys if available
            if barcode_dict_keys and match['barcode'] in barcode_dict_keys:
                barcode_detail['dict_keys'] = barcode_dict_keys[match['barcode']]
            elif barcode_dict_keys and reverse_complement(match['barcode']) in barcode_dict_keys:
                barcode_detail['dict_keys'] = barcode_dict_keys[reverse_complement(match['barcode'])]
            else:
                barcode_detail['dict_keys'] = []

            
            barcode_details.append(barcode_detail)
            
            logger.info(f"...{barcode_detail['dict_keys']} ({match['barcode']}) hit: {match['sequence']}, "
                    f"Orientation: {match['orientation']}, "
                    f"Start: {match['start']}, "
                    f"End: {match['end']}, "
                    f"Mismatches: {match['mismatches']}")

        
        return barcode_details
            
    except Exception as e:
        logger.error(f"Error in extract_region_details: {e}")
        raise

@log_errors
def find_barcode_positions(sequence, barcodes, max_mismatches=1):
    """
    Find all positions of barcodes in a sequence with tolerance for mismatches.
    
    Args:
        sequence (str): The input DNA sequence to search
        barcodes (List[str]): List of expected barcode sequences
        max_mismatches (int): Maximum allowed mismatches when matching barcodes
    
    Returns:
        List[Dict]: A list of dictionaries with details of all barcode matches
    """
    def hamming_distance(s1, s2):
        """Calculate Hamming distance between two strings."""
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))
    
    # List to store all barcode matches
    barcode_matches = []
    orientations = ['forward', 'reverse']
    
    # Iterate through all possible start positions
    for start in range(len(sequence)):
        for barcode in barcodes:
            for orientation in orientations:
                # reverse complement barcode if testing reverse strand
                if orientation == 'reverse':
                    barcode = reverse_complement(barcode)
                # Check all possible end positions for this barcode
                for end in range(start + len(barcode), len(sequence) + 1):
                    candidate = sequence[start:end]
                    
                    # Check if candidate matches barcode within max mismatches
                    if len(candidate) == len(barcode):
                        mismatches = hamming_distance(candidate, barcode)
                        if mismatches <= max_mismatches:
                            match_details = {
                                'barcode': barcode,
                                'orientation': orientation,
                                'sequence': candidate,
                                'start': start,
                                'end': end,
                                'mismatches': mismatches
                            }
                            barcode_matches.append(match_details)
    
    # Sort matches by start position
    barcode_matches.sort(key=lambda x: x['start'])
    
    return barcode_matches

@log_errors
def write_barcodes_CSV(read_pair: List, output_handler):
    """
    Write barcodes found to CSV file.
    
    Args:
        read_key (List): List of barcodes found in read pairs
        output_file (str): Output file path
    """
    for read_key in ['read1', 'read2']:
    # Get barcodes for current read
        read_barcodes = read_pair.get(read_key, {}).get('regions', [])
        read = read_pair.get(read_key, {}).get('header', str)
        for sub_barcode in read_barcodes:
            bc = sub_barcode['dict_keys']
            barcode = sub_barcode['barcode']
            orientation = sub_barcode['orientation']
            start = sub_barcode['start']
            end = sub_barcode['end']
            mm = sub_barcode['mismatches']
            output_handler.write(f"{read_key}\t{read}\t{bc}\t{barcode}\t{orientation}\t{start}\t{end}\t{mm}\n")