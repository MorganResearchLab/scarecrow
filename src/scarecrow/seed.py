#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
"""

import os
import pysam
import multiprocessing as mp
from argparse import RawTextHelpFormatter
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import FastqProcessingError, reverse_complement, count_fastq_reads, region_indices, generate_random_string
from typing import List, Dict, Optional, Set, Union
from rich.progress import track

def parser_seed(parser):
    subparser = parser.add_parser(
        "seed",
        description="""
Search fastq reads for barcodes in whitelists

Example:

scarecrow seed --fastqs R1.fastq.gz R2.fastq.gz\n\t--strands pos neg\n\t--barcodes BC1:BC1.txt BC2:BC2.txt BC3:BC3.txt\n\t--out barcode_counts.csv 
---
""",
        help="Search fastq reads for barcodes",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument(
        "--fastqs", 
        nargs="+", 
        help="Pair of FASTQ files")
    subparser.add_argument(
        "--strands", 
        metavar="strands",
        nargs="+", 
        default=['pos', 'neg'],
        help="Orientations of FASTQ files (e.g. pos neg)")
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
        "-@", "--threads",
        metavar="threads",
        help=("Number of processing threads [4]"),
        type=int,
        default=4,
    )
    return subparser

def validate_seed_args(parser, args):
    run_seed(fastqs = [f for f in args.fastqs],
             strands = [s for s in args.strands],
             barcodes = args.barcodes, 
             output_file = args.out, 
             batches = args.batch_size, 
             threads = args.threads)

@log_errors
def run_seed(fastqs, strands, barcodes, output_file, batches, threads):
    """
    Search for barcodes in fastq reads, write summary to file.
    """
    # Global logger setup
    logfile = '{}_{}.{}'.format('./scarecrow_seed', generate_random_string(), 'log')
    logger = setup_logger(logfile)
    logger.info(f"logfile: ${logfile}")

    # Originally we imported a seqspec yaml and called region_indices (scarecrow tools)
    # which would return the following structure:
    #[
    #{
    #    "R1.fastq": [  # fastq filename as key
    #        # List of RegionCoordinate objects for each region
    #        # From the YAML, these would include regions like I5, R1, cdna, BC1, etc.
    #        RegionCoordinate(
    #            region_id: str,  # e.g. "BC1", "cdna", etc.
    #            start: int,      # 0-indexed start position
    #            stop: int        # 0-indexed stop position
    #        ),
    #        # ... more regions
    #    ],
    #    "strand": "pos"  # or "neg" depending on the read
    #},
    ## ... one dict per fastq file
    #]

    # This was achieved as follows:
    
    # spec = load_spec(yaml)
    # fastq_info = region_indices(spec, fastqs)

    # However, we only use the fastq and strand information, so don't need
    # the level of complexity of the YAML file. However, for possible
    # future-compataibility, have retained the same structure when reading
    # in the fastqs and strands information:
    fastq_info = create_fastq_info(fastqs, strands)
    logger.info(f"{fastq_info}")

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
    process_read_pair_batches(fastq_info = fastq_info, batch_size = batches, output_handler = f,
                              num_workers = threads, barcodes = expected_barcodes)

    return 

def create_fastq_info(fastqs, strands):
    return [
        {fastq: [], "strand": strand}
        for fastq, strand in zip(fastqs, strands)
    ]

@log_errors
def parse_seed_arguments(barcode_args):
    """
    Parse seed arguments from command line.
    
    Args:
        barcode_args (List[str]): List of barcode arguments in format 'KEY:WHITELIST:FILE'
    
    Returns:
        Dict[str, List[str]]: Dictionary of barcodes with keys as region identifiers
    """
    logger = logging.getLogger('scarecrow')

    expected_barcodes = {}
    
    for arg in barcode_args:
        try:
            # Split the argument into key, whitelist and file path
            key, label, file_path = arg.split(':')
            
            # Read barcodes from the file
            if os.path.exists(file_path):
                barcodes = read_barcode_file(file_path)
            else:
                logger.warning(f"File not found: {file_path}")
                raise RuntimeError(f"File not found: {file_path}")
                        
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
    logger = logging.getLogger('scarecrow')

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
    output_handler: Optional[str] = None,
    num_workers: int = None,
    barcodes: Dict = None
) -> Dict[str, int]:
    """
    Process paired-end FASTQ files with improved parallel and memory-efficient processing.
    
    Args:
        fastq_info (List[Dict]): FASTQ file information
        batch_size (int): Number of read pairs per batch
        output_handler (str, optional): Output handling method
        num_workers (int, optional): Number of parallel processing workers
        barcodes (Dict, optional): Barcode information
    
    Returns:
        Dict[str, int]: Barcode count distribution
    """
    logger = logging.getLogger('scarecrow')

    # Use all available CPU cores if not specified
    if num_workers is None:
        num_workers = mp.cpu_count()
    
    # Prepare chunk arguments for parallel processing
    chunk_args = (list(fastq_info[0].keys())[0],
                  list(fastq_info[1].keys())[0],
                  batch_size, barcodes,
                  fastq_info[0]['strand'],
                  fastq_info[1]['strand']
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
            for batch in track(batches, total = len(batches)):
                for read_pair in batch:
                    try:
                        write_barcodes_CSV(read_pair, output_handler)
                    
                    except Exception as e:
                        logger.error(f"Error processing read pair: {e}")
                
                total_read_pairs += len(batch)

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
            - barcodes (Optional[Dict]): Barcode information
            - r1_strand (str): Strand information for Read 1
            - r2_strand (str): Strand information for Read 2
    
    Returns:
        List[Dict]: Processed batches of read pair information
    """
    logfile = '{}_{}.{}'.format('./scarecrow_seed_barcodes', generate_random_string(), 'log')
    logger = setup_logger(logfile)
    logger.info(f"logfile: ${logfile}")

    # Chunk arguments
    r1_file_path, r2_file_path, batch_size, barcodes, r1_strand, r2_strand = chunk_args
    
    # Test if fastq files exist
    if os.path.exists(r1_file_path) is False:
        logger.warning(f"File not found: {r1_file_path}")        
        raise RuntimeError(f"File not found: {r1_file_path}")
    if os.path.exists(r2_file_path) is False:
        logger.warning(f"File not found: {r2_file_path}")
        raise RuntimeError(f"File not found: {r2_file_path}")
    
    processed_batches = []   
    try:
        reads = count_fastq_reads(r1_file_path)
        logger.info(f"Number of read pairs: {reads}")
        with pysam.FastxFile(r1_file_path) as r1_fastq, pysam.FastxFile(r2_file_path) as r2_fastq:

            batch = []
            batch_count = 0

            while True:

                for r1_entry, r2_entry in track(zip(r1_fastq, r2_fastq), description="Processing FastQ files", total=reads):
                                   
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
                    
                    # Yield batch when full
                    if len(batch) >= batch_size:
                        processed_batches.append(batch)
                        batch = []
                        batch_count += 1
                        
                # Yield final batch if not empty
                if batch:
                    processed_batches.append(batch)

                # Move onto next read
                try:
                    r1_entry = next(r1_fastq)
                    r2_entry = next(r2_fastq)
                except StopIteration:
                    break
        
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
    logger = logging.getLogger('scarecrow')

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
    logger = logging.getLogger('scarecrow')

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
                # Deduct 1 from barcode length when adding to start position so that
                # the extracted sequence is correct length
                for end in range(start + len(barcode)-1, len(sequence)+1):
                    candidate = sequence[start:end]
                    # Check if candidate matches barcode within max mismatches
                    if len(candidate) == len(barcode):
                        #logger.info(f"query range:{start}:{end}\tcandidate:{candidate}\tbarcode:{barcode}")
                        mismatches = hamming_distance(candidate, barcode)
                        if mismatches <= max_mismatches:
                            match_details = {
                                'barcode': barcode,
                                'orientation': orientation,
                                'sequence': candidate,
                                'start': start + 1, # for 1-based start
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
    logger = logging.getLogger('scarecrow')

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