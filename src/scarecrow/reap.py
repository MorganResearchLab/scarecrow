# -*- coding: utf-8 -*-
"""
#!/usr/bin/env python3
@author: David Wragg
"""

import ast
import gzip as gz
import json
import logging
import os
import gc
import pandas as pd
import pysam
import numpy as np
import shutil
import multiprocessing as mp
from argparse import RawTextHelpFormatter
from collections import defaultdict
from functools import lru_cache
from typing import List, Tuple, Optional, Dict, Set
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string


class BarcodeMatcherOptimized:
    def __init__(self, barcode_files: Dict[str, str], mismatches: int = 1, base_quality_threshold: int = None):
        """
        Initialize matcher with JSON whitelist files
        
        Args:
            barcode_files: Dict mapping whitelist names to JSON file paths
            mismatches: Maximum number of mismatches allowed (should match JSON data)
            base_quality_threshold: Quality threshold for filtering
        """
        self.mismatches = mismatches
        self.base_quality_threshold = base_quality_threshold
        self.matchers = {}
        
        # Load and organize barcode data from JSON files
        for whitelist, json_file in barcode_files.items():
            with open(json_file) as f:
                barcode_data = json.load(f)
            
            # Create efficient lookup structures
            exact_matches = set(barcode_data["0"].keys())
            mismatch_lookup = {}
            
            # Build mismatch lookup from level 1 data if mismatches are allowed
            if self.mismatches > 0 and "1" in barcode_data:
                for variant, original in barcode_data["1"].items():
                    mismatch_lookup[variant] = original
            
            # Store precomputed data structures
            barcode_len = len(next(iter(exact_matches)))
            self.matchers[whitelist] = {
                'exact': exact_matches,
                'mismatches': mismatch_lookup,
                'length': barcode_len
            }

    @lru_cache(maxsize=1024)
    def _reverse_complement(self, sequence: str) -> str:
        """Cached reverse complement computation"""
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        return ''.join(complement.get(base, base) for base in reversed(sequence))

    def _get_sequence_with_jitter(self, full_sequence: str, start: int, end: int, jitter: int) -> List[Tuple[str, int]]:
        """Generate possible sequences with position jitter"""
        barcode_length = end - start + 1
        sequences = []
        
        min_start = max(0, start - jitter - 1)
        max_start = min(len(full_sequence) - barcode_length + 1, start + jitter)        
        sequences = [
            (full_sequence[adj_start:adj_start + barcode_length], adj_start)
            for adj_start in range(min_start, max_start)
        ]
        return sequences

    def find_match(self, sequence: str, quality_scores: str, whitelist: str, orientation: str,
                  original_start: int, original_end: int, jitter: int, base_quality: int) -> Tuple[str, int, int]:
        """
        Find best matching barcode using precomputed exact and mismatch data
        
        Returns:
            Tuple of (matched_barcode, mismatch_count, position)
        """
        # Convert whitelist string from CSV to actual whitelist key
        if isinstance(whitelist, str):
            try:
                whitelist = ast.literal_eval(whitelist)
                if isinstance(whitelist, list) and len(whitelist) == 1 and isinstance(whitelist[0], tuple):
                    whitelist = whitelist[0][0]  # Extract the actual whitelist name
            except (ValueError, SyntaxError):
                pass

        # Check if whitelist exists in matchers
        if whitelist not in self.matchers:
            return 'null', -1, original_start - 1
            
        # Apply quality filtering if needed
        if base_quality is not None:
            sequence = filter_low_quality_bases(sequence, quality_scores, base_quality)[0]
        
        matcher = self.matchers[whitelist]
        
        # Handle reverse orientation
        if orientation == 'reverse':
            sequence = self._reverse_complement(sequence)
            
        # Get possible sequences with jitter
        possible_sequences = self._get_sequence_with_jitter(sequence, original_start, original_end, jitter)
        
        # Try exact matches first
        closest_match = None
        min_distance = float('inf')

        for seq, pos in possible_sequences:
            if seq in matcher['exact']:
                distance = abs(pos - (original_start - 1))
                if distance < min_distance:
                    min_distance = distance
                    closest_match = (seq, pos)

        if closest_match:
            return closest_match[0], 0, closest_match[1]
        
        # Try one-mismatch sequences if allowed
        if self.mismatches > 0:
            best_match = 'null'
            best_position = original_start - 1
            best_score = (1, float('inf'))  # (mismatches, position_distance)
            
            for seq, pos in possible_sequences:
                if seq in matcher['mismatches']:
                    pos_distance = abs(pos - (original_start - 1))
                    current_score = (1, pos_distance)
                    
                    if current_score < best_score:
                        best_score = current_score
                        best_match = matcher['mismatches'][seq]
                        best_position = pos
                        
            if best_match != 'null':
                return best_match, 1, best_position
        
        return 'null', self.mismatches + 1, original_start - 1
    

def parser_reap(parser):
    subparser = parser.add_parser(
        "reap",
        description="""
Extract sequence range from one fastq file of a pair, and annotate the sequence header with barcode 
sequences based on predicted positions generated by scarecrow harvest.

Example:

scarecrow reap --fastqs R1.fastq.gz R2.fastq.gz\n\t--barcode_positions barcode_positions.csv\n\t--barcodes\tBC1:v1_whitelist:bc1_whitelist.txt\n\t\t\tBC2:v2_whitelist:bc2_whitelist.txt\n\t\t\tBC3:v1_whitelist:bc3_whitelist.txt\n\t--read1 0-64 \n\t--out cdna.fastq.gz
---
""",
        epilog="The --barcodes <name> must match the barcode_whitelist values in the --barcode_positions file.",
        help="Extract sequence range from fastq files",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument(
        "--fastqs", 
        nargs="+", 
        help="Pair of FASTQ files")
    subparser.add_argument(
        "-o", "--out",
        metavar="<file>",
        help=("Path to output fastq file"),
        type=str,
        default='extracted.fastq',
    )
    subparser.add_argument(
        "-p", "--barcode_positions",
        metavar="<file>",
        help=("File containing barcode positions, output by scarecrow harvest"),
        type=str,
        default=[],
    )
    subparser.add_argument(
        "-j", "--jitter",
        metavar="<int>",
        type=int,
        default=5,
        help='Barcode position jitter [5]',
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
        "-x", "--extract",
        metavar="<range>",
        help=("Sequence range to extract <read>:<range> (e.g. 1:1-64)"),
        type=str,
        default=None
    )
    subparser.add_argument(
        "-u", "--umi",
        metavar="<range>",
        help=("Sequence range to extract for UMI <read>:<range> (e.g. 2:1-10)"),
        type=str,
        default=None
    )
    subparser.add_argument(
        "-c", "--barcodes",
        metavar="<string>",
        nargs='+', 
        help='Barcode whitelist files in format <barcode_name>:<whitelist_name>:<whitelist_file>\n\t(e.g. BC1:v1:barcodes1.txt BC2:v2:barcodes2.txt ...)',
    )
    subparser.add_argument(
        "-b", "--batch_size",
        metavar="<int>",
        help=("Number of read pairs per batch to process at a time [10000]"),
        type=int,
        default=10000,
    )
    subparser.add_argument(
        "-@", "--threads",
        metavar="<int>",
        help=("Number of processing threads [1]"),
        type=int,
        default=1,
    )
    subparser.add_argument(
        "-r", "--barcode_reverse_order",
        action='store_true',
        help='Reverse retrieval order of barcodes in barcode positions file [false]'
    )
    subparser.add_argument(
        "-z", "--gzip",
        action='store_true',
        help='Compress (gzip) fastq output [false]'
    )
    subparser.add_argument(
        "-v", "--verbose",
        action='store_true',
        help='Enable verbose output [false]'
    )
    return subparser

def validate_reap_args(parser, args) -> None:
    """ 
    Validate arguments 
    """
    run_reap(fastqs = [f for f in args.fastqs], 
             barcode_positions = args.barcode_positions,
             barcode_reverse_order = args.barcode_reverse_order,
             output = args.out,
             extract = args.extract,
             umi = args.umi, 
             barcodes = args.barcodes,
             jitter = args.jitter,
             mismatches = args.mismatches,
             base_quality = args.base_quality,
             batches = args.batch_size, 
             threads = args.threads,
             verbose = args.verbose,
             gzip = args.gzip)

@log_errors
def run_reap(fastqs: List[str], 
             barcode_positions: str = None,
             barcode_reverse_order: bool = False,
             output: str = 'extracted.fastq',
             extract: str = None,
             umi: Optional[str] = None,
             barcodes: List[str] = None,
             jitter: int = 5,
             mismatches: int = 1,
             base_quality: int = None,
             batches: int = 10000,
             threads: int = 4,
             verbose: bool = False,
             gzip: bool = False) -> None:
    """
    Main function to extract sequences with barcode headers
    """    
    # Global logger setup
    logfile = '{}_{}.{}'.format('./scarecrow_reap', generate_random_string(), 'log')
    logger = setup_logger(logfile)
    logger.info(f"logfile: '{logfile}'")
    logger.info(f"fastqs: {fastqs}")
    logger.info(f"output: '{output}'")
    logger.info(f"jitter: '{jitter}'")
    logger.info(f"mismatches: '{mismatches}'")
    logger.info(f"base_quality: '{base_quality}'")

    # Extract barcodes and convert whitelist to set
    barcode_files = parse_seed_arguments(barcodes)
    if verbose:
        for whitelist, filename in barcode_files.items():
            logger.info(f"{whitelist}: {filename}")

    # Check if the output filename string ends with .gz
    if output.endswith(".gz"):
        output = output[:-3]

    # Extract sequences
    extract_sequences(
        fastq_files = [f for f in fastqs],
        barcode_positions_file = barcode_positions,
        barcode_reverse_order = barcode_reverse_order,
        barcode_files = barcode_files,
        output = output,
        extract = extract,
        umi = umi,
        jitter = jitter,
        mismatches = mismatches,
        base_quality = base_quality,
        batch_size = batches,
        threads = threads,
        verbose = verbose
    )
    
    # gzip
    if gzip:
        logger.info(f"Compressing '{output}'")
        with open(output, 'rb') as f_in, gz.open(output + ".gz", 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(output)

@log_errors
def process_read_batch(read_batch: List[Tuple], 
                      barcode_configs: List[Dict],
                      matcher: BarcodeMatcherOptimized,
                      read_range: Tuple[int, int],
                      read_index: int,
                      umi_index: int,
                      umi_range: Tuple[int, int],
                      base_quality: int = None,
                      jitter: int = 5,
                      verbose: bool = False) -> List[str]:
    """
    Modified process_read_batch to handle jitter and improved matching
    """
    logger = logging.getLogger('scarecrow')
    output_entries = []
    read_count = len(read_batch) # Count reads in this batch
    
    for reads in read_batch:
        barcodes = []
        positions = []
        mismatches = []
        for config in barcode_configs:
            seq = reads[config['file_index']].sequence
            start, end = config['start'], config['end']            

            whitelist = ast.literal_eval(config['whitelist'])[1]
            if isinstance(whitelist, list) and len(whitelist) == 1 and isinstance(whitelist[0], tuple):
                whitelist = whitelist[0]

            if verbose:
                print(f"Read: {reads[config['file_index']].name} {reads[config['file_index']].comment}")
                print(f"Sequence: {seq}")
                print(f"Looking for barcode in range {start}-{end} with jitter {jitter}")

            if whitelist in matcher.matchers.keys():
                matched_barcode, mismatch_count, adj_position = matcher.find_match(
                    seq, reads[config['file_index']].quality, whitelist, config['orientation'], 
                    start, end, jitter, base_quality)

                if verbose:
                    print(f"Matched barcode: {matched_barcode} with {mismatch_count} mismatches at position {adj_position}")
                
                barcodes.append(matched_barcode)
                positions.append(str(adj_position + 1))  # Convert to 1-based position for output
                mismatches.append(str(mismatch_count))
            else:
                barcodes.append('null')
                positions.append(str(start))
                mismatches.append('NA')

        # Create output entry with original and adjusted positions
        source_entry = reads[read_index]
        # Apply base quality filtering to extracted sequence
        if base_quality is not None:
            filtered_seq, filtered_qual = filter_low_quality_bases(
                source_entry.sequence[read_range[0]:read_range[1]], 
                source_entry.quality[read_range[0]:read_range[1]], base_quality)
        else:
            filtered_seq = source_entry.sequence[read_range[0]:read_range[1]]
            filtered_qual = source_entry.quality[read_range[0]:read_range[1]]
        
        header = f"@{source_entry.name} {source_entry.comment}"
        header += f" barcodes={('_').join(barcodes)}"
        header += f" positions={('_').join(positions)}"
        header += f" mismatches={('_').join(mismatches)}"

        if umi_index is not None:
            umi_seq = reads[umi_index].sequence[umi_range[0]:umi_range[1]]
            header += f" UMI={umi_seq}"

        #output_entries.append(f"{header}\n{extract_seq}\n+\n{extract_qual}\n")
        output_entries.append(f"{header}\n{filtered_seq}\n+\n{filtered_qual}\n")


    return output_entries, read_count

@log_errors
def extract_sequences(
    fastq_files: List[str] = None,
    barcode_positions_file: str = None,
    barcode_reverse_order: bool = False,
    barcode_files: Dict[str, str] = None,
    output: str = 'extracted.fastq.gz',
    extract: str = None,
    umi: Optional[str] = None,
    jitter: int = 5,
    mismatches: int = 1,
    base_quality: int = None,
    batch_size: int = 100000,
    threads: Optional[int] = None,
    verbose: bool = False
) -> None:
    """
    Modified to use JSON whitelist files with mismatch parameter
    """
    logger = logging.getLogger('scarecrow')

    # Initialize configurations
    barcode_positions = pd.read_csv(barcode_positions_file)
    if barcode_reverse_order:
        barcode_positions = barcode_positions[::-1].reset_index(drop=True)
    barcode_configs = prepare_barcode_configs(barcode_positions, jitter)

    # Extract range
    extract_index, extract_range = extract.split(':')
    extract_index = int(extract_index)-1
    extract_range = parse_range(extract_range)
    logger.info(f"FASTQ sequence range to extract: '{extract}'")

    # UMI range
    if umi is not None:
        umi_index, umi_range = umi.split(':')
        umi_index = int(umi_index)-1
        umi_range = parse_range(umi_range)
        logger.info(f"UMI sequence range to extract: '{umi}'")
    else:
        umi_index = None
        umi_range = None

    # Create matcher with JSON whitelist files
    logger.info(f"Generating barcode matcher from JSON whitelists")
    matcher = BarcodeMatcherOptimized(
        barcode_files = barcode_files,
        mismatches = mismatches,
        base_quality_threshold = base_quality
    )
       
    # Process files with minimal overhead
    if threads is None:
        threads = min(mp.cpu_count() - 1, 8)
    else:
        threads = min(threads, mp.cpu_count())
    logger.info(f"Using {threads} threads")

    # List of files generated
    files = []
    total_reads = 0 # Initialize total read counter

    def write_and_clear_results(jobs):
        """Retrieve results from completed jobs, write to file, and free memory."""
        nonlocal total_reads  # Access the outer scope counter
        for idx, job in enumerate(jobs):
            results, batch_count = job.get()
            total_reads += batch_count  # Update total count
            outfile = output + "_" + str(idx)
            if outfile not in files:
                files.append(outfile)
                if os.path.exists(outfile):
                    os.remove(outfile)
            with open(outfile, 'a') as out_fastq:
                out_fastq.writelines(results)
            del results                
        jobs.clear()
        gc.collect()
    
    def combine_results_chunked(files, output, chunk_size=1024*1024):
        with open(output, 'w') as out_fastq:
            for fastq_file in files:
                with open(fastq_file, 'r') as in_fastq:
                    while chunk := in_fastq.read(chunk_size):
                        out_fastq.write(chunk)
                os.remove(fastq_file)

    logger.info(f"Processing reads")
    with pysam.FastqFile(fastq_files[0]) as r1, \
         pysam.FastqFile(fastq_files[1]) as r2:
        
        # Create batches efficiently
        read_pairs = zip(r1, r2)        

        pool = mp.Pool(threads)
        batch = []
        jobs = []
        counter = 0       

        for reads in read_pairs:            
            batch.append(reads)            
            if len(batch) >= batch_size:
                jobs.append(pool.apply_async(worker_task, args=((batch, barcode_configs, matcher, 
                                                                 extract_range, extract_index, 
                                                                 umi_index, umi_range, 
                                                                 base_quality, jitter, verbose),)))
                counter += len(batch)                

                # Write results and free memory if jobs exceed a threshold
                if len(jobs) >= threads:
                    write_and_clear_results(jobs)
                    logger.info(f"Processed {counter} reads")
                
                # Clear the current batch
                batch = []                
        
        # Process remaining reads
        if batch:
            jobs.append(pool.apply_async(worker_task, args=((batch, barcode_configs, matcher, 
                                                                 extract_range, extract_index, 
                                                                 umi_index, umi_range, 
                                                                 base_quality, jitter, verbose),)))
            
        write_and_clear_results(jobs)
        logger.info(f"Processed {counter} reads")
        
        # Close pools
        pool.close()
        pool.join()

    # Combine results
    logger.info(f"Combining results: {files}")
    combine_results_chunked(files, output)

    # Log final count
    logger.info(f"Total reads processed: {total_reads}")



def parse_range(range_str: str) -> Tuple[int, int]:
    """
    Parse range string
    """
    start, end = map(int, range_str.split('-'))
    start = max(0, start -1)
    return (start, end)

def parse_seed_arguments(barcodes: List[str]) -> Dict[str, str]:
    """
    Modified to handle JSON whitelist files instead of text files
    Format: <barcode_name>:<whitelist_name>:<whitelist_file.json>
    Returns: Dict mapping whitelist names to JSON file paths
    """
    barcode_files = {}
    
    if barcodes:
        for barcode in barcodes:
            name, whitelist, filename = barcode.split(':')
            if not os.path.exists(filename):
                raise FileNotFoundError(f"Barcode whitelist file not found: {filename}")
            barcode_files[whitelist] = filename
            
    return barcode_files

@lru_cache(maxsize=1024)
def filter_low_quality_bases(sequence: str, quality: str, threshold: int) -> Tuple[str, str]:
    """
    Efficiently filter low-quality bases with minimal overhead
    """
    # Pre-compute quality scores to avoid repeated calculation
    qual_scores = [ord(q) - 33 for q in quality]    
    # Use list comprehension for efficient filtering
    filtered_seq = ''.join(
        base if score >= threshold else 'N'
        for base, score in zip(sequence, qual_scores)
    )    
    return filtered_seq, quality

def prepare_barcode_configs(positions: pd.DataFrame, jitter: int) -> List[Dict]:
    """
    Prepare barcode configurations
    """
    return [{
        'index': idx,
        'file_index': 0 if row['read'] == 'read1' else 1,
        'start': row['start'],
        'end': row['end'],
        'orientation': row['orientation'],
        'whitelist': row['barcode_whitelist']
    } for idx, row in positions.iterrows()]

def worker_task(args):
    entries, count = process_read_batch(*args)
    return entries, count
