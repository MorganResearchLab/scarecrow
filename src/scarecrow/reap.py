# -*- coding: utf-8 -*-
"""
#!/usr/bin/env python3
@author: David Wragg
"""

import ast
import gzip as gz
import itertools
import logging
import os
import gc
import pandas as pd
import pysam
import shutil
import multiprocessing as mp
from argparse import RawTextHelpFormatter
from collections import Counter
import csv
from functools import lru_cache
from typing import List, Tuple, Optional, Dict
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string


class BarcodeMatcherOptimized:
    def __init__(self, barcode_files: Dict[str, str], mismatches: int = 1, base_quality_threshold: int = None, verbose: bool = False):
        """
        Initialize matcher with text files containing barcodes (one per line)
        """
        self.mismatches = mismatches
        self.base_quality_threshold = base_quality_threshold
        self.verbose = verbose
        self.matchers = {}
        
        # Get logger but ensure it's properly configured first
        self.logger = setup_worker_logger()
        self.logger.info(f"Initializing matcher with whitelists: {list(barcode_files.keys())}")
        
        # Load and process barcode data from text files
        for whitelist, txt_file in barcode_files.items():
            self.logger.info(f"Loading whitelist '{whitelist}' from {txt_file}")
            
            # Read barcodes and generate mismatch dictionary
            with open(txt_file, 'r') as f:
                barcodes = [line.strip() for line in f if line.strip()]
            
            # Create exact matches set
            exact_matches = set(barcodes)
            
            # Generate mismatch lookup
            mismatch_lookup = {}
            if self.mismatches > 0:
                for barcode in barcodes:
                    mismatches = generate_mismatches(barcode, self.mismatches)
                    # We only use 1-mismatch sequences for now
                    if 1 in mismatches:
                        mismatch_lookup.update(mismatches[1])
            
            barcode_len = len(next(iter(exact_matches)))
            self.matchers[whitelist] = {
                'exact': exact_matches,
                'mismatches': mismatch_lookup,
                'length': barcode_len
            }
            
            if self.verbose:
                self.logger.info(f"Loaded {len(exact_matches)} exact matches and {len(mismatch_lookup)} variants for whitelist '{whitelist}'")
                self.logger.info(f"{whitelist} exact matches:\n{exact_matches}\n")
                self.logger.info(f"{whitelist} mismatches:\n{mismatch_lookup}\n")

    def _parse_whitelist_identifier(self, whitelist: str) -> str:
        """
        Parse whitelist identifier from various formats
        """
        original = whitelist
        
        try:
            if isinstance(whitelist, str):
                parsed = ast.literal_eval(whitelist)
                
                if isinstance(parsed, list) and len(parsed) == 1:
                    if isinstance(parsed[0], tuple):
                        whitelist = parsed[0][0]
                    else:
                        whitelist = parsed[0]
                elif isinstance(parsed, tuple):
                    whitelist = parsed[0]
                else:
                    whitelist = parsed
        except (ValueError, SyntaxError) as e:
            self.logger.debug(f"Could not parse whitelist identifier '{original}': {e}")
            pass
        
        if self.verbose:
            self.logger.info(f"Parsed whitelist identifier '{original}' -> '{whitelist}'")
        return whitelist

    def find_match(self, sequence: str, quality_scores: str, whitelist: str, orientation: str,
                  original_start: int, original_end: int, jitter: int, base_quality: int) -> Tuple[str, int, int]:
        """
        Find best matching barcode using precomputed exact and mismatch data
        Returns: Tuple of (barcode, mismatch_count, position)
        Always returns a tuple, using 'null' for no match
        """
        # Parse whitelist identifier
        whitelist = self._parse_whitelist_identifier(whitelist)
        
        # Check if whitelist exists in matchers
        if whitelist not in self.matchers:
            self.logger.warning(f"Whitelist '{whitelist}' not found in available whitelists: {list(self.matchers.keys())}")
            return 'NNNNNNNN', -1, original_start - 1
            
        # Apply quality filtering if needed
        if base_quality is not None:
            sequence = filter_low_quality_bases(sequence, quality_scores, base_quality)[0]
        
        # Generate string for null match of length equal to exact barcode on whitelist
        matcher = self.matchers[whitelist]
        
        # Handle reverse orientation
        if orientation == 'reverse':
            sequence = self._reverse_complement(sequence)
        
        if self.verbose:
            self.logger.debug(f"Looking for barcode in sequence: {sequence}")
            self.logger.debug(f"Using whitelist '{whitelist}' with {len(matcher['exact'])} exact matches")
        
        # Get possible sequences with jitter
        possible_sequences = self._get_sequence_with_jitter(sequence, original_start, original_end, jitter)
        
        # First pass: Look for exact matches with best position
        exact_matches = []
        for seq, pos in possible_sequences:
            if seq in matcher['exact']:
                pos_distance = abs(pos - (original_start - 1))
                exact_matches.append((seq, pos_distance, pos))
                if self.verbose:
                    self.logger.debug(f"Found exact match: {seq} at position {pos}")
        
        # If we found exact matches, return the one closest to the expected position
        if exact_matches:
            exact_matches.sort(key=lambda x: x[1])
            match = exact_matches[0]
            if self.verbose:
                self.logger.debug(f"Selected best exact match: {match[0]} at position {match[2]}")
            return match[0], 0, match[2]
            
        # If no exact match was found, check mismatch lookup
        if self.mismatches > 0:
            for seq, pos in possible_sequences:
                if seq in matcher['mismatches']:
                    original = matcher['mismatches'][seq]
                    if self.verbose:
                        self.logger.debug(f"Found mismatch variant: {seq} -> {original} at position {pos}")
                    return original, 1, pos
        
        # No match found
        if self.verbose:
            self.logger.debug("No match found")
        null_match = "N" * len(next(iter(matcher['exact'])))
        return null_match, -1, original_start - 1

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
        
        for adj_start in range(min_start, max_start):
            seq = full_sequence[adj_start:adj_start + barcode_length]
            if len(seq) == barcode_length:
                sequences.append((seq, adj_start))
                
        return sequences

# Add a new class to track statistics
class BarcodeStats:
    def __init__(self):
        self.mismatch_counts = Counter()
        self.position_counts = Counter()
        
    def update(self, mismatches: List[str], positions: List[str]):
        # Convert mismatches to integers, handling 'NA' values as -1
        mismatch_values = [int(m) if m != 'NA' else -1 for m in mismatches]
        
        # If any mismatch is negative, sum only the negative values
        if any(m < 0 for m in mismatch_values):
            total_mismatches = sum(m for m in mismatch_values if m < 0)
        else:
            total_mismatches = sum(mismatch_values)
            
        self.mismatch_counts[total_mismatches] += 1
        
        # Track each position
        for pos in positions:
            if pos != '0':  # Exclude invalid positions
                self.position_counts[int(pos)] += 1
    
    def write_stats(self, output_prefix: str):
        # Write mismatch counts
        with open(f"{output_prefix}_mismatch_stats.csv", 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['mismatches', 'count'])
            for mismatches, count in sorted(self.mismatch_counts.items()):
                writer.writerow([mismatches, count])
                
        # Write position counts
        with open(f"{output_prefix}_position_stats.csv", 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['position', 'count'])
            for position, count in sorted(self.position_counts.items()):
                writer.writerow([position, count])

def generate_mismatches(sequence: str, max_mismatches: int) -> Dict[int, Dict[str, str]]:
    """Generate all possible sequences with up to max_mismatches mismatches"""
    bases = {'A', 'C', 'G', 'T'}
    mismatch_dict = {i: {} for i in range(1, max_mismatches + 1)}
    seq_length = len(sequence)
    
    for num_mismatches in range(1, max_mismatches + 1):
        for positions in itertools.combinations(range(seq_length), num_mismatches):
            for replacements in itertools.product(bases, repeat=num_mismatches):
                # Convert to list for mutation
                mutated_seq = list(sequence)
                valid_mutation = False
                
                # Apply mutations
                for pos, new_base in zip(positions, replacements):
                    if mutated_seq[pos] != new_base:  # Ensure we introduce a mismatch
                        mutated_seq[pos] = new_base
                        valid_mutation = True
                
                # Only add if we actually introduced a mismatch
                if valid_mutation:
                    mutated_str = ''.join(mutated_seq)
                    if mutated_str not in mismatch_dict[num_mismatches]:
                        mismatch_dict[num_mismatches][mutated_str] = sequence

    return mismatch_dict
    

def parser_reap(parser):
    subparser = parser.add_parser(
        "reap",
        description="""
Extract sequence range from one fastq file of a pair, and annotate the sequence header with barcode 
sequences based on predicted positions generated by scarecrow harvest.

Example:

scarecrow reap --fastqs R1.fastq.gz R2.fastq.gz\n\t--barcode_positions barcode_positions.csv\n\t--barcodes\tBC1:v1_whitelist:bc1_whitelist.txt\n\t\t\tBC2:v2_whitelist:bc2_whitelist.txt\n\t\t\tBC3:v1_whitelist:bc3_whitelist.txt\n\t--read1 0-64 \n\t--out extracted_sequences \n\t--out_sam
---
""",
        epilog="The --barcodes <name> must match the barcode_whitelist values in the --barcode_positions file.",
        help="Extract sequence range from fastq files",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument(
        "--fastqs", 
        nargs="+", 
        required=True,
        help="Pair of FASTQ files")
    subparser.add_argument(
        "-o", "--out",
        metavar="<file>",
        help=("Prefix for output"),
        type=str,
        default='extracted',
    )
    # Add a mutually exclusive group
    out_format = subparser.add_mutually_exclusive_group(required = False)
    out_format.add_argument(
        "--out_sam", 
        action='store_true',
        help='Write output to SAM format [default]'
    )
    out_format.add_argument(
        "--out_fastq", 
        action='store_true',
        help='Write output to FASTQ format'
    )    
    subparser.add_argument(
        "-p", "--barcode_positions",
        metavar="<file>",
        help=("File containing barcode positions, output by scarecrow harvest"),
        type=str,
        required=True,
        default=[],
    )
    subparser.add_argument(
        "-j", "--jitter",
        metavar="<int>",
        type=int,
        default=5,
        help='Barcode position jitter [3]',
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
        required=True,
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
        required=True,
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
             FASTQ = args.out_fastq,
             SAM = args.out_sam,
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
             jitter: int = 3,
             mismatches: int = 1,
             base_quality: int = None,
             batches: int = 10000,
             threads: int = 4,
             FASTQ: bool = False,
             SAM: bool = True,
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

    # Default output to SAM format if not specified
    if FASTQ is False and SAM is False:
        SAM = True

    # Append suffix to output prefix
    logger.info(f"FASTQ = {FASTQ} SAM = {SAM}")
    if FASTQ:
        outfile = f'{output}.fastq'        
    if SAM:
        outfile = f'{output}.sam'
    logger.info(f"Results will be written to {outfile}")

    # Extract sequences
    extract_sequences(
        fastq_files = [f for f in fastqs],
        barcode_positions_file = barcode_positions,
        barcode_reverse_order = barcode_reverse_order,
        barcode_files = barcode_files,
        output = outfile,
        extract = extract,
        umi = umi,
        jitter = jitter,
        mismatches = mismatches,
        base_quality = base_quality,
        batch_size = batches,
        threads = threads,
        FASTQ = FASTQ,
        SAM = SAM,
        verbose = verbose
    )
    
    # gzip
    if FASTQ and gzip:
        logger.info(f"Compressing '{outfile}'")
        with open(outfile, 'rb') as f_in, gz.open(outfile + ".gz", 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(outfile)

@log_errors
def process_read_batch(read_batch: List[Tuple], 
                       barcode_configs: List[Dict],
                       matcher: BarcodeMatcherOptimized,
                       read_range: Tuple[int, int],
                       read_index: int,
                       umi_index: int,
                       umi_range: Tuple[int, int],
                       stats: BarcodeStats,
                       base_quality: int = None,
                       jitter: int = 3,
                       FASTQ: bool = False,
                       SAM: bool = True,
                       verbose: bool = False) -> List[str]:
    """
    Modified process_read_batch to handle jitter and improved matching.
    Now includes original barcodes (CR), barcode qualities (CY), 
    matched barcodes (CB), UMI sequence (UR), and UMI qualities (UY).
    Now tracks barcode statistics.
    """
    logger = setup_worker_logger()
    if verbose:
        logger.info(f"Processing batch of {len(read_batch)} reads")
    output_entries = []
    read_count = len(read_batch)
    
    for i, reads in enumerate(read_batch):
        if i % 100000 == 0:
            logger.debug(f"Processing read {i}/{read_count}")
            
        original_barcodes = []  # CR
        barcode_qualities = []  # CY
        matched_barcodes = []  # CB
        positions = []
        mismatches = []

        for config in barcode_configs:
            seq = reads[config['file_index']].sequence
            qual = reads[config['file_index']].quality
            start, end = config['start'], config['end']
            
            # Extract original barcode and its quality scores
            original_barcode = seq[start-1:end]
            barcode_quality = qual[start-1:end]
            original_barcodes.append(original_barcode)
            barcode_qualities.append(barcode_quality)

            whitelist = ast.literal_eval(config['whitelist'])[1]
            if isinstance(whitelist, list) and len(whitelist) == 1 and isinstance(whitelist[0], tuple):
                whitelist = whitelist[0]

            if verbose:
                logger.info(f"Read: {reads[config['file_index']].name} {reads[config['file_index']].comment}")
                logger.info(f"Sequence: {seq}")
                logger.info(f"Looking for barcode in range {start}-{end} with jitter {jitter}")

            if whitelist in matcher.matchers.keys():                
                matched_barcode, mismatch_count, adj_position = matcher.find_match(
                    seq, reads[config['file_index']].quality, whitelist, config['orientation'], 
                    start, end, jitter, base_quality)

                if verbose:
                    logger.info(f"Matched barcode: {matched_barcode} with {mismatch_count} mismatches at position {adj_position}")
                
                matched_barcodes.append(matched_barcode)
                positions.append(str(adj_position + 1))  # Convert to 1-based position
                mismatches.append(str(mismatch_count))
            else:
                logger.warning(f"Whitelist {whitelist} not found in matcher")
                matched_barcodes.append('null')
                positions.append(str(start))
                mismatches.append('NA')

        # Update statistics
        stats.update(mismatches, positions)

        # Create output entry
        source_entry = reads[read_index]
        
        # Apply base quality filtering to extracted sequence
        if base_quality is not None:
            filtered_seq, filtered_qual = filter_low_quality_bases(
                source_entry.sequence[read_range[0]:read_range[1]], 
                source_entry.quality[read_range[0]:read_range[1]], base_quality)
        else:
            filtered_seq = source_entry.sequence[read_range[0]:read_range[1]]
            filtered_qual = source_entry.quality[read_range[0]:read_range[1]]
        
        #  Output FASTQ
        if FASTQ:
            # Build the header with all components
            header = f"@{source_entry.name} {source_entry.comment}"
            header += f" CR={('_').join(original_barcodes)}"
            header += f" CY={('_').join(barcode_qualities)}"
            header += f" CB={('_').join(matched_barcodes)}"
            header += f" XP={('_').join(positions)}"
            header += f" XM={('_').join(mismatches)}"

            # Add UMI information if specified
            if umi_index is not None:
                umi_seq = reads[umi_index].sequence[umi_range[0]:umi_range[1]]
                umi_qual = reads[umi_index].quality[umi_range[0]:umi_range[1]]
                header += f" UR={umi_seq}"
                header += f" UY={umi_qual}"

            output_entries.append(f"{header}\n{filtered_seq}\n+\n{filtered_qual}\n")
        
        # Output SAM
        if SAM:
            qname = source_entry.name  # QNAME
            flag = "4"  # FLAG (4 means unaligned)
            rname = "*"  # RNAME (unaligned, so *)
            pos = "0"  # POS (unaligned, so 0)
            mapq = "255"  # MAPQ (255 means unavailable)
            cigar = "*"  # CIGAR (unaligned, so *)
            rnext = "*"  # RNEXT (unaligned, so *)
            pnext = "0"  # PNEXT (unaligned, so 0)
            tlen = "0"  # TLEN (unaligned, so 0)
            seq = filtered_seq  # SEQ
            qual = filtered_qual  # QUAL

            # Optional tags
            tags = []
            tags.append(f"CR:Z:{'_'.join(original_barcodes)}")
            tags.append(f"CY:Z:{'_'.join(barcode_qualities)}")
            tags.append(f"CB:Z:{'_'.join(matched_barcodes)}")
            tags.append(f"XP:Z:{'_'.join(positions)}")
            tags.append(f"XM:Z:{'_'.join(mismatches)}")

            if umi_index is not None:
                umi_seq = reads[umi_index].sequence[umi_range[0]:umi_range[1]]
                umi_qual = reads[umi_index].quality[umi_range[0]:umi_range[1]]
                tags.append(f"UR:Z:{umi_seq}")
                tags.append(f"UY:Z:{umi_qual}")

            # Combine all fields into a SAM line
            sam_line = "\t".join([qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual] + tags) + "\n"
            output_entries.append(sam_line)

    if verbose:
        logger.info(f"Completed processing batch of {read_count} reads")
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
    jitter: int = 3,
    mismatches: int = 1,
    base_quality: int = None,
    batch_size: int = 100000,
    threads: Optional[int] = None,
    FASTQ: bool = False,
    SAM: bool = True,
    verbose: bool = False
) -> None:
    """
    Modified to track stats
    """
    logger = setup_worker_logger()
    logger.info("Starting sequence extraction")

    # Initialize master stats tracker
    master_stats = BarcodeStats()

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
    logger.info(f"Generating barcode matcher")
    matcher = BarcodeMatcherOptimized(
        barcode_files = barcode_files,
        mismatches = mismatches,
        base_quality_threshold = base_quality,
        verbose = verbose
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
            # Get batch data
            results, batch_count, batch_stats = job.get()
            
            # Update total read count
            total_reads += batch_count  

            # Update master stats
            master_stats.mismatch_counts.update(batch_stats.mismatch_counts)
            master_stats.position_counts.update(batch_stats.position_counts)

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
                                                                 base_quality, jitter, FASTQ, SAM,
                                                                 verbose),)))
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
                                                                 base_quality, jitter, FASTQ, SAM,
                                                                 verbose),)))
            
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

    # Write final statistics after all processing is complete
    master_stats.write_stats(output)
    logger.info(f"Barcode statistics written to {output}_mismatch_stats.csv and {output}_position_stats.csv")


def parse_range(range_str: str) -> Tuple[int, int]:
    """
    Parse range string
    """
    start, end = map(int, range_str.split('-'))
    start = max(0, start -1)
    return (start, end)

def parse_seed_arguments(barcodes: List[str]) -> Dict[str, str]:
    """
    Modified to handle JSON whitelist files
    Format: <barcode_name>:<whitelist_name>:<whitelist_file.json>
    Returns: Dict mapping whitelist names to JSON file paths
    """
    logger = logging.getLogger('scarecrow')
    barcode_files = {}
    
    if barcodes:
        for barcode in barcodes:
            name, whitelist, filename = barcode.split(':')
            logger.info(f"Processing barcode definition: {barcode}")
            logger.info(f"  Name: {name}")
            logger.info(f"  Whitelist: {whitelist}")
            logger.info(f"  File: {filename}")
            
            if not os.path.exists(filename):
                raise FileNotFoundError(f"Barcode whitelist file not found: {filename}")
            barcode_files[whitelist] = filename
    
    logger.info(f"Parsed barcode whitelists: {barcode_files}")
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

def setup_worker_logger(log_file: str = None):
    """Configure logger for worker processes with file output"""
    logger = logging.getLogger('scarecrow')
    if not logger.handlers:  # Only add handlers if none exist
        # Create formatters
        formatter = logging.Formatter('%(asctime)s - %(processName)s - %(name)s - %(levelname)s - %(message)s')
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
        
        # File handler if log_file provided
        if log_file:
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
        
        logger.setLevel(logging.INFO)
    return logger

def worker_task(args):
    """Worker function that ensures logger is configured for the process"""
    logger = setup_worker_logger()
    (batch, barcode_configs, matcher, 
     extract_range, extract_index, 
     umi_index, umi_range, 
     base_quality, jitter, FASTQ, SAM, verbose) = args
    
    stats = BarcodeStats()  # Create stats tracker for this batch

    if verbose:
        logger.info(f"Worker process started")
    try:
        entries, count = process_read_batch(batch, barcode_configs, matcher,
            extract_range, extract_index, umi_index, umi_range,
            stats, base_quality, jitter, FASTQ, SAM, verbose)
        if verbose:
            logger.info(f"Processed batch of {count} reads")
        return entries, count, stats
    except Exception as e:
        if verbose:
            logger.error(f"Error processing batch: {str(e)}", exc_info=True)
        raise