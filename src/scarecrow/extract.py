#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
"""

import pysam
from seqspec import Assay
from seqspec.Region import RegionCoordinate
from seqspec.utils import load_spec
from seqspec.seqspec_print import run_seqspec_print
from seqspec.seqspec_index import get_index_by_primer, format_kallisto_bus
from argparse import RawTextHelpFormatter
from typing import List, Dict, Set, Optional
import logging
from scarecrow.logger import log_errors, setup_logger

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
        help="Extract sequence element from fastqs",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument("yaml", help="Sequencing specification yaml file")
    subparser.add_argument("fastqs", nargs="+", help="List of FASTQ files")
    subparser.add_argument(
        "-o", "--out",
        metavar="out.fastq",
        help=("Path to output fastq file"),
        type=str,
        default="./cDNA.fq",
    )
    subparser.add_argument(
        "-r", "--header_regions",
        metavar="region_id",
        help=("List of elements to include in sequence header"),
        nargs="*",
        type=str,
        default=[],
    )
    subparser.add_argument(
        "-t", "--target",
        metavar="target",
        help=("Target element to extract sequence of [\"cdna\"]"),
        type=str,
        default="cdna",
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

def validate_extract_args(parser, args):
    run_extract(yaml = args.yaml, fastqs = [f for f in args.fastqs], target = args.target,
                output_file = args.out, batches = args.batch_size, regions = args.header_regions, 
                logfile = args.logfile, threads = args.threads, max_batches = args.max_batches)

def run_extract(yaml, fastqs, output_file, target, batches, max_batches, regions, threads, logfile):
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
                                 barcodes = None, target = target)

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


@log_errors
def process_paired_fastq_batches(
    fastq_info: List[Dict],
    batch_size: int = 10000,
    max_batches: Optional[int] = None,
    output_handler: Optional[str] = None,
    region_ids: Optional[List[str]] = None,
    num_workers: int = None,
    barcodes: Dict = None,
    target: str = None
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
    logger = logging.getLogger('scarecrow')

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
    not_found_regions = set()
    
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
                            # Sequence extraction and barcode tracking
                            if barcodes is None:
                                extracted_sequences = safe_extract_sequences(
                                    [read_pair],
                                    region_ids,
                                    not_found_tracker=not_found_regions
                                )
                                
                                # Barcode key generation
                                barcode_key = []
                                for region_id in region_ids:
                                    if region_id in extracted_sequences:
                                        barcode_key.append(
                                            extracted_sequences[region_id]['sequence']
                                        )
                                barcode_key = "_".join(
                                    [str(element) for element in barcode_key]
                                )
                                
                                # Update barcode counts
                                barcode_counts[barcode_key] = barcode_counts.get(
                                    barcode_key, 0
                                ) + 1
                                
                                # Optional output handling
                                if output_handler:
                                    write_cDNA_fastq(
                                        read_pair, 
                                        extracted_sequences, 
                                        region_ids, 
                                        target, 
                                        output_handler
                                    )
                            else:
                                if output_handler:
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
def safe_extract_sequences(data, region_ids=None, verbose=False, not_found_tracker=None):
    """
    Safely extract sequences with additional error handling
    """
    logger = logging.getLogger('scarecrow')
    try:
        return extract_sequences(data, region_ids, verbose, not_found_tracker)
    except Exception as e:
        if verbose:
            print(f"Unexpected error in extraction: {e}")
        return {}
    
@log_errors
def extract_sequences(data, region_ids=None, verbose=False, not_found_tracker=None):
    """
    Extract sequences for specified region_ids from read1 and read2
    
    Args:
    data (list or None): JSON data containing read1 and read2
    region_ids (list): List of region_ids to extract sequences for
    verbose (bool): Print additional information about extraction process
    not_found_tracker (set): Optional set to track regions not found
    
    Returns:
    dict: Extracted sequences
    """
    logger = logging.getLogger('scarecrow')

    # Input validation
    if data is None:
        if verbose:
            print("Error: Input data is None")
        return {}
    
    # Ensure region_ids is a list
    if region_ids is None:
        region_ids = []
    
    # Initialize results dictionary
    sequences = {}
    
    # Tracks which region_ids were not found
    batch_not_found_regions = set(region_ids)
    
    # Ensure data is a list and has at least one element
    if not isinstance(data, list) or len(data) == 0:
        if verbose:
            print("Error: Data is not a non-empty list")
        return {}
    
    # Iterate through reads (read1 and read2)
    for read_key in ['read1', 'read2']:
        try:
            # Get regions for current read
            regions = data[0].get(read_key, {}).get('regions', [])
            
            # Find sequences for specified region_ids
            for region in regions:
                region_id = region.get('region_id')
                
                # Check if this region is in our desired region_ids
                if region_id in region_ids:
                    # Remove from not found list
                    if region_id in batch_not_found_regions:
                        batch_not_found_regions.remove(region_id)
                    
                    # Store the sequence
                    sequences[region_id] = {
                        'sequence': region.get('sequence', ''),
                        'qualities': region.get('qualities', ''),
                        'read': read_key
                    }
        
        except Exception as e:
            if verbose:
                print(f"Error processing {read_key}: {e}")
    
    # If tracking not found regions, update the tracker
    if not_found_tracker is not None:
        not_found_tracker.update(batch_not_found_regions)
    
    # Optionally warn about regions not found in this batch
    if batch_not_found_regions and verbose:
        print("Regions not found in this batch:")
        for region in batch_not_found_regions:
            print(region)
    
    return sequences


@log_errors
def extract_region_seq(
    entry: pysam.FastxRecord,
    regions: Optional[List[RegionCoordinate]] = None,
    rev: bool = False
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
        if regions is None:
            raise ValueError("Must provide region details for extraction")
        
        # Prepare sequence
        full_sequence = entry.sequence
        full_qualities = entry.quality
        
        # Extract by regions if provided
        region_details = []
        for region_coord in regions:
            # Validate region coordinates
            if region_coord.start < 0 or region_coord.stop > len(full_sequence):
                logger.warning(f"Invalid region coordinates for {region_coord.name}")
                continue
            # Reverse sequence to align with region positions from seqspec
            if rev:
                full_sequence = full_sequence[::-1]
                full_qualities = full_qualities[::-1]                       
            region_sequence = full_sequence[region_coord.start:region_coord.stop]
            region_qualities = full_qualities[region_coord.start:region_coord.stop] if full_qualities else ''
            
            region_details.append({
                'region_id': region_coord.region_id,
                'region_type': region_coord.region_type,
                'name': region_coord.name,
                'sequence_type': region_coord.sequence_type,
                'sequence': region_sequence,
                'qualities': region_qualities,
                'start': region_coord.start,
                'stop': region_coord.stop,
                'min_len': region_coord.min_len,
                'max_len': region_coord.max_len
            })
        
        return region_details
                
    except Exception as e:
        logger.error(f"Error in extract_region_details: {e}")
        raise


@log_errors
def write_cDNA_fastq(read_pair: Dict, extracted_sequences: Dict, region_ids: List, target: str, output_handler):
    """
    Write processed read pair to output file.
    
    Args:
        read_pair (Dict): Read pair information
        extracted_sequences (Dict): Extracted sequences
        output_file (str): Output file path
    """
    logger = logging.getLogger('scarecrow')

    # Implementation depends on specific output requirements
    header = [read_pair['read1']['header']]
    for region_id in region_ids:
        if region_id in extracted_sequences:
            header.append(extracted_sequences[region_id]['sequence'])
    # Extract cDNA
    """
    An issue here is that it is looking for 'cdna' so if the region is labeled cDNA for example it will not find it.
    Solution is to lowercase the region, probably when running region_indices in scarecrow_extract.py
    """
    cdna_seq = {safe_extract_sequences([read_pair], [target])[target]['sequence']}
    cdna_qual = {safe_extract_sequences([read_pair], [target])[target]['qualities']}
    # Write to file
    output_handler.write(f"@{"_".join([str(element) for element in header])}\n")
    output_handler.write(f"{cdna_seq}\n")
    output_handler.write("+\n")
    output_handler.write(f"{cdna_qual}\n")
