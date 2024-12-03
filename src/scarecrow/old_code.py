#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
"""

from scarecrow.fastq_logging import logger, log_errors
import pysam
import os
from typing import List, Dict, Generator, Optional, Set, Union
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
import resource
from seqspec.Region import RegionCoordinate


class FastqProcessingError(Exception):
    """Custom exception for FASTQ processing errors."""
    pass

def get_memory_usage() -> float:
    """
    Get current memory usage of the process.
    
    Returns:
        float: Memory usage in MB
    """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024

def reverse_complement(seq):
    """
    Short function to reverse complement a sequence
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

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
def parse_barcode_arguments(barcode_args):
    """
    Parse barcode arguments from command line.
    
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
def batch_process_paired_fastq(
    fastq_info: List[Dict], 
    batch_size: int = 10000, 
    max_batches: Optional[int] = None,
    barcodes: List[Dict] = None,
) -> Generator[List[Dict], None, None]:
    """
    Process paired-end FASTQ files in memory-efficient batches.
    
    Args:
        fastq_info (List[Dict]): FASTQ file paths and region information
        batch_size (int): Number of read pairs per batch
        max_batches (int, optional): Limit on number of batches
    
    Yields:
        List[Dict]: Batches of read pair information
    
    Raises:
        FastqProcessingError: If file reading fails
    """
    logger.info(f"Starting batch processing. Initial memory: {get_memory_usage():.2f} MB")
    
    try:
        # Unpack FASTQ file information
        r1_file_path = list(fastq_info[0].keys())[0]
        r1_regions = fastq_info[0][r1_file_path]
        r1_strand = fastq_info[0]['strand']
        
        r2_file_path = list(fastq_info[1].keys())[0]
        r2_regions = fastq_info[1][r2_file_path]
        r2_strand = fastq_info[1]['strand']
        
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
                
                # Extract region details
                r1_region_details = extract_region_details(r1_entry, r1_regions, rev = False, expected_barcodes = barcodes)
                r2_region_details = extract_region_details(r2_entry, r2_regions, rev = True, expected_barcodes = barcodes)
                
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
                    logger.info(f"Batch {batch_count + 1} memory: {get_memory_usage():.2f} MB")
                    yield batch
                    batch = []
                    batch_count += 1
                    
                    # Optional batch limit
                    if max_batches and batch_count >= max_batches:
                        break
            
            # Yield final batch
            if batch:
                yield batch
        
    except IOError as e:
        logger.error(f"File reading error: {e}")
        raise FastqProcessingError(f"Unable to process FASTQ files: {e}")
    
@log_errors
def process_paired_fastq_batches(
    fastq_info: List[Dict], 
    batch_size: int = 10000, 
    max_batches: Optional[int] = None,
    output_file: Optional[str] = None,
    region_ids: Optional[List[str]] = None,
    num_workers: int = 4,
    barcodes: Union[Dict[str, List[str]], List[str], None] = None,
) -> Dict[str, int]:
    """
    Process paired-end FASTQ files with parallel and memory-efficient processing.
    
    Args:
        fastq_info (List[Dict]): FASTQ file information
        batch_size (int): Number of read pairs per batch
        max_batches (int, optional): Limit on number of batches
        output_file (str, optional): Output FASTQ file path
        region_ids (List[str], optional): Regions to extract
        num_workers (int): Number of parallel processing workers
    
    Returns:
        Dict[str, int]: Barcode count distribution
    """
    barcode_counts = {}
    total_read_pairs = 0
    not_found_regions = set()

    logger.info(f"Expected barcodes")
    for barcode in barcodes:
        logger.info(f"{barcode}: {barcodes[barcode]}")

    try:
        # Parallel processing setup
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            # Prepare batches
            batches = list(batch_process_paired_fastq(
                fastq_info, 
                batch_size=batch_size, 
                max_batches=max_batches,
                barcodes=barcodes
            ))
            
            if output_file:
                print(f"\033[32m\nProcessing cDNA and writing to {output_file}\033[0m")
                output_handler = open(output_file, 'w') 
            else:
                output_handler = None
            
            # Progress tracking
            with tqdm(total=len(batches), desc="Processing Batches") as pbar:
                for batch in batches:
                    # Process batch
                    for read_pair in batch:
                        try:
                            # Extract sequences
                            extracted_sequences = safe_extract_sequences(
                                [read_pair], 
                                region_ids, 
                                not_found_tracker=not_found_regions
                            )
                            
                            # Barcode tracking logic                            
                            barcode_key = []
                            for region_id in region_ids:
                                if region_id in extracted_sequences:
                                    barcode_key.append(extracted_sequences[region_id]['sequence'])
                            barcode_key = "_".join([str(element) for element in barcode_key])
                            barcode_counts[barcode_key] = barcode_counts.get(barcode_key, 0) + 1
                            
                            # Optional file writing
                            if output_file and barcodes is None:
                                write_output(read_pair, extracted_sequences, region_ids, output_handler)

                        except Exception as e:
                            logger.error(f"Error processing read pair: {e}")
                    
                    total_read_pairs += len(batch)
                    pbar.update(1)
        
        # Report and visualization
        if barcodes is None:
            report_processing_results(
                total_read_pairs, 
                not_found_regions, 
                barcode_counts, 
                output_file
            )

    except Exception as e:
        logger.error(f"Processing failed: {e}")

    return barcode_counts

def report_processing_results(
    total_read_pairs: int, 
    not_found_regions: Set[str], 
    barcode_counts: Dict[str, int], 
    output_file: Optional[str]
):
    """
    Generate comprehensive processing report.
    
    Args:
        total_read_pairs (int): Total number of read pairs processed
        not_found_regions (Set[str]): Regions not found during processing
        barcode_counts (Dict[str, int]): Barcode distribution
        output_file (str, optional): Output file path
    """
    logger.info(f"Total read pairs processed: {total_read_pairs}")
    
    if not_found_regions:
        logger.warning("Regions not found:")
        for region in not_found_regions:
            logger.warning(region)
    
    # Convert barcode distribution to log-friendly format
    barcode_summary = sorted(barcode_counts.items(), key=lambda x: x[1], reverse=True)
    logger.info("Barcode Distribution (10 most frequent umi + barcode combinations):")
    for barcode, count in barcode_summary[:10]:  # Log top 10 barcodes
        logger.info(f"Barcode: {barcode}, Count: {count}")
    
    logger.info(f"Total unique UMI + barcode combinations: {len(barcode_counts)}")
    logger.info(f"Min barcode count: {min(barcode_counts.values()) if barcode_counts else 0}")
    logger.info(f"Max barcode count: {max(barcode_counts.values()) if barcode_counts else 0}")
    
    # Write barcode counts to file
    if output_file:
        with open(f"{output_file}.barcode_counts.csv", 'w') as f:
            f.write("umi_barcodes,Count\n")  # CSV header
            for barcode, count in sorted(barcode_counts.items(), key=lambda x: x[1], reverse=True):
                f.write(f"{barcode},{count}\n")
    
    # Print total read pairs processed to console
    print(f"\033[32m\nTotal read pairs processed\033[0m: {total_read_pairs}")

    # Create binned ASCII histogram of barcode distribution and output to console
    create_binned_ascii_histogram(barcode_counts)

@log_errors
def write_output(read_pair: Dict, extracted_sequences: Dict, region_ids: List, output_handler):
    """
    Write processed read pair to output file.
    
    Args:
        read_pair (Dict): Read pair information
        extracted_sequences (Dict): Extracted sequences
        output_file (str): Output file path
    """
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
    cdna_seq = {safe_extract_sequences([read_pair], ['cdna'])['cdna']['sequence']}
    cdna_qual = {safe_extract_sequences([read_pair], ['cdna'])['cdna']['qualities']}
    # Write to file
    output_handler.write(f"@{"_".join([str(element) for element in header])}\n")
    output_handler.write(f"{cdna_seq}\n")
    output_handler.write("+\n")
    output_handler.write(f"{cdna_qual}\n")

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
def safe_extract_sequences(data, region_ids=None, verbose=False, not_found_tracker=None):
    """
    Safely extract sequences with additional error handling
    """
    try:
        return extract_sequences(data, region_ids, verbose, not_found_tracker)
    except Exception as e:
        if verbose:
            print(f"Unexpected error in extraction: {e}")
        return {}

def create_binned_ascii_histogram(counts, num_bins=10):
    """
    Create a binned ASCII histogram of barcode counts
    
    :param counts: Dictionary of counts to visualize
    :param num_bins: Number of bins to use for distribution
    """
    # If no counts, return early
    if not counts:
        print("No data to create histogram.")
        return
    
    # Calculate bin boundaries
    count_values = list(counts.values())
    min_count = min(count_values)
    max_count = max(count_values)
    
    # Create logarithmic bins (log scale is often better for count distributions)
    import numpy as np
    bins = np.logspace(np.log10(max(1, min_count)), np.log10(max_count), num=num_bins+1)
    
    # Bin the counts
    binned_counts = [0] * num_bins
    for count in count_values:
        # Find which bin the count belongs to
        bin_index = np.digitize(count, bins) - 1
        # Ensure it's within the last bin if at the edge
        bin_index = min(bin_index, num_bins - 1)
        binned_counts[bin_index] += 1
    
    # Find max for scaling
    max_bin_count = max(binned_counts)
    
    # Create ASCII histogram
    print("\033[32m\nBarcode Distribution\033[0m")
    print("\033[34mBin Ranges \033[0m(\033[34mRead Counts\033[0m) | Histogram")
    
    # Height of the histogram
    height = 20
    
    for i, (bin_count, bin_left) in enumerate(zip(binned_counts, bins[:-1])):
        # Scale the bar
        bar_length = int((bin_count / max_bin_count) * height)
        bar = '#' * bar_length
        
        # Format bin range (use log10 for readability)
        bin_range = f"{bin_left:.0f}-{bins[i+1]:.0f}"
        
        # Print the histogram line
        print(f"{bin_range:20} | {bar:20} ({bin_count})")
    
    # Print total unique barcodes
    print(f"\n\033[32mTotal unique barcodes\033[0m: {len(counts)}")
    print(f"\033[32mMin count\033[0m: {min_count}")
    print(f"\033[32mMax count\033[0m: {max_count}")            

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
    
    # Iterate through all possible start positions
    for start in range(len(sequence)):
        for barcode in barcodes:
            # Check all possible end positions for this barcode
            for end in range(start + len(barcode), len(sequence) + 1):
                candidate = sequence[start:end]
                
                # Check if candidate matches barcode within max mismatches
                if len(candidate) == len(barcode):
                    mismatches = hamming_distance(candidate, barcode)
                    if mismatches <= max_mismatches:
                        match_details = {
                            'barcode': barcode,
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
def extract_region_details(
    entry: pysam.FastxRecord,
    regions: Optional[List[RegionCoordinate]] = None,
    expected_barcodes: Optional[Union[Dict[str, List[str]], List[str]]] = None,
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
    try:
        # Validate input
        if regions is None and expected_barcodes is None:
            raise ValueError("Must provide either regions or barcodes for extraction")
        
        # Prepare sequence
        full_sequence = entry.sequence
        full_qualities = entry.quality
        
        # Extract by regions if provided
        if regions is not None and expected_barcodes is None:
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
        
        # Find barcodes if regions are not provided
        if expected_barcodes is not None:
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
            
            # Paired-end reads are both 5' to 3', but on opposite strands (Read1 =  forward, Read2 = reverse)
            # 5' ---Read1--->         <---Read2--- 5'
            # So read 2 needs reverse complement as barcodes are on forward strand
            if rev:
                full_sequence = reverse_complement(full_sequence)

            # Find barcode positions
            barcode_matches = find_barcode_positions(full_sequence, barcodes_to_search)

            logger.info(f">{entry.name} {entry.comment}")
            logger.info(f"{full_sequence}")
            
            # Prepare barcode match details
            barcode_details = []
            for match in barcode_matches:
                barcode_detail = {
                    'barcode': match['barcode'],
                    'sequence': match['sequence'],
                    'start': match['start'],
                    'end': match['end'],
                    'mismatches': match['mismatches']
                }
                
                # Add all matching dictionary keys if available
                if barcode_dict_keys and match['barcode'] in barcode_dict_keys:
                    barcode_detail['dict_keys'] = barcode_dict_keys[match['barcode']]
                
                barcode_details.append(barcode_detail)
                
                logger.info(f"...{barcode_detail['dict_keys']} ({match['barcode']}) hit: {match['sequence']}, "
                        f"Start: {match['start']}, "
                        f"End: {match['end']}, "
                        f"Mismatches: {match['mismatches']}")

            
            return barcode_details
        
        # This should never be reached due to initial validation
        raise ValueError("No extraction method selected")
    
    except Exception as e:
        logger.error(f"Error in extract_region_details: {e}")
        raise