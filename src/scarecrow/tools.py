#!/usr/bin/env python3
"""
@author: David Wragg
"""

from seqspec import Assay
from seqspec.Region import RegionCoordinate
from seqspec.seqspec_index import get_index_by_primer
import os
import resource
from scarecrow.fastq_logging import logger, log_errors
from typing import List, Dict, Optional, Generator, Set, Union
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
import pysam

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
def region_indices(spec: Assay, fastqs):
    """
    Identify library elements contained in sequencing reads
    """
    logger.info(f"Library elements identified by seqspec.get_index_by_primer")
    indices = []
    for fastq in fastqs:
        index = get_index_by_primer(spec, "rna", fastq)
        indices.append(index)
    for index in indices:
        for file, regions in index.items():
            if file in fastqs:
                logger.info(f"{file}")                
                for region in regions:
                    logger.info(f"{region.region_id}: {region.start}-{region.stop}") 

    return indices

@log_errors
def process_paired_fastq_batches(
    fastq_info: List[Dict], 
    batch_size: int = 10000, 
    max_batches: Optional[int] = None,
    output_handler: Optional[str] = None,
    region_ids: Optional[List[str]] = None,
    num_workers: int = 4,
    barcodes: Dict = None,
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

    try:
        # Parallel processing setup
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            # Prepare batches (this is outside the progress bar but takes a long time on a large fastq file)
            batches = list(batch_process_paired_fastq(
                fastq_info, 
                batch_size = batch_size, 
                max_batches = max_batches,
                barcodes = barcodes
            ))
            
            # Open output file for writing
            #output_handler = open(output_file, 'w') if output_file else None


            # Progress tracking
            with tqdm(total=len(batches), desc="Processing Batches") as pbar:
                for batch in batches:
                    # Process batch
                    for read_pair in batch:
                        try:

                            # Logic for processing region element sequence extraction
                            if barcodes is None:

                                # Extract sequences
                                extracted_sequences = safe_extract_sequences(
                                    [read_pair], 
                                    region_ids, 
                                    not_found_tracker=not_found_regions
                                )

                                # Barcode tracking logic for when writing cDNA file
                                barcode_key = []
                                for region_id in region_ids:
                                    if region_id in extracted_sequences:
                                        barcode_key.append(extracted_sequences[region_id]['sequence'])
                                barcode_key = "_".join([str(element) for element in barcode_key])
                                barcode_counts[barcode_key] = barcode_counts.get(barcode_key, 0) + 1
                                if output_handler:
                                    write_cDNA_fastq(read_pair, extracted_sequences, region_ids, output_handler)
                            
                            else:
                                if output_handler:
                                        write_barcodes_CSV(read_pair, output_handler)

                        except Exception as e:
                            logger.error(f"Error processing read pair: {e}")
                    
                    total_read_pairs += len(batch)
                    pbar.update(1)
        
        # Report and visualization
        #if barcodes is None:
        #    report_processing_results(
        #        total_read_pairs, 
        #        not_found_regions, 
        #        barcode_counts, 
        #        output_file
        #    )

    except Exception as e:
        logger.error(f"Processing failed: {e}")

    return barcode_counts

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
def batch_process_paired_fastq(
    fastq_info: List[Dict], 
    batch_size: int = 10000, 
    max_batches: Optional[int] = None,
    barcodes: Dict = None,
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
                
                # Extract region details if barcodes not provided
                if barcodes:
                    r1_region_details = extract_barcodes(r1_entry, rev = False, expected_barcodes = barcodes)
                    r2_region_details = extract_barcodes(r2_entry, rev = True, expected_barcodes = barcodes)
                else:
                    r1_region_details = extract_region_seq(r1_entry, r1_regions, rev = False)
                    r2_region_details = extract_region_seq(r2_entry, r2_regions, rev = True)
                
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
def write_cDNA_fastq(read_pair: Dict, extracted_sequences: Dict, region_ids: List, output_handler):
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
            start = sub_barcode['start']
            end = sub_barcode['end']
            mm = sub_barcode['mismatches']
            output_handler.write(f"{read_key},{read},{bc},{barcode},{start},{end},{mm}\n")
