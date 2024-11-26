import pysam
import os
import json
from typing import List, Dict, Generator, Any
from tqdm import tqdm

def extract_region_details(entry: pysam.FastxRecord, 
                            regions: List) -> List[Dict]:
    """
    Extract region-specific details from a single FASTQ entry
    
    :param entry: pysam FastxRecord
    :param regions: List of RegionCoordinate objects
    :return: List of region details
    """
    # Extract full sequence and quality
    full_sequence = entry.sequence
    full_qualities = entry.quality
    
    # Extract region-specific information
    region_details = []
    for region_coord in regions:
        start = region_coord.start
        stop = region_coord.stop
        
        # Extract sequence for the specific region
        region_sequence = full_sequence[start:stop]
        
        # Extract quality scores for the specific region
        region_qualities = full_qualities[start:stop] if full_qualities else ''
        
        region_details.append({
            'region_id': region_coord.region_id,
            'region_type': region_coord.region_type,
            'name': region_coord.name,
            'sequence_type': region_coord.sequence_type,
            'sequence': region_sequence,
            'qualities': region_qualities,
            'start': start,
            'stop': stop,
            'min_len': region_coord.min_len,
            'max_len': region_coord.max_len
        })
    
    return region_details

def batch_process_paired_fastq(fastq_info: List[Dict], 
                                batch_size: int = 10000, 
                                max_batches: int = None) -> Generator[List[Dict], None, None]:
    """
    Process paired-end FASTQ files in memory-efficient batches.
    
    :param fastq_info: List of dictionaries with paired FASTQ file paths and RegionCoordinate objects
    :param batch_size: Number of read pairs to process in each batch
    :param max_batches: Optional limit on number of batches to process
    :yield: Batches of extracted region information for read pairs
    """
    # Prepare file paths and regions
    r1_file_path = list(fastq_info[0].keys())[0]
    r1_regions = fastq_info[0][r1_file_path]
    r1_strand = fastq_info[0]['strand']
    
    r2_file_path = list(fastq_info[1].keys())[0]
    r2_regions = fastq_info[1][r2_file_path]
    r2_strand = fastq_info[1]['strand']
    
    try:
        # Open both FASTQ files
        with pysam.FastxFile(r1_file_path) as r1_fastq, \
             pysam.FastxFile(r2_file_path) as r2_fastq:
            
            batch = []
            batch_count = 0
            
            while True:
                # Try to read a pair of reads
                try:
                    r1_entry = next(r1_fastq)
                    r2_entry = next(r2_fastq)
                except StopIteration:
                    # Reached end of file
                    break
                
                # Extract region details for both reads
                r1_region_details = extract_region_details(r1_entry, r1_regions)
                r2_region_details = extract_region_details(r2_entry, r2_regions)
                
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
                
                # Add to batch
                batch.append(read_pair_info)
                
                # Yield batch when it reaches batch size
                if len(batch) >= batch_size:
                    yield batch
                    batch = []
                    batch_count += 1
                    
                    # Optional: limit number of batches processed
                    if max_batches and batch_count >= max_batches:
                        break
            
            # Yield any remaining reads in the final batch
            if batch:
                yield batch
    
    except IOError as e:
        print(f"Error reading files: {e}")

def process_paired_fastq_batches(fastq_info: List[Dict], 
                                  batch_size: int = 10000, 
                                  max_batches: int = None, 
                                  output_file: str = None,
                                  region_ids: List = None) -> List[Dict]:
    """
    Process paired-end FASTQ files in batches with optional file output.
    
    :param fastq_info: List of dictionaries with paired FASTQ file paths and RegionCoordinate objects
    :param batch_size: Number of read pairs to process in each batch
    :param max_batches: Optional limit on number of batches to process
    :param output_file: Optional file to write results
    :return: List of processed batches
    """
    all_results = []
    
    # Open output file if specified
    output_handler = open(output_file, 'w') if output_file else None
    
    try:
        # Process batches
        for batch_index, batch in tqdm(enumerate(batch_process_paired_fastq(
            fastq_info, 
            batch_size = batch_size, 
            max_batches = max_batches
        )), desc="Progress", total = max_batches, unit=""):
            # Optional: process or store each batch
            #print(f"Processing batch {batch_index + 1}")
            
            # Write to file if output_file is specified
            if output_handler:
                write_cDNA(batch, region_ids, output_handler)
                #json.dump(batch, output_handler)
                #output_handler.write('\n')  # Newline between batches for easier parsing
            
            # Accumulate results (optional, can be memory intensive for large files)            
            all_results.extend(batch)

    
    except Exception as e:
        print(f"Error processing batches: {e}")
    
    finally:
        # Close output file
        if output_handler:
            output_handler.close()
    
    return all_results




def extract_sequences(data, region_ids=None, verbose=False):
    """
    Extract sequences for specified region_ids from read1 and read2
    
    Args:
    data (list or None): JSON data containing read1 and read2
    region_ids (list): List of region_ids to extract sequences for
    verbose (bool): Print additional information about extraction process
    
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
    not_found_regions = list(region_ids)
    
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
                    not_found_regions = [r for r in not_found_regions if r != region_id]
                    
                    # Store the sequence
                    sequences[region_id] = {
                        'sequence': region.get('sequence', ''),
                        'qualities': region.get('qualities', ''),
                        'read': read_key
                    }
        
        except Exception as e:
            if verbose:
                print(f"Error processing {read_key}: {e}")
    
    # Optionally warn about regions not found
    if not_found_regions and verbose:
        print("The following region_ids were not found:")
        for region in not_found_regions:
            print(region)
    
    return sequences



# Defensive wrapper function for additional safety
def safe_extract_sequences(data, region_ids=None, verbose=False):
    """
    Safely extract sequences with additional error handling
    """
    try:
        return extract_sequences(data, region_ids, verbose)
    except Exception as e:
        if verbose:
            print(f"Unexpected error in extraction: {e}")
        return {}


def write_cDNA(batch, region_ids, output_handler):
    # Generate cDNA fastq output
    header = []
    header.append(batch[0]['read1']['header'])
    #region_ids = ['UMI', 'Round_1_BC', 'Round_2_BC', 'Round_3_BC']
    extracted_sequences = safe_extract_sequences(batch, region_ids, verbose=True)
    for region_id in region_ids:
        if region_id in extracted_sequences:
            info = extracted_sequences[region_id]
            #print(f"\033[0m{region_id} \t\033[32m{info['sequence']} \t\033[0m({info['read']})")
            header.append(info['sequence'])
    delim="_"
    output_handler.write("\n@")
    output_handler.write(delim.join([str(element) for element in header]))
    output_handler.write("\n")
    output_handler.write(safe_extract_sequences(batch, ['cDNA'], verbose=True)['cDNA']['sequence'])
    output_handler.write("\n+\n")
    output_handler.write(safe_extract_sequences(batch, ['cDNA'], verbose=True)['cDNA']['qualities'])
    output_handler.write("\n")

            