import pysam
import os
import json
from typing import List, Dict, Generator, Any


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
                                  output_file: str = None) -> List[Dict]:
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
        for batch_index, batch in enumerate(batch_process_paired_fastq(
            fastq_info, 
            batch_size=batch_size, 
            max_batches=max_batches
        )):
            # Optional: process or store each batch
            print(f"Processing batch {batch_index + 1}")
            
            # Write to file if output_file is specified
            if output_handler:
                json.dump(batch, output_handler)
                output_handler.write('\n')  # Newline between batches for easier parsing
            
            # Accumulate results (optional, can be memory intensive for large files)
            all_results.extend(batch)
    
    except Exception as e:
        print(f"Error processing batches: {e}")
    
    finally:
        # Close output file
        if output_handler:
            output_handler.close()
    
    return all_results