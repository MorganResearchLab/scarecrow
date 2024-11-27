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
    :param region_ids: List of region IDs to extract
    :return: List of processed batches
    """
    #all_results = []

    print(f"\033[32m\nProcessing cDNA and writing to {output_file}\033[0m")
    
    # Track unique regions not found across all batches
    not_found_regions = set()

    # Track total read pairs processed
    total_read_pairs = 0

    # Track unique barcodes and their counts
    barcode_counts = {}
    
    # Open output file if specified
    output_handler = open(output_file, 'w') if output_file else None
    output_handler_barcodes = open(output_file + ".barcode.counts", 'w') if output_file else None
    
    try:
        # Create a generator to count total batches if max_batches not specified
        def count_total_batches(fastq_info, batch_size):
            with pysam.FastxFile(list(fastq_info[0].keys())[0]) as r1_fastq, \
                 pysam.FastxFile(list(fastq_info[1].keys())[0]) as r2_fastq:
                try:
                    total = sum(1 for _ in zip(r1_fastq, r2_fastq)) // batch_size
                except Exception:
                    total = None
            return total
        
        # Determine total batches
        total_batches = max_batches if max_batches is not None else count_total_batches(fastq_info, batch_size)

        # Process batches
        for batch_index, batch in tqdm(enumerate(batch_process_paired_fastq(
            fastq_info, 
            batch_size = batch_size, 
            max_batches = max_batches
        )), desc="Progress: ", total=total_batches, unit="%", unit_scale=True):
            # Write to file if output_file is specified
            if output_handler:
                # Pass the not_found_regions tracker to extract_sequences
                safe_extract_sequences(batch, region_ids, verbose = False, not_found_tracker = not_found_regions)

                # Custom write function to track barcodes
                for read_pair in batch:
                    # Prepare header with region sequences
                    header = [read_pair['read1']['header']]
                    extracted_sequences = safe_extract_sequences([read_pair], region_ids)
                    for region_id in region_ids:
                        if region_id in extracted_sequences:
                            header.append(extracted_sequences[region_id]['sequence'])
                    
                    # Create barcode string
                    delim = "_"
                    barcode = delim.join([str(element) for element in header[1:]])
                    
                    # Track barcode count
                    barcode_counts[barcode] = barcode_counts.get(barcode, 0) + 1
                    
                    # Write to file
                    output_handler.write(f"@{delim.join([str(element) for element in header])}\n")
                    output_handler.write(f"{safe_extract_sequences([read_pair], ['cDNA'])['cDNA']['sequence']}\n")
                    output_handler.write("+\n")
                    output_handler.write(f"{safe_extract_sequences([read_pair], ['cDNA'])['cDNA']['qualities']}\n")
            

            
            # Accumulate results and count read pairs
            #all_results.extend(batch)
            total_read_pairs += len(batch)

    except Exception as e:
        print(f"Error processing batches: {e}")
    
    finally:
        # Close output file
        if output_handler:
            output_handler.close()
    
    # Report unique regions not found after processing all batches
    if not_found_regions:
        print(f"\033[32m\nUnique regions not found across all batches\033[0m")
        for region in not_found_regions:
            print(f"{region}")

    # Print total read pairs processed
    print(f"\033[32m\nTotal read pairs processed: \033[0m{total_read_pairs}")

    # Create binned ASCII histogram of barcode distribution
    print("\033[32m\nBarcode Distribution:\033[0m")
    create_binned_ascii_histogram(barcode_counts)

    # Output barcode counts to file
    for barcode, count in sorted(barcode_counts.items(), key=lambda x: x[1], reverse=True):
        output_handler_barcodes.write(f"{barcode},{count}")

    return barcode_counts



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
    print("\nBarcode Count Distribution:")
    print("Bin Ranges (Read Counts) | Histogram")
    
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
    print(f"\nTotal unique barcodes: {len(counts)}")
    print(f"Min count: {min_count}")
    print(f"Max count: {max_count}")    