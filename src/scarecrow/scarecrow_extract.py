#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 09:37:12 2024

@author: s14dw4
"""

from seqspec import Assay
from seqspec.utils import load_spec, map_read_id_to_regions
from argparse import RawTextHelpFormatter
import gzip
from typing import Iterator, Tuple, Union
from pathlib import Path

def parser_extract(parser):
    subparser = parser.add_parser(
        "extract",
        description="""
Extract cDNA sequence from fastq files

Examples:
scarecrow extract spec.yaml fastqs -v -o ~/path/to/output
---
""",
        help="Extract cDNA from fastqs",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument("yaml", help="Sequencing specification yaml file")
    subparser.add_argument("fastqs", nargs="+", help="List of FASTQ files")
    subparser.add_argument(
        "-o",
        metavar="OUT",
        help=("Path to output cDNA fastq files"),
        type=str,
        default=None,
    )
    subparser.add_argument(
        "-p",
        help=("Paired-end sequencing"),
        action='store_true'
    )
    subparser.add_argument(
        "-v",
        help=("Verbose output"),
        action='store_true'
    )
    return subparser


def validate_extract_args(parser, args):
    run_extract(yaml = args.yaml, fastqs = [f for f in args.fastqs], 
                outdir = args.o, paired = args.p, verbose = args.v)


def run_extract(yaml, fastqs, outdir, paired, verbose):
    """
    This function works by:
        1) Going through each modality in the specification
        2) Finding all reads for that modality
        3) Using map_read_id_to_regions to get the associated regions
        4) Collecting the length and sequence information for each region
    """

    # import seqspec.yaml
    spec = load_spec(yaml)
    
    # extract regions from each read
    read_regions = extract_read_regions_info(spec, verbose)
    
    # process fastq files
    process_fastq_files(fastqs, paired, read_regions)


    if outdir:
        # This needs to be edited still
        print("")
    else:
        print("")
    return spec


def extract_read_regions_info(spec: Assay, verbose: bool = True):
    """
    Extract information about regions and their lengths for each Read in the Assay
    """
    results = []
    
    # Iterate through each modality in the assay
    for modality in spec.list_modalities():
        # Get all reads for this modality
        reads = spec.get_seqspec(modality)
        
        for read in reads:
            # Get the regions associated with this read
            try:
                _, regions = map_read_id_to_regions(spec, modality, read.read_id)
                
                # Collect information for this read
                read_info = {
                    'read_id': read.read_id,
                    'modality': modality,
                    'regions': [{
                        'region_id': region.region_id,
                        'region_type': region.region_type,
                        'min_length': region.min_len,
                        'max_length': region.max_len,
                        'sequence': region.sequence
                    } for region in regions]
                }
                results.append(read_info)
                
            except IndexError as e:
                print(f"Warning: Could not map regions for read {read.read_id}: {str(e)}")
    
    if verbose:
        # Print the results in a readable format
        for read in results:
            print(f"\033[32m\nRead: {read['read_id']} (Modality: {read['modality']})\033[0m")
            print("Regions:")
            for region in read['regions']:
                print(f"  - {region['region_id']} ({region['region_type']}):")
                print(f"    Min length: {region['min_length']}")
                print(f"    Max length: {region['max_length']}")
                if region['sequence']:
                    print(f"    Sequence: {region['sequence']}")
        
    return results


def read_fastq_buffered(filename: str, is_paired_end: bool = True, buffer_size: int = 100000) -> Iterator[Union[Tuple[str, str, str], Tuple[Tuple[str, str, str], Tuple[str, str, str]]]]:
    """
    Memory-efficient FASTQ reader supporting single and paired-end reads
    
    Args:
        filename (str): Path to FASTQ file (can be gzipped)
        is_paired_end (bool): Whether the reads are paired-end
        buffer_size (int): Number of bytes to buffer
        
    Yields:
        For single-end: (header, sequence, quality)
        For paired-end: ((header1, seq1, qual1), (header2, seq2, qual2))
    """
    
    # Check if file is gzipped
    is_gzipped = filename.endswith('.gz')
    opener = gzip.open if is_gzipped else open
    
    with opener(filename, 'rt') as f:
        if is_paired_end:
            # For paired-end, read 8 lines at a time
            while True:
                # Read first read
                header1 = f.readline().strip()
                if not header1:  # End of file
                    break
                
                # Ensure header starts with '@'
                if not header1.startswith('@'):
                    continue
                
                seq1 = f.readline().strip()
                plus1 = f.readline().strip()
                qual1 = f.readline().strip()
                
                # Read second read
                header2 = f.readline().strip()
                if not header2.startswith('@'):
                    continue
                
                seq2 = f.readline().strip()
                plus2 = f.readline().strip()
                qual2 = f.readline().strip()
                
                # Yield only if both reads are valid
                if header1 and seq1 and qual1 and header2 and seq2 and qual2:
                    yield ((header1[1:], seq1, qual1), (header2[1:], seq2, qual2))
        
        else:
            # For single-end, read 4 lines at a time
            while True:
                header = f.readline().strip()
                if not header:  # End of file
                    break
                
                # Ensure header starts with '@'
                if not header.startswith('@'):
                    continue
                
                seq = f.readline().strip()
                plus = f.readline().strip()
                qual = f.readline().strip()
                
                # Yield only if all components are valid
                if header and seq and qual:
                    yield (header[1:], seq, qual)
                    

# Example using the buffered reader with batch processing
def process_fastq_in_batches(fastq_file: str, paired: bool, regions_info, batch_size: int = 1000):
    """
    Process FASTQ reads in batches for better performance.
    """
    batch = []
    read_count = 0
    
    print(f"\033[32mProcessing {fastq_file}\033[0m")
    
    for read in read_fastq_buffered(filename = fastq_file, is_paired_end = paired):
        if paired:
            # For paired-end: read is a tuple of two read tuples
            header1, seq1, qual1 = read[0]
            header2, seq2, qual2 = read[1]
            batch.append((header1, seq1, qual1))
            batch.append((header2, seq2, qual2))
        else:
            # For single-end: read is a tuple of (header, seq, qual)
            header, seq, qual = read
            batch.append((header, seq, qual))
        
        read_count += 1
        
        if len(batch) >= batch_size:
            # Process the batch
            process_batch(batch, regions_info)
            batch = []
            
            # Print progress
            print(f"Processed {read_count:,} reads...", end='\r')
    
    # Process any remaining reads
    if batch:
        process_batch(batch, regions_info)
    
    print(f"\nCompleted processing {read_count:,} reads")
       

def process_batch(batch, regions_info):
    """
    Logic function to process a batch of reads.
    """
    for header, seq, qual in batch:
            # Check state of header
            parsed_header = parse_fastq_header(header)
            if not parsed_header['is_standard']:
                print(f"Warning: Non-standard header: {header}")
        
            # Extract each region based on the expected lengths
            current_pos = 0
            extracted_regions = {}
            for region in regions_info['regions']:
                region_length = region['min_length']  # Or use some logic to determine length
                region_sequence = seq[current_pos:current_pos + region_length]
                region_quality = qual[current_pos:current_pos + region_length]
                    
                extracted_regions[region['region_id']] = {
                    'sequence': region_sequence,
                    'quality': region_quality,
                    'type': region['region_type']
                }
                
                current_pos += region_length

            # Do something with the extracted regions
            # For example, print the first read's regions:
            if len(extracted_regions) > 0:
                print(f"\nExample read regions for {header}:")
                for region_id, region_data in extracted_regions.items():
                    print(f"{region_id}: {region_data['sequence']}")
                break  # Just show first read as example
        

def process_fastq_files(fastqs: list, paired, read_regions):
    """
    Process FASTQ reads using their region specifications
    """
    try:
        for file in fastqs:
            regions_info = match_fastq_to_regions(read_regions, Path(file).name)
            print(f"\nProcessing {file}")
            print(f"Modality: {regions_info['modality']}")
            process_fastq_in_batches(file, paired, regions_info, batch_size=1000)                
    
    except ValueError as e:
        print(f"Error: {e}")
    
    except Exception as e:
        print(f"Unexpected error processing file: {e}")


def match_fastq_to_regions(read_regions: list, fastq_file: str) -> dict:
    """
    Check FASTQ file is present among read_regions
    """
    # First find the read_regions entry matching fastq_file
    matching_read = None
    
    for read in read_regions:
        if read['read_id'] == fastq_file:  # or whatever identifier matches your FASTQ file
            matching_read = read
            break    
    
    if not matching_read:
        raise ValueError(f"No matching read configuration found for {fastq_file}")
        
    return matching_read        


def parse_fastq_header(header):
    """
    Safely parse FASTQ headers with fallback mechanism
    """
    # Skip '+' lines or empty headers
    if header == '+' or not header:
        return {
            'original_header': header,
            'is_standard': False
        }
    
    # Remove '@' prefix if present
    if header.startswith('@'):
        header = header[1:]
    
    # Split standard Illumina header
    parts = header.split(':')
    
    # If header doesn't match expected format, return as non-standard
    if len(parts) < 7:
        return {
            'original_header': header,
            'is_standard': False
        }
    
    return {
        'original_header': header,
        'is_standard': True,
        'instrument_id': parts[0],
        'flowcell_lane': parts[1],
        'tile': parts[2],
        'x_coord': parts[3],
        'y_coord': parts[4],
        'read_number': parts[5],
        'filtered_flag': parts[6]
    }


'''
Full header example:
VH00582:1:AAATJF3HV:1:1101:31183:40756 1:N:0:GATCAG:
Components typically include:

Instrument ID (VH00582)
Flowcell lane (1)
Flowcell tile (AAATJF3HV)
X coordinate (1)
Y coordinate (1101)
Read number (31183)
Filtered flag (40756)
Paired-end information (1:N:0:GATCAG)

1 indicates first or second read in pair
N likely means not filtered
0 might be a control flag
GATCAG could be an index/barcode sequence
'''