#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 09:37:12 2024

@author: David Wragg
"""

from seqspec import Assay
from seqspec.utils import load_spec, map_read_id_to_regions
from seqspec.seqspec_print import run_seqspec_print
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
        1) Iterate through each modality in the specification
        2) Finding all reads for that modality
        3) Use map_read_id_to_regions to get the associated regions
        4) Collect the length and sequence information for each region

        process_fastq_files(fastq, paired, verbose, read_regions)
            => match_fastq_to_regions
            => process_fastq_in_batches
            => process_batch
    """
    
    # Run seqspec print to get format
    print(f"\033[32m\nseqspec print {yaml}:\033[0m\n")
    run_seqspec_print(yaml, fmt="library-ascii", o = None)

    # import seqspec.yaml
    spec = load_spec(yaml)

    # extract regions from each read
    regions = extract_regions(spec)
    print(f"\033[32m\nregion coordinates:\033[0m\n")
    for region in regions:
        print(f"{region.region_id} : {region.start}-{region.stop}")

    read_regions = extract_read_regions_info(spec)
    
    # process fastq files
    process_fastq_files(fastqs, paired, verbose, read_regions)

    if outdir:
        # This needs to be edited still
        print("")
    else:
        print("")
    return spec



def extract_regions(spec: Assay, verbose: bool = True):
    """
    Extract subregion start stop positions
    """
    from seqspec.Region import project_regions_to_coordinates

    cuts = []
    # Iterate through each modality in the assay
    for modality in spec.list_modalities():
        # Get all reads for this modality
        libspec = spec.get_libspec(modality)
        leaves = libspec.get_leaves()
        cuts = project_regions_to_coordinates(leaves)
    return cuts



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



def process_fastq_in_batches(r1_file: str, r2_file: str, regions_info, verbose: bool = False, batch_size: int = 1000):
    """
    Process paired FASTQ reads in batches
    """
    batch = []
    read_count = 0
    
    print(f"\033[32mProcessing: \033[0m{r1_file} \033[32mand \033[0m{r2_file}")
    
    for (header1, seq1, qual1), (header2, seq2, qual2) in read_fastq_pairs(r1_file, r2_file):
        batch.append(((header1, seq1, qual1), (header2, seq2, qual2)))
        read_count += 1
        
        if len(batch) >= batch_size:
            # Process the batch
            process_batch(batch, regions_info, verbose)
            batch = []
            
            # Print progress
            print(f"Processed {read_count:,} read pairs...", end='\r')
    
    # Process any remaining reads
    if batch:
        process_batch(batch, regions_info)
    
    print(f"\nCompleted processing {read_count:,} read pairs")


def read_fastq_pairs(r1_file: str, r2_file: str) -> Iterator[Tuple[Tuple[str, str, str], Tuple[str, str, str]]]:
    """
    Read paired FASTQ files using pysam
    
    Args:
        r1_file (str): Path to R1 FASTQ file
        r2_file (str): Path to R2 FASTQ file
        
    Yields:
        ((header1, seq1, qual1), (header2, seq2, qual1))
    """
    import pysam
    
    with pysam.FastxFile(r1_file) as r1, pysam.FastxFile(r2_file) as r2:
        for read1, read2 in zip(r1, r2):
            # Yield headers, sequences and qualities for each pair
            yield ((read1.name, read1.sequence, read1.quality), 
                   (read2.name, read2.sequence, read2.quality))



def extract_subregion(regions_info, seq, qual):
    """ 
    Extract subregions from region
    """
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
    
    return extracted_regions


def process_batch(batch, regions_info, verbose):
    """
    Process a batch of read pairs
    """
    for (header1, seq1, qual1), (header2, seq2, qual2) in batch:
        # Parse headers
        parsed_header1 = parse_fastq_header(header1)
        parsed_header2 = parse_fastq_header(header2)
        
        if not parsed_header1['is_standard'] or not parsed_header2['is_standard']:
            print(f"Warning: Non-standard header in pair: {header1} or {header2}")
            continue
            
        # Process read pairs...
        if len(regions_info) > 1:            
            if verbose:
                print(f"\n\033[32mHeader 1: \033[0m{header1}")
                subregions = extract_subregion(regions_info[0], seq1, qual1)
                if len(subregions) > 0:
                    for region_id, region_data in subregions.items():
                        print(f"{region_id}: {region_data['sequence']}")

                print(f"\n\033[32mHeader 2: \033[0m{header2}")
                subregions = extract_subregion(regions_info[1], seq2, qual2)
                if len(subregions) > 0:
                    for region_id, region_data in subregions.items():
                        print(f"{region_id}: {region_data['sequence']}")



def process_fastq_files(fastqs: list, paired: bool, verbose: bool, read_regions):
    """
    Process FASTQ reads using their region specifications
    """
    try:
        if paired:
            # Assume fastqs are ordered R1, R2
            for i in range(0, len(fastqs), 2):
                r1_file = fastqs[i]
                r2_file = fastqs[i + 1]                
                regions_info_r1 = match_fastq_to_regions(read_regions, Path(r1_file).name)
                regions_info_r2 = match_fastq_to_regions(read_regions, Path(r2_file).name)
                process_fastq_in_batches(r1_file, r2_file, [regions_info_r1, regions_info_r2], verbose, batch_size=1000)
        else:
            for file in fastqs:
                regions_info = match_fastq_to_regions(read_regions, Path(file).name)
                print(f"\nProcessing {file}")
                print(f"Modality: {regions_info['modality']}")
                process_fastq_in_batches(file, None, regions_info, verbose, batch_size=1000)
    
    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"Unexpected error processing files: {e}")



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