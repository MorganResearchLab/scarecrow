#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 09:37:12 2024

@author: s14dw4
"""

from seqspec import Assay
from seqspec.utils import load_spec, map_read_id_to_regions
from argparse import RawTextHelpFormatter


def parser_extract(parser):
    subparser = parser.add_parser(
        "extract",
        description="""
Extract cDNA sequence from fastq files

Examples:
scarecrow extract spec.yaml -o ~/path/to/output
---
""",
        help="Extract cDNA from fastqs",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument("yaml", help="Sequencing specification yaml file")
    subparser.add_argument(
        "-o",
        metavar="OUT",
        help=("Path to output cDNA fastq files"),
        type=str,
        default=None,
    )
    subparser.add_argument(
        "-v",
        help=("Verbose output"),
        action='store_true'
    )
    return subparser


def validate_extract_args(parser, args):
    run_extract(yaml = args.yaml, outdir = args.o, verbose = args.v)


def run_extract(yaml, outdir, verbose):
    ''' This function works by:
        1) Going through each modality in the specification
        2) Finding all reads for that modality
        3) Using map_read_id_to_regions to get the associated regions
        4) Collecting the length and sequence information for each region
    '''

    # import seqspec.yaml
    spec = load_spec(yaml)
    
    # extract regions from each read
    read_regions = extract_read_regions_info(spec, verbose)

    if outdir:
        # This needs to be edited still
        print("")
    else:
        print("")
    return spec


def extract_read_regions_info(spec: Assay, verbose = True):
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


