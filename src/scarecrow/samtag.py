#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Low-memory SAM/FASTQ tag processing for large genomic files
"""

import os
import sys
import logging
import pysam
import argparse
from typing import Dict, Generator, Optional

def create_argument_parser():
    """Create argument parser for the script."""
    parser = argparse.ArgumentParser(
        description="Low-memory SAM/FASTQ tag processing",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-f", "--fastq", 
        required=True, 
        help="Path to input FASTQ file"
    )
    parser.add_argument(
        "-s", "--sam", 
        required=True, 
        help="Path to input SAM/BAM file"
    )
    parser.add_argument(
        "-o", "--output", 
        help="Path to output SAM/BAM file (default: input filename with .tagged suffix)"
    )
    parser.add_argument(
        "-m", "--max-read-buffer", 
        type=int, 
        default=10000, 
        help="Maximum number of reads to buffer (default: 10000)"
    )
    return parser

def setup_logging():
    """Configure logging."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger(__name__)

def stream_fastq_tags(fastq_path: str) -> Generator[Dict[str, Optional[str]], None, None]:
    """
    Stream tags from FASTQ file with minimal memory usage.
    
    Yields dictionaries of tags for each read, preserving memory efficiency.
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Extracting tags from FASTQ: {fastq_path}")
    
    tag_names = ["CR", "CY", "CB", "XP", "XM", "UR", "UY"]
    
    with open(fastq_path, 'r') as fq:
        while True:
            # Read header line
            header = fq.readline().strip()
            if not header:
                break
            
            # Skip sequence, '+', and quality lines
            fq.readline()
            fq.readline()
            fq.readline()
            
            # Process only FASTQ headers starting with '@'
            if header.startswith('@'):
                # Extract read name (first part of header)
                read_name = header.split()[0][1:]
                
                # Extract tags
                tags = {}
                for item in header.split():
                    if '=' in item:
                        key, value = item.split('=')
                        tags[key] = value
                
                # Prepare tag dictionary with specified tag names
                yield {
                    'read_name': read_name,
                    **{tag: tags.get(tag) for tag in tag_names}
                }

def process_sam_with_tags(
    input_path: str, 
    output_path: str, 
    tag_generator: Generator[Dict[str, Optional[str]], None, None],
    max_read_buffer: int = 10000
):
    """
    Process SAM/BAM file with extremely low memory usage.
    
    Streams reads, looks up tags, and writes to output with minimal memory overhead.
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Processing SAM/BAM: {input_path}")
    
    # Prepare tag lookup
    tag_cache: Dict[str, Dict[str, Optional[str]]] = {}
    tag_names = ["CR", "CY", "CB", "XP", "XM", "UR", "UY"]
    
    # Open input and output files
    with pysam.AlignmentFile(input_path, 'rb') as infile, \
         pysam.AlignmentFile(output_path, 'wb', header=infile.header) as outfile:
        
        # Stream tags and build a small, rotating cache
        for tag_entry in tag_generator:
            read_name = tag_entry.pop('read_name')
            tag_cache[read_name] = tag_entry
            
            # Limit cache size
            if len(tag_cache) > max_read_buffer:
                oldest_read = min(tag_cache.keys())
                del tag_cache[oldest_read]
        
        # Reset file pointer
        infile.reset()
        
        # Process reads
        for read in infile:
            # Check if read name exists in current tag cache
            if read.query_name in tag_cache:
                tags = tag_cache[read.query_name]
                
                # Add tags if they exist
                for tag_name in tag_names:
                    tag_value = tags.get(tag_name)
                    if tag_value:
                        read.set_tag(tag_name, tag_value, value_type='Z')
            
            # Write read to output
            outfile.write(read)
        
        logger.info(f"Completed processing: {output_path}")

def main():
    """Main execution function."""
    # Parse arguments
    parser = create_argument_parser()
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging()
    
    try:
        # Determine output path
        if not args.output:
            base, ext = os.path.splitext(args.sam)
            args.output = f"{base}.tagged{ext}"
        
        # Stream tags from FASTQ
        tag_generator = stream_fastq_tags(args.fastq)
        
        # Process SAM/BAM with tags
        process_sam_with_tags(
            args.sam, 
            args.output, 
            tag_generator,
            max_read_buffer=args.max_read_buffer
        )
        
        logger.info("Processing completed successfully.")
    
    except Exception as e:
        logger.error(f"Error processing files: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()