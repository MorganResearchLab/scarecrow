#!/usr/bin/env python3
"""
@author: David Wragg
"""

import resource
import gzip
import logging
from scarecrow.logger import log_errors
from typing import List, Dict, Set
import string, random

class FastqProcessingError(Exception):
    """Custom exception for FASTQ processing errors."""
    pass

def generate_random_string(n: int = 8):
    # Define the character set to use: letters, digits, or customize as needed
    characters = string.ascii_letters + string.digits  # Uppercase, lowercase, and digits
    random_string = ''.join(random.choices(characters, k=n))  # Generate a random string of length n
    return random_string

def get_memory_usage() -> float:
    """
    Get current memory usage of the process.
    
    Returns:
        float: Memory usage in MB
    """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024

def count_fastq_reads(filename):
    opener = gzip.open if filename.endswith('.gz') else open
    with opener(filename, 'rt') as f:
        return sum(1 for line in f) // 4
    
def reverse_complement(seq):
    """
    Short function to reverse complement a sequence
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq))

@log_errors
def region_indices(spec: Assay, fastqs):
    """
    Identify library elements contained in sequencing reads
    """
    logger = logging.getLogger('scarecrow')

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

