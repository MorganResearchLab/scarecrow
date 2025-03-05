#!/usr/bin/env python3
"""
@author: David Wragg
"""

import gzip
import logging
import os
import psutil
import random
import resource
import string
import sys
from functools import lru_cache
from scarecrow.logger import log_errors


class FastqProcessingError(Exception):
    """Custom exception for FASTQ processing errors."""

    pass


def generate_random_string(n: int = 8):
    # Define the character set to use: letters, digits, or customize as needed
    characters = (
        string.ascii_letters + string.digits
    )  # Uppercase, lowercase, and digits
    random_string = "".join(
        random.choices(characters, k=n)
    )  # Generate a random string of length n
    return random_string


def get_memory_usage() -> float:
    """
    Get current memory usage of the process.

    Returns:
        float: Memory usage in MB
    """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def get_process_memory_usage():
    # Get the current process
    main_process = psutil.Process(os.getpid())

    # Get all child processes
    child_processes = main_process.children(recursive=True)

    # Sum memory usage (RSS) of the main process and all child processes
    total_rss = main_process.memory_info().rss  # Main process RSS
    for child in child_processes:
        try:
            total_rss += child.memory_info().rss  # Add child process RSS
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            # Handle cases where the child process no longer exists or access is denied
            pass

    # Convert to MB and GB
    total_rss_mb = total_rss / (1024**2)  # Convert to MB
    total_rss_gb = total_rss / (1024**3)  # Convert to GB
    return total_rss_mb, total_rss_gb


def count_fastq_reads(file):
    opener = gzip.open if file.endswith(".gz") else open
    with opener(file, "rt") as f:
        return sum(1 for line in f) // 4


@lru_cache(maxsize=1024)
def reverse_complement(seq):
    """
    Short function to reverse complement a sequence
    """
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(complement[base] for base in reversed(seq))


@log_errors
def read_barcode_file(file):
    """
    Read barcode sequences from a text file.

    Args:
        file_path (str): Path to the barcode file

    Returns:
        List[str]: List of unique barcode sequences
    """
    logger = logging.getLogger("scarecrow")

    try:
        with open(file, "r") as f:
            # Read lines, strip whitespace, remove empty lines
            barcodes = [line.strip() for line in f if line.strip()]

        # Remove duplicates while preserving order
        unique_barcodes = list(dict.fromkeys(barcodes))

        if not unique_barcodes:
            logger.warning(f"No barcodes found in file: {file}")

        return unique_barcodes

    except FileNotFoundError:
        raise FileNotFoundError(f"Barcode whitelist file not found: {file}")

    except Exception as e:
        logger.error(f"Error reading barcode file {file}: {e}")
        sys.exit(1)


@log_errors
def parse_seed_arguments(barcode_args):
    """
    Parse seed arguments from command line.

    Args:
        barcode_args (List[str]): List of barcode arguments in format 'KEY:WHITELIST:FILE'

    Returns:
        Dict[str, List[str]]: Dictionary of barcodes with keys as region identifiers
    """
    logger = logging.getLogger("scarecrow")

    expected_barcodes = {}

    for arg in barcode_args:
        try:
            # Split the argument into key, whitelist and file path
            key, label, file = arg.split(":")

            # Read barcodes from the file
            barcodes = read_barcode_file(file)

            # Store barcodes in dict under str(key:label)
            expected_barcodes[f"{key}:{label}"] = barcodes
            logger.info(
                f"Loaded {len(barcodes)} barcodes for barcode '{key}' from whitelist '{label}' file '{file}'"
            )

        except ValueError:
            logger.error(
                f"Invalid barcode argument format: {arg}. Use 'INDEX:NAME:FILE'"
            )

    return expected_barcodes
