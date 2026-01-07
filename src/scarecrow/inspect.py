# -*- coding: utf-8 -*-
"""
#!/usr/bin/env python3
@author: David Wragg
"""

import os
import gzip
import logging
import pickle
import sys
from pathlib import Path
import numpy as np
from argparse import ArgumentTypeError, RawTextHelpFormatter
from scarecrow import __version__
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string

def parser_inspect(parser):
    subparser = parser.add_parser(
        "inspect",
        description="""
Inspect a pickle file and return a summary.

Example:

scarecrow inspect --in BC1.seed.pkl.gz
---
""",
        help="Inspects a pickle file and returns a summary of the data structure",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument(
        "-i",
        "--in",
        dest="infile",
        metavar="<file>",
        help=("Pickle file generated via scarecrow seed or scarecrow encode using seed or kmer indexing"),
        type=str,
        required=True,
        default=None,
    )
    return subparser


def validate_inspect_args(parser, args) -> None:
    """
    Validate arguments
    """
    # Global logger setup
    logfile = "{}_{}.{}".format(
        "./scarecrow_inspect", generate_random_string(), "log"
    )
    logger = setup_logger(logfile)
    logger.info(f"scarecrow version {__version__}")
    logger.info(f"logfile: '{logfile}'")

    run_inspect(input_file=args.infile)


@log_errors
def run_inspect(input_file: str = None) -> None:
    """
    Main function to inspect pickle file
    """
    logger = logging.getLogger("scarecrow")

    try:
        # Validate input file
        if not isinstance(input_file, str):
            raise TypeError("Input file path must be a string")

        input_path = Path(input_file)
        if not input_path.exists():
            raise FileNotFoundError(f"Input file does not exist: {input_file}")
        if not os.access(input_file, os.R_OK):
            raise PermissionError(f"Input file is not readable: {input_file}")

        # Validate file extension
        if input_path.suffixes[-2:] == ['.pkl', '.gz']:
            inspect_pkl(input_file)

        else:
            raise ValueError("Input file extension not .pkl.gz")

    except (TypeError, ValueError, FileNotFoundError, PermissionError) as e:
        logger.error(f"Validation error: {str(e)}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error during validation: {str(e)}")
        raise


def preview_list(arr, name="Array", limit=3):
    """
    Returns a string preview of the first `limit` elements of a list/array.
    """
    arr_len = len(arr)
    if arr_len == 0:
        return "[]"
    preview = arr[:limit]
    suffix = "..." if arr_len > limit else ""
    return f"{preview}{suffix} (total {arr_len})"

def preview_array(arr, name="Array", limit=3):
    """
    Returns a string preview of the first `limit` elements of a numpy array.
    """
    if arr is None or len(arr) == 0:
        return "[]"
    arr_len = len(arr)
    preview = arr[:limit].tolist()
    suffix = "..." if arr_len > limit else ""
    return f"{preview}{suffix} (total {arr_len})"

@log_errors
def inspect_pkl(pickle_file: str = None):
    """
    Read pickle file
    """
    logger = logging.getLogger("scarecrow")
    with gzip.open(pickle_file, "rb") as f:
        data = pickle.load(f)

    logger.info(f"Top-level keys: {list(data.keys())}")
    logger.info(f"Index type: {data.get('index_type')}")
    logger.info(f"Mismatches: {data.get('index_mismatches')}")
    logger.info(f"Barcode length: {data.get('barcode_length')}")
    logger.info(f"kmer length: {data.get('kmer_length')}")
    logger.info(f"Pickle version: {data.get('version')}")

    automata = data.get("automata", {})
    logger.info(f"Whitelists in automata: {list(automata.keys())}")

    kmer_index = data.get("kmer_index", {})
    for whitelist_key, idx_data in kmer_index.items():
        logger.info(f"\nWhitelist: {whitelist_key}")
        logger.info(f"  Number of barcodes: {len(idx_data.get('barcodes', []))}")
        logger.info(f"  Index type stored: {idx_data.get('index_type')}")
        logger.info(f"  Has seed_index?: {'Yes' if idx_data.get('seed_index') else 'No'}")

        # Preview automaton/trie
        trie_preview = preview_list(list(automata.get(whitelist_key, {}).keys())
                                    if isinstance(automata.get(whitelist_key), dict)
                                    else list(automata.get(whitelist_key, [])), limit=3)
        logger.info(f"  Automaton (first 3 keys/values): {trie_preview}")

        index_type = idx_data.get("index_type")

        # Index-specific preview
        if index_type == "kmer":
            kmer_array = idx_data.get("kmer_array", [])
            kmer_preview = preview_array(np.array(kmer_array), limit=3)
            logger.info(
                f"  k-mer array (first 3 values): {kmer_preview}"
            )

        elif index_type == "seed":
            seed_index = idx_data.get("seed_index") or {}

            seed_lengths = sorted(seed_index.keys())
            preview_lengths = preview_list(seed_lengths, limit=3)
            logger.info(
                f"  Seed index lengths (first 3): {preview_lengths}"
            )

            if seed_lengths:
                L = seed_lengths[0]
                seeds = list(seed_index[L].keys())
                seed_preview = preview_list(seeds, limit=3)
                logger.info(
                    f"  Example seeds (length {L}): {seed_preview}"
                )

        else:
            logger.info("  No index data available for preview")
