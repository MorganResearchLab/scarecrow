#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
"""

import ast
import os
import gzip
import logging
import pickle
from ahocorasick import Automaton
from argparse import RawTextHelpFormatter
from itertools import combinations, product
from typing import List, Dict, Set
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string, reverse_complement, parse_seed_arguments

class BarcodeMatcher:
    def __init__(self, barcode_sequences: Dict[str, Set[str]], pickle_file: str = None, mismatches: int = 0):
        """
        Base class for barcode matching strategies.
        """
        self.barcode_sequences = barcode_sequences
        self.mismatches = mismatches

class BarcodeMatcherAhoCorasick(BarcodeMatcher):
    def __init__(self, barcode_sequences: Dict[str, Set[str]], pickle_file: str = None, jitter: int = 0, mismatches: int = 0):
        """
        Initialize Aho-Corasick based barcode matcher.
        """
        super().__init__(barcode_sequences, mismatches)
        # Initialize automata dictionary only once
        if not hasattr(self, 'automata'):
            self.automata = {}

        self.barcode_info = {}
        self.jitter = jitter
        self.logger = setup_worker_logger()

        if not pickle_file is None and os.path.exists(pickle_file):
            self._load_trie(pickle_file, mismatches = mismatches)

        else:
            self._build_trie(barcode_sequences)                        
            self._save_trie(pickle_file)

    def _build_trie(self, barcode_sequences: Dict[str, Set[str]]):
        """
        Constructs the Aho-Corasick trie from barcode sequences and their single-N mismatch variants.
        """
        total_barcodes = sum(len(v) for v in barcode_sequences.values())
        self.logger.info(f"Building Aho-Corasick trie for {total_barcodes} barcodes")

        for whitelist_key, sequences in barcode_sequences.items():
            automaton = Automaton()
            for i, seq in enumerate(sequences, start=1):
                if i % 1000000 == 0 or i == total_barcodes:
                    self.logger.info(f"Progress {100*(i/total_barcodes):.2f}%")

                # Add the original sequence to the automaton
                automaton.add_word(seq, (whitelist_key, seq, 0))
                
        automaton.make_automaton()
        self.automata[str(whitelist_key)] = automaton
        self.logger.info(f"Finished building trie for whitelist: {whitelist_key}")
        
    def _load_trie(self, pickle_file: str, mismatches: int):
        """
        Loads a pre-built Aho-Corasick trie from a pickle file.
        """
        self.logger.info(f"Before loading, automata keys: {list(self.automata.keys())}")
        self.logger.info(f"Loading Aho-Corasick trie from '{pickle_file}'")
        with gzip.open(pickle_file, "rb") as f:
            loaded_data = pickle.load(f)
        
        # Extract the automaton object from the loaded data
        if isinstance(loaded_data, dict):
            # Assuming the dictionary contains a single automaton under a key
            whitelist_key = list(loaded_data.keys())[0]
            automaton = loaded_data[whitelist_key]
        else:
            # If the loaded data is not a dictionary, assume it is the automaton itself
            whitelist_key = "default_key"  # Use a meaningful key if possible
            automaton = loaded_data

        self.automata[whitelist_key] = automaton
        self.mismatches = mismatches
                
        self.logger.info(f"After loading, automata keys: {list(self.automata.keys())}\n")


    def _save_trie(self, pickle_file: str):
        """
        Saves the current Aho-Corasick trie to a pickle file.
        """
        self.logger.info(f"Pickling and compressing Aho-Corasick trie to '{pickle_file}'")
        with gzip.open(pickle_file, "wb") as f:
            pickle.dump(self.automata, f, protocol=pickle.HIGHEST_PROTOCOL)
        self.logger.info(f"Pickling complete")

    def _generate_query_variants(self, sequence: str, mismatches: int) -> List[str]:
        """
        Generates all possible variants of the query sequence with up to 'n' mismatches,
        treating 'N' as a wildcard.
        """
        bases = ['A', 'C', 'G', 'T']
        variants = set()
        sequence_length = len(sequence)

        # Iterate through mismatch counts (1 to n mismatches)
        for mismatch_count in range(1, mismatches + 1):
            # Get all combinations of mismatch_count positions in the sequence
            for mismatch_positions in combinations(range(sequence_length), mismatch_count):
                # For each position, create all possible base replacements
                for replacements in product(bases, repeat=mismatch_count):
                    # Create a mutable list of characters from the sequence
                    variant_list = list(sequence)

                    # Replace specified positions with replacement bases
                    for pos, base in zip(mismatch_positions, replacements):
                        if base != sequence[pos]:  # Avoid generating the original sequence
                            variant_list[pos] = base

                    # Convert list back to a string and add to the set
                    variants.add(''.join(variant_list))

        return list(variants)

    def _search_with_n_mismatches(self, sequence: str, mismatches: int, orientation: str, automaton: Automaton) -> List[Dict]:
        """
        Searches for matches with up to 'n' mismatches by generating query variants
        and checking them against the trie.
        """
        matches = []
        query_variants = self._generate_query_variants(sequence, mismatches)
        #self.logger.info(f"query variants: {query_variants}")
        for variant in query_variants:
            # Use the automaton's iter method to search for matches
            for end_index, (whitelist_key, variant_seq, n) in automaton.iter(variant):                
                start_index = end_index - len(variant_seq) + 1
                #self.logger.info(f"-> variant_seq: {variant_seq} start_index: {start_index} end_index: {end_index} n: {n} mismatches: {mismatches}")
                matches.append({
                    'barcode': variant_seq,
                    'whitelist': whitelist_key,
                    'orientation': orientation,
                    'start': start_index + 1,
                    'end': end_index + 1,
                    'mismatches': mismatches  # Since we're generating n-mismatch variants
                })

        return matches

    def find_matches(self, sequence: List[tuple[str,int]], whitelist_key: str, orientation: str, original_start: int) -> List[Dict]:
        """
        Find all matching barcodes in a sequence using Aho-Corasick.
        If no exact match is found, searches for matches with up to 'n' mismatches.
        """
        matches = []

        if whitelist_key not in self.automata:
            self.logger.error(f"No automaton found for whitelist: {whitelist_key}")
            return matches
        
        automaton = self.automata[whitelist_key]

        # Iterate over possible sequences and their start positions
        for seq, start_pos in sequence:
            #self.logger.info(f"-> seq {seq} start_pos: {start_pos}")

            # Reverse complement sequence if required
            if orientation == 'reverse':
                seq = reverse_complement(seq)

            # Check for exact matches in the automaton
            for end_index, (wl_key, original_seq, mismatches) in automaton.iter(seq):
                match_start = start_pos + (end_index - len(original_seq)) 
                match_dist = abs((match_start + 1) - original_start)
                #self.logger.info(f"match_start: {match_start} start_pos: {original_start} match_dist: {match_dist}")
                matches.append({
                    'barcode': original_seq,
                    'whitelist': wl_key,
                    'orientation': orientation,
                    'start': match_start + 1, 
                    'end': match_start + len(original_seq), # this was - 1 but removed to get correct end for seed
                    'mismatches': mismatches,
                    'distance': match_dist  
                })
            #if matches:
                #self.logger.info(f"{matches}")

        # If no exact matches found, search with up to 'n' mismatches
        if not matches and self.mismatches > 0:
            for n in range(1, self.mismatches + 1):
                #self.logger.info(f"Checking {n} mismatches")
                for seq, start_pos in sequence:
                    n_mismatch_matches = self._search_with_n_mismatches(seq, n, orientation, automaton)
                    for match in n_mismatch_matches:
                        match_start = start_pos + match['start'] - 1 # Adjust start position
                        match_dist = abs(match_start - original_start)  # Distance from expected start 
                        matches.append({
                            'barcode': match['barcode'],
                            'whitelist': match['whitelist'],
                            'orientation': orientation,
                            'start': match_start, 
                            'end': match_start + len(match['barcode']), # this was - 1 but removed to get correct end for seed
                            'mismatches': match['mismatches'],
                            'distance': match_dist
                        })
                        #self.logger.info(f"\n\t{seq}\tstart_pos: {start_pos} oroginal_start: {original_start} match_start: {match_start} match_dist: {match_dist} mismatches: {match['mismatches']}")

        return sorted(matches, key=lambda x: x['start'])

def parser_encode(parser):
    subparser = parser.add_parser(
        "encode",
        description="""
Generate Aho-Corasick Trie.

Example:

scarecrow encode --barcodes whitelist.txt --trie
---
""",
        epilog="If the out file exists then this will be loaded rather than re-generated.",
        help="Generate Aho-Corasick trie and pickle to compressed file",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument(
        "-b", "--barcodes",
        metavar="<file>",
        help=("Path to barcode text file"),
        type=str,
        default=None,
    )
    out_format = subparser.add_mutually_exclusive_group(required = True)
    out_format.add_argument(
        "-t", "--trie",
        action='store_true',
        help='Encode whitelist as an Aho-Corasick trie [true]'
    )
    return subparser

def validate_encode_args(parser, args):
    run_encode(barcodes = args.barcodes,
             out_trie = args.trie,)
    
@log_errors
def run_encode(barcodes: str = None,
               out_trie: bool = True) -> None:
    """
    Function to encode whitelist in a format suitable for efficient barcode matching
    """
    # Setup logging
    logfile = f'./scarecrow_encode_{generate_random_string()}.log'
    logger = setup_logger(logfile)
    logger.info(f"logfile: '{logfile}'")
    
    # Parse barcode whitelists
    key, label, file_path = barcodes.split(':')
    logger.info(f"Parsing '{key}' '{label}' in {file_path}")
    if os.path.exists(file_path):
        expected_barcodes = parse_seed_arguments([barcodes])

        # Generate Aho-Corasick Trie
        if out_trie:
            matcher = BarcodeMatcherAhoCorasick(
                    barcode_sequences={k: set(v) for k, v in expected_barcodes.items()},
                    pickle_file =  f'{file_path}.trie.gz'
                )
    else:
        logger.info(f"{file_path} not found")

def setup_worker_logger(log_file: str = None):
    """Configure logger for worker processes with file output"""
    logger = logging.getLogger('scarecrow')
    if not logger.handlers:  # Only add handlers if none exist
        # Create formatters
        formatter = logging.Formatter('%(asctime)s - %(processName)s - %(name)s - %(levelname)s - %(message)s')
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
        
        # File handler if log_file provided
        if log_file:
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
        
        logger.setLevel(logging.INFO)
    return logger   