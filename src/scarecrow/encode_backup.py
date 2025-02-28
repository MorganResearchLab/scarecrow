#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
"""

import os
import gzip
import logging
import pickle
from ahocorasick import Automaton
from argparse import RawTextHelpFormatter
from itertools import combinations, product
from collections import defaultdict
from typing import List, Dict, Set, Tuple
from scarecrow import __version__
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string, reverse_complement, parse_seed_arguments

class CompactKmerIndex:
    def __init__(self, k):
        self.k = k
        self.barcodes = []  # List to store unique barcodes
        self.kmer_to_barcode_indices = {}  # Key: integer-encoded k-mer, Value: list of indices in self.barcodes

    def add_barcode(self, barcode):
        """Add a barcode to the index."""
        barcode_index = len(self.barcodes)
        self.barcodes.append(barcode)
        
        # Add all k-mers in the barcode to the index
        for i in range(len(barcode) - self.k + 1):
            kmer = barcode[i:i + self.k]
            kmer_int = self._kmer_to_int(kmer)
            if kmer_int not in self.kmer_to_barcode_indices:
                self.kmer_to_barcode_indices[kmer_int] = []
            self.kmer_to_barcode_indices[kmer_int].append(barcode_index)

    def _kmer_to_int(self, kmer):
        """Convert a k-mer to an integer."""
        encoding = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        result = 0
        for base in kmer:
            result = (result << 2) | encoding.get(base, 0)  # Default to 0 for unknown bases
        return result

    def find_matches(self, sequence):
        """Find approximate matches for a sequence using the k-mer index."""
        matches = set()
        for i in range(len(sequence) - self.k + 1):
            kmer = sequence[i:i + self.k]
            kmer_int = self._kmer_to_int(kmer)
            if kmer_int in self.kmer_to_barcode_indices:
                for barcode_index in self.kmer_to_barcode_indices[kmer_int]:
                    matches.add(self.barcodes[barcode_index])
        return matches
    
class BarcodeMatcher:
    def __init__(self, barcode_sequences: Dict[str, Set[str]], pickle_file: str = None, mismatches: int = 0):
        """
        Base class for barcode matching strategies.
        """
        self.barcode_sequences = barcode_sequences
        self.mismatches = mismatches

class BarcodeMatcherAhoCorasick(BarcodeMatcher):
    def __init__(self, barcode_sequences: Dict[str, Set[str]], pickle_file: str = None, kmer_length: int = 4, jitter: int = 0, mismatches: int = 0):
        """
        Initialize Aho-Corasick based barcode matcher.
        """
        super().__init__(barcode_sequences, mismatches)
        # Initialize automata dictionary only once
        if not hasattr(self, 'automata'):
            self.automata = {}

        self.barcode_info = {}
        self.kmer_index = {}  # k-mer index for approximate matching
        self.kmer_length = kmer_length
        self.jitter = jitter
        self.logger = setup_worker_logger()

        if not pickle_file is None and os.path.exists(pickle_file):
            self._load_trie(pickle_file, mismatches = mismatches)

        else:
            self._build_trie(barcode_sequences)
            self._build_kmer_index(barcode_sequences)  # Build k-mer index for approximate matching
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
        self.logger.info(f"Finished building trie for whitelist: '{whitelist_key}'")
        
    def _build_kmer_index(self, barcode_sequences: Dict[str, Set[str]]):
        """
        Build a k-mer index for approximate matching.
        """
        self.logger.info(f"Building k-mer index for approximate matching using k-mer size {self.kmer_length}")
        for whitelist_key, sequences in barcode_sequences.items():
            kmer_index = defaultdict(set)
            for seq in sequences:
                for i in range(len(seq) - self.kmer_length + 1):
                    kmer = seq[i:i + self.kmer_length]
                    kmer_index[kmer].add(seq)
            self.kmer_index[whitelist_key] = kmer_index
        self.logger.info("Finished building k-mer index")

    def _load_trie(self, pickle_file: str, mismatches: int):
        """
        Loads a pre-built Aho-Corasick trie from a pickle file.
        """
        #self.logger.info(f"Before loading, automata keys: {list(self.automata.keys())}")
        self.logger.info(f"Loading trie and k-mer index from '{pickle_file}'")
        try:
            # Load automata and k-mer index
            with gzip.open(pickle_file, "rb") as f:
                loaded_data = pickle.load(f)

            # Extract the automaton object from the loaded data
            automaton = loaded_data.get('automata', {})
            whitelist_key = list(automaton.keys())[0]
            self.automata[whitelist_key] = automaton[whitelist_key]

            # Extract k-mer index and convert to compact format
            kmer_index_data = loaded_data.get('kmer_index', {})
            self.kmer_index = CompactKmerIndex(k=self.kmer_length)
            for kmer_int, barcode_indices in kmer_index_data.items():
                self.kmer_index.kmer_to_barcode_indices[kmer_int] = barcode_indices
            self.kmer_index.barcodes = loaded_data.get('barcodes', [])

            # Retrieve the k-mer length from the k-mer index
            first_kmer = next(iter(self.kmer_index.kmer_to_barcode_indices[whitelist_key].keys()))
            self.kmer_length = len(first_kmer)            

            self.mismatches = mismatches
            self.logger.info(f"Loaded automata for whitelists: {list(self.automata.keys())}")
            self.logger.info(f"Loaded k-mer indices for whitelists: {list(self.kmer_index.kmer_to_barcode_indices.keys())}")
            self.logger.info(f"Retrieved k-mer length: {self.kmer_length}")

        except:
            self.logger.info(f"Error loading automata and k-mer index from {pickle_file}")
            raise ImportError       

    def _save_trie(self, pickle_file: str):
        """
        Save the Aho-Corasick trie and k-mer index to a pickle file.
        """
        self.logger.info(f"Pickling and compressing Aho-Corasick trie and k-mer index to '{pickle_file}'")
        data_to_save = {
            'automata': self.automata,
            'kmer_index': self.kmer_index
        }
        with gzip.open(pickle_file, "wb") as f:
            pickle.dump(data_to_save, f, protocol=pickle.HIGHEST_PROTOCOL)
        self.logger.info("Pickling complete")

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

    def _find_approximate_matches(self, sequence: str, whitelist_key: str, k: int = 4, max_mismatches: int = 1) -> List[str]:
        """
        Find approximate matches using the k-mer index.
        k = k-mer size
        """
        if whitelist_key not in self.kmer_index.kmer_to_barcode_indices:
            self.logger.warning(f"No k-mer index found for whitelist: {whitelist_key}")
            return []

        # Find candidate barcodes using the k-mer index
        candidates = self.kmer_index.find_matches(sequence)

        # Filter candidates by Hamming distance
        matches = []
        for candidate in candidates:
            dist = self._hamming_distance(sequence, candidate)
            if dist <= max_mismatches:
                matches.append(candidate)
        return matches


    def _hamming_distance(self, s1: str, s2: str) -> int:
        """
        Calculate the Hamming distance between two sequences.
        """
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))

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
        exact_match = False

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
                exact_match = True
            #if matches:
                #self.logger.info(f"{matches}")

            # If no exact matches found, search with the k-mer index
            if not exact_match and self.mismatches > 0:
                approximate_matches = self._find_approximate_matches(seq, whitelist_key, k=self.kmer_length, max_mismatches=self.mismatches)
                #if approximate_matches:
                    #self.logger.info(f"Approximate matches (m={self.mismatches} for {seq}): {approximate_matches}")
                for match in approximate_matches:
                    match_start = start_pos
                    match_dist = abs(match_start - original_start)
                    match_mismatches = self._hamming_distance(seq, match)
                    matches.append({
                        'barcode': match,
                        'whitelist': whitelist_key,
                        'orientation': orientation,
                        'start': match_start,
                        'end': match_start + len(match),
                        'mismatches': match_mismatches,
                        'distance': match_dist
                    })
        
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
                help='Barcode whitelist files in format <barcode_name>:<whitelist_name>:<whitelist_file>\n\t(e.g. BC1:v1:barcodes1.txt)',
        type=str,
        default=None,
    )
    out_format = subparser.add_mutually_exclusive_group(required = True)
    out_format.add_argument(
        "-t", "--trie",
        action='store_true',
        help='Encode whitelist as an Aho-Corasick trie [true]'
    )
    subparser.add_argument(
        "-k", "--kmer_length",
        metavar="<int>",
        help=("K-mer length for building k-mer index for approximate matching [4]"),
        type=int,
        default=4,
    )
    subparser.add_argument(
        "-f", "--force_overwrite",
        action='store_true',
        help='Force overwrite of existing trie file if it exists [false]'
    )
    return subparser

def validate_encode_args(parser, args):
    """ 
    Validate arguments 
    """
    # Global logger setup
    logfile = '{}_{}.{}'.format('./scarecrow_sam2fastq', generate_random_string(), 'log')
    logger = setup_logger(logfile)
    logger.info(f"scarecrow version {__version__}")
    logger.info(f"logfile: '{logfile}'")
     
    run_encode(barcodes = args.barcodes,
               out_trie = args.trie,
               force_overwrite = args.force_overwrite,
               kmer_length = args.kmer_length)
    
@log_errors
def run_encode(barcodes: str = None,
               out_trie: bool = True,
               force_overwrite: bool = False,
               kmer_length: int = 4) -> None:
    """
    Function to encode whitelist in a format suitable for efficient barcode matching
    """
    logger = logging.getLogger('scarecrow')    
    
    # Parse barcode whitelists
    key, label, file_path = barcodes.split(':')
    logger.info(f"Parsing '{key}' '{label}' in {file_path}")
    if os.path.exists(file_path):
        expected_barcodes = parse_seed_arguments([barcodes])

        # Generate Aho-Corasick Trie
        pickle_file = f'{file_path}.k{kmer_length}.trie.gz'
        if out_trie:
            if force_overwrite and os.path.exists(pickle_file):
                logger.info(f"Removing existing file '{pickle_file}'")
                os.remove(pickle_file)

            matcher = BarcodeMatcherAhoCorasick(
                    barcode_sequences={k: set(v) for k, v in expected_barcodes.items()},
                    pickle_file = pickle_file,
                    kmer_length = kmer_length
                )
    else:
        logger.info(f"{file_path} not found")
    
    logger.info("Finished!")

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