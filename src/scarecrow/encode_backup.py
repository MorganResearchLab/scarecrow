#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
"""

import os
import gzip
import logging
import pickle
import numpy as np
from ahocorasick import Automaton
from argparse import RawTextHelpFormatter
from itertools import combinations, product
from collections import defaultdict
from pympler import asizeof
from typing import List, Dict, Set, Tuple
from scarecrow import __version__
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string, reverse_complement, parse_seed_arguments

class CompactKmerIndex:
    def __init__(self, k: int):
        self.k = k  # k-mer size
        self.barcodes = []  # List to store unique barcodes
        self.kmer_to_barcode_indices = {}  # Temporary dictionary for fast lookups during build
        self.kmer_array = None  # Will store k-mers as a numpy array after finalization
        self.barcode_indices = None  # Will store barcode indices as a numpy array after finalization

    def add_barcode(self, barcode: str):
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

    def finalize(self):
        """
        Convert the temporary dictionary to numpy arrays for memory efficiency.
        This should be called after all barcodes have been added.
        """
        # Convert k-mer dictionary to numpy arrays
        self.kmer_array = np.array(list(self.kmer_to_barcode_indices.keys()), dtype=np.uint16)
        self.barcode_indices = np.array(list(self.kmer_to_barcode_indices.values()), dtype=object)

        # Free the temporary dictionary
        del self.kmer_to_barcode_indices

    def _kmer_to_int(self, kmer: str) -> int:
        """Convert a k-mer to an integer using 2 bits per nucleotide."""
        encoding = {'A': 0b00, 'C': 0b01, 'G': 0b10, 'T': 0b11}
        result = 0
        for base in kmer:
            result = (result << 2) | encoding.get(base, 0)  # Default to 0 for unknown bases
        return result

    def find_matches(self, sequence: str, max_mismatches: int = 1) -> Set[str]:
        """
        Find approximate matches for a sequence using the k-mer index.
        """
        matches = set()
        sequence_length = len(sequence)

        # Encode the entire sequence into a bit representation
        sequence_int = self._kmer_to_int(sequence)

        # Iterate through all k-mers in the sequence
        for i in range(sequence_length - self.k + 1):
            # Extract the k-mer from the sequence
            kmer_int = (sequence_int >> (2 * (sequence_length - self.k - i))) & ((1 << (2 * self.k)) - 1)

            # Find candidate barcodes using the k-mer index
            if kmer_int in self.kmer_array:
                idx = np.where(self.kmer_array == kmer_int)[0][0]
                for barcode_index in self.barcode_indices[idx]:
                    barcode = self.barcodes[barcode_index]
                    # Calculate Hamming distance between the sequence and the barcode
                    if self._hamming_distance(sequence, barcode) <= max_mismatches:
                        matches.add(barcode)
        return matches

    def _hamming_distance(self, s1: str, s2: str) -> int:
        """
        Calculate the Hamming distance between two sequences.
        """
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))
    
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
        self.kmer_index = {} 
        self.kmer_length = kmer_length
        self.jitter = jitter
        self.logger = setup_worker_logger()

        if not pickle_file is None and os.path.exists(pickle_file):
            self._load_trie(pickle_file, mismatches = mismatches)

        else:
            self._build_trie(barcode_sequences)
            self._build_kmer_index(barcode_sequences)
            self._save_trie(pickle_file)
        
        self.report_memory_usage()

    def report_memory_usage(self):
        """
        Report the memory usage of the trie and k-mer index using pympler.asizeof.
        This method is called from `self` and has access to all class attributes.
        """
        trie_memory = 0
        kmer_index_memory = 0

        # Calculate memory usage of the trie
        for whitelist_key, trie in self.automata.items():
            trie_memory += asizeof.asizeof(trie)

        # Calculate memory usage of the k-mer index
        if hasattr(self, 'kmer_index') and self.kmer_index:
            for whitelist_key, kmer_index in self.kmer_index.items():
                # Memory usage of kmer_array
                kmer_index_memory += asizeof.asizeof(kmer_index.kmer_array)
                # Memory usage of barcode_indices
                kmer_index_memory += asizeof.asizeof(kmer_index.barcode_indices)
                # Memory usage of barcodes
                kmer_index_memory += asizeof.asizeof(kmer_index.barcodes)

        # Print memory usage
        self.logger.info(f"Trie memory usage: {trie_memory / 1024 / 1024:.2f} MB")
        self.logger.info(f"K-mer index memory usage: {kmer_index_memory / 1024 / 1024:.2f} MB")

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

                # Encode the sequence into a compact string
                encoded_seq = encode_dna_to_compact_string(seq)
                # Add the encoded sequence to the automaton
                automaton.add_word(encoded_seq, (whitelist_key, seq, 0))
                
        automaton.make_automaton()
        self.automata[str(whitelist_key)] = automaton
        self.logger.info(f"Finished building trie for whitelist: '{whitelist_key}'")
        
    def _build_kmer_index(self, barcode_sequences: Dict[str, Set[str]]):
        """
        Build a k-mer index for approximate matching.
        """
        self.logger.info(f"Building k-mer index for approximate matching using k-mer size {self.kmer_length}")
        for whitelist_key, sequences in barcode_sequences.items():
            kmer_index = CompactKmerIndex(k=self.kmer_length)
            for seq in sequences:
                kmer_index.add_barcode(seq)
            kmer_index.finalize()  # Convert to numpy arrays after all barcodes are added
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

            # Load k-mer index
            self.kmer_length = loaded_data.get('kmer_length', {})
            kmer_index_data = loaded_data.get('kmer_index', {})
            for whitelist_key, kmer_index_data in kmer_index_data.items():
                kmer_index = CompactKmerIndex(k=self.kmer_length)
                kmer_index.barcodes = kmer_index_data.get('barcodes', [])
                kmer_index.kmer_array = np.array(kmer_index_data.get('kmer_array', []), dtype=np.uint16)
                kmer_index.barcode_indices = kmer_index_data.get('barcode_indices', [])  # Already a list of lists
                self.kmer_index[whitelist_key] = kmer_index

            self.mismatches = mismatches
            self.logger.info(f"Loaded Aho-Corasick trie for whitelists: {list(self.automata.keys())}")
            self.logger.info(f"Loaded k-mer indices for whitelists: {list(self.kmer_index.keys())}")
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
            'kmer_index': {
                whitelist_key: {
                    'barcodes': kmer_index.barcodes,
                    'kmer_array': kmer_index.kmer_array.tolist(),  # Convert numpy array to list
                    'barcode_indices': kmer_index.barcode_indices  # Already a list of lists
                }
                for whitelist_key, kmer_index in self.kmer_index.items()
            },
            'kmer_length': self.kmer_length
        }
        with gzip.open(pickle_file, "wb") as f:
            pickle.dump(data_to_save, f, protocol=pickle.HIGHEST_PROTOCOL)
        self.logger.info("Pickling complete")

    def output_kmer_index_details(self, num_keys=5, num_values=2):
        """
        Output details about the k-mer index for the first whitelist.
        """
        if not self.kmer_index:
            self.logger.info("K-mer index is empty.")
            return

        # Get the first whitelist and its k-mer index
        first_whitelist_key = next(iter(self.kmer_index))
        kmer_index = self.kmer_index[first_whitelist_key]

        self.logger.info(f"\nK-mer Index Details for Whitelist: '{first_whitelist_key}'")
        self.logger.info(f"Total number of k-mers: {len(kmer_index.kmer_array)}")
        self.logger.info(f"Total number of whitelists: {len(self.kmer_index)}")

        # Output the first few k-mers and their values
        for i in range(min(num_keys, len(kmer_index.kmer_array))):
            kmer_int = kmer_index.kmer_array[i]
            barcode_indices = kmer_index.barcode_indices[i]

            # Decode the k-mer integer back to a string
            kmer = int_to_kmer(kmer_int, kmer_index.k)

            # Get the first few barcodes
            barcodes = [kmer_index.barcodes[idx] for idx in barcode_indices[:num_values]]
            total_values = len(barcode_indices)

            # Log the k-mer, first few barcodes, and total count of barcodes
            self.logger.info(f"  K-mer: {kmer}")
            self.logger.info(f"    First {num_values} barcodes: {barcodes}")
            self.logger.info(f"    Total barcodes: {total_values}")

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
        if whitelist_key not in self.kmer_index:
            self.logger.warning(f"No k-mer index found for whitelist: {whitelist_key}")
            return []

        # Get the CompactKmerIndex instance for the given whitelist_key
        kmer_index = self.kmer_index[whitelist_key]

        # Find candidate barcodes using the k-mer index
        candidates = kmer_index.find_matches(sequence)

        # Filter candidates by Hamming distance
        matches = []
        for candidate in candidates:
            dist = hamming_distance(sequence, candidate)
            if dist <= max_mismatches:
                matches.append(candidate)
        return matches

    def find_matches(self, sequence: List[tuple[str,int]], whitelist_key: str, orientation: str, original_start: int) -> List[Dict]:
        """
        Find all matching barcodes in a sequence using Aho-Corasick.
        If no exact match is found, searches for approximate matches up to 'n' mismatches.
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

            # Encode the sequence into a compact string
            encoded_seq = encode_dna_to_compact_string(seq)

            # Check for exact matches in the automaton
            for end_index, (wl_key, original_seq, mismatches) in automaton.iter(encoded_seq):
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
                    match_mismatches = hamming_distance(seq, match)
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
    logger.info(f"Parsing '{key}' '{label}' in '{file_path}'")
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
    
    if matcher:
        # Output k-mer index
        matcher.output_kmer_index_details(num_keys=1, num_values=3)


    logger.info("Finished!")

def int_to_kmer(kmer_int, k):
    """
    Convert an integer-encoded k-mer back to a string.
    """
    encoding = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    kmer = []
    for _ in range(k):
        kmer.append(encoding[kmer_int & 0b11])  # Extract the last 2 bits
        kmer_int >>= 2  # Shift right by 2 bits
    return ''.join(reversed(kmer))

def hamming_distance(s1, s2):
    """
    Calculate the Hamming distance between two sequences.
    """
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def encode_dna_to_compact_string(sequence):
    """
    Encodes a DNA sequence into a compact string representation using 2 bits per nucleotide.
    """
    nucleotide_to_char = {'A': '0', 'C': '1', 'G': '2', 'T': '3'}
    encoded = []
    for nucleotide in sequence:
        encoded.append(nucleotide_to_char.get(nucleotide, '0'))  # Default to '0' for unknown bases
    return ''.join(encoded)

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



