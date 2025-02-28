#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Wragg
"""

import os
import gzip
import logging
import pickle
import random
from argparse import RawTextHelpFormatter
from itertools import combinations, product
from collections import defaultdict
from typing import List, Dict, Set, Tuple
from scarecrow import __version__
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string, reverse_complement, parse_seed_arguments

class CompactTrieNode:
    def __init__(self):
        self.children = {}  # Key: byte, Value: index of child node
        self.is_end = False  # Marks the end of a barcode

class CompactTrie:
    def __init__(self):
        self.nodes = [CompactTrieNode()]  # Root node at index 0

    def insert(self, sequence):
        """Insert a sequence into the trie."""
        node_index = 0  # Start at the root
        for byte in sequence.encode('ascii'):  # Convert to bytes
            node = self.nodes[node_index]
            if byte not in node.children:
                # Add a new node
                new_node_index = len(self.nodes)
                self.nodes.append(CompactTrieNode())
                node.children[byte] = new_node_index
            node_index = node.children[byte]
        self.nodes[node_index].is_end = True  # Mark end of sequence

    def find_matches(self, sequence):
        """Find all matches in a sequence."""
        matches = []
        for i in range(len(sequence)):
            node_index = 0  # Start at the root
            for j in range(i, len(sequence)):
                byte = sequence[j].encode('ascii')[0]  # Convert to byte
                if byte not in self.nodes[node_index].children:
                    break  # No match
                node_index = self.nodes[node_index].children[byte]
                if self.nodes[node_index].is_end:
                    matches.append((i-1, j, sequence[i:j+1]))  # Match found
        return matches

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

class BarcodeMatcherTrie(BarcodeMatcher):
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
        Constructs the custom trie from barcode sequences
        """
        total_barcodes = sum(len(v) for v in barcode_sequences.values())
        self.logger.info(f"Building custom trie for {total_barcodes} barcodes")

        for whitelist_key, sequences in barcode_sequences.items():
            # Initialize a custom trie for this whitelist
            custom_trie = CompactTrie()

            # Insert all sequences into the custom trie
            for i, seq in enumerate(sequences, start=1):
                if i % 1000000 == 0 or i == total_barcodes:
                    self.logger.info(f"Progress {100*(i/total_barcodes):.2f}%")

                # Add the original sequence to the custom trie
                custom_trie.insert(seq)

            # Store the custom trie in the automata dictionary
            self.automata[str(whitelist_key)] = custom_trie
            self.logger.info(f"Finished building custom trie for whitelist: '{whitelist_key}'")
        
    def _build_kmer_index(self, barcode_sequences: Dict[str, Set[str]]):
        """
        Build a k-mer index for approximate matching using the CompactKmerIndex class.
        """
        self.logger.info(f"Building k-mer index for approximate matching using k-mer size {self.kmer_length}")

        # Initialize the kmer_index as a dictionary to store k-mer indices per whitelist
        self.kmer_index = {}

        # Add barcodes to the k-mer index for each whitelist
        for whitelist_key, sequences in barcode_sequences.items():
            # Initialize a CompactKmerIndex for this whitelist
            kmer_index = CompactKmerIndex(k = self.kmer_length)

            # Add barcodes to the k-mer index
            for seq in sequences:
                kmer_index.add_barcode(seq)

            # Store the k-mer index for this whitelist
            self.kmer_index[whitelist_key] = kmer_index

        self.logger.info("Finished building k-mer index")

    def _load_trie(self, pickle_file: str, mismatches: int):
        """
        Loads a pre-built custom trie and k-mer index from a pickle file.
        """
        self.logger.info(f"Loading custom trie and k-mer index from '{pickle_file}'")
        try:
            # Load data from the pickle file
            with gzip.open(pickle_file, "rb") as f:
                loaded_data = pickle.load(f)

            # Load the custom trie
            self.automata = loaded_data.get('automata', {})
            if not self.automata:
                raise ValueError("No automata found in the pickle file.")

            # Load the k-mer index in a compact format
            kmer_index_data = loaded_data.get('kmer_index', {})
            if not kmer_index_data:
                raise ValueError("No k-mer index found in the pickle file.")

            # Initialize the kmer_index as a dictionary to store k-mer indices per whitelist
            self.kmer_index = {}

            # Load the k-mer index for each whitelist
            for whitelist_key, kmer_index_data_for_whitelist in kmer_index_data.items():
                # Initialize the CompactKmerIndex with the correct k-mer length
                self.kmer_length = kmer_index_data_for_whitelist.get('kmer_length')
                kmer_index = CompactKmerIndex(k = self.kmer_length)

                # Populate the k-mer index
                for kmer_int, barcode_indices in kmer_index_data_for_whitelist.get('kmer_to_barcode_indices', {}).items():
                    kmer_index.kmer_to_barcode_indices[kmer_int] = barcode_indices
                kmer_index.barcodes = kmer_index_data_for_whitelist.get('barcodes', [])

                # Store the k-mer index for this whitelist
                self.kmer_index[whitelist_key] = kmer_index

            self.mismatches = mismatches
            self.logger.info(f"Extracted custom trie(s) for whitelist(s): {list(self.automata.keys())}")
            self.logger.info(f"Extracted k-mer indices for whitelist(s): {list(self.kmer_index.keys())}")
            self.logger.info(f"Retrieved k-mer length: {self.kmer_length}")

        except Exception as e:
            self.logger.error(f"Error loading custom trie and k-mer index from {pickle_file}: {str(e)}")
            raise ImportError

    def _save_trie(self, pickle_file: str):
        """
        Save the custom trie and k-mer index to a pickle file.
        """
        self.logger.info(f"Pickling and compressing custom trie and k-mer index to '{pickle_file}'")
        try:
            # Prepare data to save
            kmer_index_data = {}
            for whitelist_key, kmer_index in self.kmer_index.items():
                kmer_index_data[whitelist_key] = {
                    'kmer_to_barcode_indices': kmer_index.kmer_to_barcode_indices,
                    'barcodes': kmer_index.barcodes,
                    'kmer_length': kmer_index.k  # Save the k-mer length
                }

            data_to_save = {
                'automata': self.automata,
                'kmer_index': kmer_index_data
            }

            # Save data to a compressed pickle file
            with gzip.open(pickle_file, "wb") as f:
                pickle.dump(data_to_save, f, protocol=pickle.HIGHEST_PROTOCOL)
            
            self.logger.info("Pickling complete")

        except Exception as e:
            self.logger.error(f"Error saving custom trie and k-mer index to {pickle_file}: {str(e)}")
            raise

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

    def _find_approximate_matches(self, sequence: str, whitelist_key: str, k: int = 4, max_mismatches: int = 1) -> List[str]:
        """
        Find approximate matches using the k-mer index.
        k = k-mer size
        """
        if whitelist_key not in self.kmer_index:
            self.logger.warning(f"No k-mer index found for whitelist: {whitelist_key}")
            return []

        # Get the k-mer index for this whitelist
        kmer_index = self.kmer_index[whitelist_key]

        # Find candidate barcodes using the k-mer index
        candidates = kmer_index.find_matches(sequence)

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
        
        # Use the custom trie for matching
        trie = self.automata[whitelist_key]
        exact_match = False

        # Iterate over possible sequences and their start positions
        for seq, start_pos in sequence:
            if orientation == 'reverse':
                seq = reverse_complement(seq)

            # Find matches in the custom trie
            trie_matches = trie.find_matches(seq)
            for match_start, match_end, matched_seq in trie_matches:
                match_dist = abs((start_pos + match_start) - original_start)
                matches.append({
                    'barcode': matched_seq,
                    'whitelist': whitelist_key,
                    'orientation': orientation,
                    'start': start_pos + match_start + 1,
                    'end': start_pos + match_end,
                    'mismatches': 0,  # Exact match
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

        # Generate custom trie and kmer index
        pickle_file = f'{file_path}.k{kmer_length}.trie.gz'
        if out_trie:
            if force_overwrite and os.path.exists(pickle_file):
                logger.info(f"Removing existing file '{pickle_file}'")
                os.remove(pickle_file)

            matcher = BarcodeMatcherTrie(
                    barcode_sequences={k: set(v) for k, v in expected_barcodes.items()},
                    pickle_file = pickle_file,
                    kmer_length = kmer_length
                )
    else:
        logger.info(f"{file_path} not found")

    if matcher:
        # Print a manageable subset of the structure
        output_trie_and_kmer_index(matcher, max_depth=matcher.kmer_length, sample_size=5, max_barcodes=2)
    
    logger.info("Finished!")


def output_trie(trie, max_depth=3, max_children=5,  num_base_nodes=3):
    """
    Print the structure of the trie up to a specified depth, summarizing large subtrees.
    """
    if not trie:
        print("Trie is empty.")
        return
    
    _output_trie_node(trie, 0, "", max_depth, depth=0)


def _output_trie_node(trie, node_index, prefix, max_depth, depth):
    """
    Recursively print a node in the trie, following a single path to the last level.
    """
    node = trie.nodes[node_index]
    byte_value = chr(list(node.children.keys())[0]) if node.children else "End"
    print(prefix + "└── " + f"{byte_value} (depth {depth})")

    # Stop if we've reached the penultimate level
    if depth >= max_depth - 1:
        # Print all nodes at the penultimate level
        children = list(node.children.items())
        if children:
            for byte, child_index in children:
                child_node = trie.nodes[child_index]
                child_byte_value = chr(list(child_node.children.keys())[0]) if child_node.children else "End"
                print(prefix + "    └── " + f"{child_byte_value} (depth {depth + 1})")
        return

    # Update the prefix for child nodes
    new_prefix = prefix + "    "

    # Print children, following a single path
    children = list(node.children.items())
    if children:
        # Follow the first child
        first_child_byte, first_child_index = children[0]
        _output_trie_node(trie, first_child_index, new_prefix, max_depth, depth + 1)

        # Print the number of unprinted children
        if len(children) > 1:
            print(new_prefix + f"├── ... ({len(children) - 1} more children)")

def output_kmer_index(kmer_indices, sample_size=10, max_barcodes=3):
    """
    Print a random sample of k-mers and their associated barcodes for each whitelist,
    limiting the number of barcodes displayed.
    """
    if not kmer_indices:
        print("K-mer index is empty.")
        return

    for whitelist_key, kmer_index in kmer_indices.items():
        print(f"K-mer Index for Whitelist: {whitelist_key}")
        print(f"Random sample of {sample_size} k-mers (max {max_barcodes} barcodes per k-mer):")

        if not kmer_index.kmer_to_barcode_indices:
            print("  No k-mers found for this whitelist.")
            continue

        # Get a random sample of k-mers
        kmer_samples = random.sample(
            list(kmer_index.kmer_to_barcode_indices.items()),
            min(sample_size, len(kmer_index.kmer_to_barcode_indices))
        )

        for kmer_int, barcode_indices in kmer_samples:
            # Decode the k-mer integer back to a string
            kmer = int_to_kmer(kmer_int, kmer_index.k)
            barcodes = [kmer_index.barcodes[idx] for idx in barcode_indices[:max_barcodes]]
            if len(barcode_indices) > max_barcodes:
                barcodes.append(f"... ({len(barcode_indices) - max_barcodes} more)")
            print(f"  K-mer: {kmer} -> Barcodes: {barcodes}")
        print()

def int_to_kmer(kmer_int, k):
    """
    Convert an integer-encoded k-mer back to a string.
    """
    encoding = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    kmer = []
    for _ in range(k):
        kmer.append(encoding[kmer_int & 0b11])
        kmer_int >>= 2
    return ''.join(reversed(kmer))

def output_trie_and_kmer_index(matcher, max_depth=3, max_children=5, num_base_nodes=3, sample_size=10, max_barcodes=3):
    """
    Print a manageable subset of the trie and k-mer index structure.
    """
    print("")
    for whitelist_key, trie in matcher.automata.items():
        print(f"Whitelist: {whitelist_key}")
        print("Trie Structure (single path to last level):")
        output_trie(trie, max_depth=max_depth)
        print()

    # Print k-mer indices for each whitelist
    if hasattr(matcher, 'kmer_index') and matcher.kmer_index:
        output_kmer_index(matcher.kmer_index, sample_size=sample_size, max_barcodes=max_barcodes)
    else:
        print("No k-mer index found.")

    print("")

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

