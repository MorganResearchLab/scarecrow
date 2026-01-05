# -*- coding: utf-8 -*-
"""
#!/usr/bin/env python3
@author: David Wragg

Revision uses seed indexing as a replacement for k-mer indexing for the approximate matching.

1.	Exact matching with the trie (Aho–Corasick)
2.	Approximate matching using multi-seed indexing (mismatch-aware)
	•	Divide each barcode into d+1 disjoint seeds (d = user-defined max mismatches).
	•	Index each seed separately (length-bucketed).
	•	Lookup candidate barcodes by seed matches.
	•	Union candidate sets and verify with bit-packed Hamming.
	•	Guarantees that at least one seed hits if the barcode is within d mismatches.
3.	Bit-packed verification
	•	Encode barcodes as integers (uint32/64) depending on length.
	•	Query sequence → XOR with candidate barcode → popcount → Hamming distance.
	•	Extremely fast and works with any d.


Original k-mer approach worked by chopping whitelist barcode into overlapping k-mers. Look-up table maps kmer to list of barcodes that contain it.
When a read barcode fails exact matching, all k-mers are extracted for the read barcode, union all barcode lists for those kmers, check Hamming distance.
Problem: single kmer can occur in thousands of barcodes, problem is amplified when also considering jitter and mismatches.

Seed index splits whitelist barcode into a few seeds. Number and placement of seeds is chosen so that if two barcodes differ by < m mismatches, at least one
seed must match exactly. Index maps seed to barcodes containing seed at that position. When a read barcode fails exact matching, extract same seeds from the read,
lookup only barcodes sharing at least one seed, check those candidates with Hamming distance. This is faster because each seed lookup returns far fewer barcodes.
This approach is more 'aligner-like' rather than brute-force.

'To improve the performance and correctness of approximate barcode matching, we replaced the original k-mer–based lookup with a seed index–based approach.
The previous implementation indexed all overlapping k-mers from each barcode, which could generate large candidate sets when common k-mers were shared across many barcodes,
leading to unnecessary comparisons. The new seed index partitions each barcode into a small number of carefully chosen seeds derived from the barcode length and the maximum
number of allowed mismatches. This guarantees that any barcode within the mismatch threshold will share at least one exact seed with the query sequence, while substantially
reducing the number of candidate barcodes that must be evaluated. Final barcode assignment is still determined using full Hamming distance, ensuring identical matching
behavior while providing more predictable and scalable performance.'
"""

import os
import gzip
import logging
import pickle
import zlib
import numpy as np
from ahocorasick import Automaton
from argparse import RawTextHelpFormatter
from itertools import combinations, product
from collections import namedtuple, defaultdict
from pympler import asizeof
from typing import List, Dict, Set
from scarecrow import __version__
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import (
    generate_random_string,
    reverse_complement,
    parse_seed_arguments,
)

# Define a namedtuple for match results
MatchResult = namedtuple(
    "MatchResult",
    ["barcode", "whitelist", "orientation", "start", "end", "mismatches", "distance"],
)


class SeedIndex:
    """
    Mismatch-aware seed index for variable-length barcodes.
    Each barcode is split into max_mismatches + 1 disjoint seeds (pigeonhole principle),
    ensuring at least one exact seed match if barcode is within max_mismatches.
    """

    def __init__(self, max_mismatches: int):
        self.mismatches = max_mismatches
        self.barcodes = []  # original barcodes
        # seed_index[length][seed_id][seed_string] -> set of barcode indices
        self.seed_index = {}

    def add_barcode(self, barcode: str):
        idx = len(self.barcodes)
        self.barcodes.append(barcode)
        L = len(barcode)
        d = self.mismatches
        seed_count = d + 1
        seed_len = L // seed_count

        self.seed_index.setdefault(L, {})
        for seed_id in range(seed_count):
            start = seed_id * seed_len
            end = start + seed_len if seed_id < seed_count - 1 else L
            seed = barcode[start:end]
            self.seed_index[L].setdefault(seed_id, {}).setdefault(seed, set()).add(idx)

    def build_from_barcodes(self, barcodes, logger=None):
        total_barcodes = len(barcodes)
        report_interval = max(1_000_000, total_barcodes // 100)  # log every 1M or 1% of total
        for i, bc in enumerate(barcodes):
            self.add_barcode(bc)
            # Report progress every 10% or for very small sets every barcode
            if logger and (i % report_interval == 0 or i == total_barcodes):
                percent = 100 * i / total_barcodes
                logger.info(f"Progress: {percent:.2f}%")


    def query(self, sequence: str):
        """
        Return all barcodes within max_mismatches of the input sequence.
        """
        L = len(sequence)
        if L not in self.seed_index:
            return []

        d = self.mismatches
        seed_count = d + 1
        seed_len = L // seed_count
        candidates = set()

        for seed_id in range(seed_count):
            start = seed_id * seed_len
            end = start + seed_len if seed_id < seed_count - 1 else L
            seed = sequence[start:end]
            candidates.update(self.seed_index[L][seed_id].get(seed, set()))

        # Bit-packed Hamming verification
        query_bits = _encode_dna_to_bits(sequence)
        matches = []
        for idx in candidates:
            barcode_bits = _encode_dna_to_bits(self.barcodes[idx])
            if _hamming_distance_bits(query_bits, barcode_bits) <= d:
                matches.append(self.barcodes[idx])
        return matches



class CompactKmerIndex:
    def __init__(self, k: int | None = None):
        self.k = k
        self.barcodes = []  # List to store unique barcodes
        self.kmer_to_barcode_indices = {}  # Temporary dictionary for fast lookups during build
        self.kmer_array = None  # Will store k-mers as a numpy array after finalization
        self.barcode_indices = (
            None  # Will store barcode indices as a numpy array after finalization
        )

    def add_barcode(self, barcode: str):
        """Add a barcode to the index."""
        barcode_index = len(self.barcodes)
        self.barcodes.append(barcode)

        # Add all k-mers in the barcode to the index
        for i in range(len(barcode) - self.k + 1):
            kmer = barcode[i : i + self.k]
            kmer_int = self._kmer_to_int(kmer)

            # Add the barcode index to the k-mer's list of indices
            if kmer_int not in self.kmer_to_barcode_indices:
                self.kmer_to_barcode_indices[kmer_int] = []
            self.kmer_to_barcode_indices[kmer_int].append(barcode_index)

    def finalize(self):
        """
        Convert the temporary dictionary to numpy arrays for memory efficiency.
        This should be called after all barcodes have been added.
        """
        # Convert k-mer dictionary to numpy arrays
        if self.k > 8:
            self.kmer_array = np.array(
                list(self.kmer_to_barcode_indices.keys()), dtype=np.uint32
            )
        else:
            self.kmer_array = np.array(
                list(self.kmer_to_barcode_indices.keys()), dtype=np.uint16
            )

        self.barcode_indices = np.array(
            list(self.kmer_to_barcode_indices.values()), dtype=object
        )

        # Free the temporary dictionary
        del self.kmer_to_barcode_indices

    def _kmer_to_int(self, kmer: str) -> int:
        """Convert a k-mer to an integer using 2 bits per nucleotide."""
        encoding = {"A": 0b00, "C": 0b01, "G": 0b10, "T": 0b11}
        result = 0
        for base in kmer:
            result = (result << 2) | encoding.get(
                base, 0
            )  # Default to 0 for unknown bases
        return result

    def query(self, sequence: str, max_mismatches: int = None) -> Set[str]:
        """
        Return all barcodes within max_mismatches of the input sequence.
        Uses k-mer seed lookup followed by bit-packed Hamming verification.
        """
        if max_mismatches is None:
            max_mismatches = 1

        matches = set()

        if self.kmer_array is None or self.kmer_array.size == 0 or len(self.barcodes) == 0:
            return matches

        # Determine effective k if not set
        k = self.k
        n_seeds = max_mismatches + 1
        seq_len = len(sequence)
        if k is None or k <= 0:
            k = seq_len // n_seeds

        # Encode full sequence into bits
        seq_bits = _encode_dna_to_bits(sequence)

        seen_candidates = set()
        for i in range(seq_len - k + 1):
            kmer = sequence[i : i + k]
            kmer_int = self._kmer_to_int(kmer)
            # Find matching barcodes for this k-mer
            idx_matches = np.where(self.kmer_array == kmer_int)[0]
            for idx in idx_matches:
                for barcode_index in self.barcode_indices[idx]:
                    if barcode_index in seen_candidates:
                        continue
                    candidate = self.barcodes[barcode_index]
                    barcode_bits = _encode_dna_to_bits(candidate)
                    if _hamming_distance_bits(seq_bits, barcode_bits) <= max_mismatches:
                        matches.add(candidate)
                    seen_candidates.add(barcode_index)

        return matches


    def _hamming_distance(self, s1: str, s2: str) -> int:
        """
        Calculate the Hamming distance between two sequences.
        """
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))


class BarcodeMatcher:
    def __init__(
        self,
        barcode_sequences: Dict[str, Set[str]],
        pickle_file: str = None,
        mismatches: int = 0,
    ):
        """
        Base class for barcode matching strategies.
        """
        self.barcode_sequences = barcode_sequences
        self.mismatches = mismatches


class BarcodeMatcherAhoCorasick(BarcodeMatcher):
    __slots__ = [
        "barcode_sequences",
        "mismatches",
        "automata",
        "index",
        "kmer_index",
        "kmer_length",
        "jitter",
        "logger",
    ]

    def __init__(
        self,
        barcode_sequences: Dict[str, Set[str]],
        pickle_file: str = None,
        index_type: str = None,
        kmer_length: int | None = None,
        jitter: int = 0,
        mismatches: int = 0,
    ):
        """
        Initialize Aho-Corasick based barcode matcher.
        """
        super().__init__(barcode_sequences, mismatches)
        # Initialize automata dictionary only once
        if not hasattr(self, "automata"):
            self.automata = {}

        self.index_type = index_type
        self.kmer_index = {}
        self.kmer_length = kmer_length
        self.jitter = jitter
        self.mismatches = mismatches
        self.logger = setup_worker_logger()

        if pickle_file is not None and os.path.exists(pickle_file):
            self._load_trie(pickle_file, mismatches = mismatches)

        else:
            self._build_trie(barcode_sequences)
            self._build_index(barcode_sequences, index_type = self.index_type)
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
        if hasattr(self, "kmer_index") and self.kmer_index:
            for whitelist_key, kmer_index in self.kmer_index.items():
                # Memory usage of the temporary dictionary (if not finalized)
                if hasattr(kmer_index, "kmer_to_barcode_indices"):
                    kmer_index_memory += asizeof.asizeof(
                        kmer_index.kmer_to_barcode_indices
                    )
                # Memory usage of kmer_array (if finalized)
                if hasattr(kmer_index, "kmer_array"):
                    kmer_index_memory += asizeof.asizeof(kmer_index.kmer_array)
                # Memory usage of barcode_indices (if finalized)
                if hasattr(kmer_index, "barcode_indices"):
                    kmer_index_memory += asizeof.asizeof(kmer_index.barcode_indices)
                # Memory usage of barcodes
                if hasattr(kmer_index, "barcodes"):
                    kmer_index_memory += asizeof.asizeof(kmer_index.barcodes)

        # Print memory usage
        self.logger.info(f"trie memory usage: {trie_memory / 1024 / 1024:.2f} MB")
        self.logger.info(f"{self.index_type} index memory usage: {kmer_index_memory / 1024 / 1024:.2f} MB")

    def _build_trie(self, barcode_sequences: Dict[str, Set[str]]):
        """
        Constructs the Aho-Corasick trie from barcode sequences
        """
        total_barcodes = sum(len(v) for v in barcode_sequences.values())
        self.logger.info(f"Building Aho-Corasick trie for {total_barcodes} barcodes")

        for whitelist_key, sequences in barcode_sequences.items():
            automaton = Automaton()
            for i, seq in enumerate(sequences, start=1):
                if i % 1000000 == 0 or i == total_barcodes:
                    self.logger.info(f"Progress: {100 * (i / total_barcodes):.2f}%")

                # Encode the sequence into a compact string
                encoded_seq = encode_dna_to_compact_string(seq)
                # Add the encoded sequence to the automaton
                automaton.add_word(encoded_seq, (whitelist_key, seq, 0))

        automaton.make_automaton()
        self.automata[str(whitelist_key)] = automaton
        self.logger.info(f"Finished building trie for whitelist: '{whitelist_key}'")


    def _build_seed_index(self, barcode_sequences: Dict[str, Set[str]]):
        """Build SeedIndex for approximate matching."""
        self.logger.info(
            f"Building seed index for approximate matching with max_mismatches={self.mismatches}"
        )
        for whitelist_key, sequences in barcode_sequences.items():
            self.logger.info(f"Building seed index for whitelist: '{whitelist_key}'")
            seed_index = SeedIndex(max_mismatches=self.mismatches)
            seed_index.build_from_barcodes(list(sequences), logger=self.logger)
            self.kmer_index[whitelist_key] = seed_index
        self.logger.info("Finished building seed index")


    def _build_compact_kmer_index(self, barcode_sequences: Dict[str, Set[str]]):
        """
        Build CompactKmerIndex for approximate matching using k-mers.
        k is derived from barcode length and mismatch budget:
            k = floor(L / (mismatches + 1))
        """
        d = self.mismatches
        self.logger.info(
            f"Building k-mer index for approximate matching (mismatches={d})"
        )

        for whitelist_key, sequences in barcode_sequences.items():
            self.logger.info(f"Building kmer index for whitelist: '{whitelist_key}'")
            if not sequences:
                raise ValueError(f"Empty whitelist: {whitelist_key}")

            # Ensure uniform barcode length
            lengths = {len(seq) for seq in sequences}
            if len(lengths) != 1:
                raise ValueError(
                    f"All barcodes in whitelist '{whitelist_key}' must have the same length "
                    f"for k-mer indexing. Found lengths: {sorted(lengths)}"
                )
            L = lengths.pop()
            k = max(1, L // (d + 1))
            self.logger.info(
                f"Derived k-mer length for whitelist '{whitelist_key}': "
                f"L={L}, mismatches={d} → k={k}"
            )
            kmer_index = CompactKmerIndex(k=k)

            total_seqs = len(sequences)
            report_interval = max(1_000_000, total_seqs // 100)

            for i, seq in enumerate(sequences, start=1):
                kmer_index.add_barcode(seq)
                if self.logger and (i % report_interval == 0 or i == total_seqs):
                    percent = 100 * i / total_seqs
                    self.logger.info(f"Progress: {percent:.2f}%")

            kmer_index.finalize()
            self.kmer_index[whitelist_key] = kmer_index

        self.logger.info("Finished building k-mer index")


    def _build_index(self, barcode_sequences: Dict[str, Set[str]], index_type: str = "seed"):
        """
        Dispatch function to build either seed index or k-mer index.
        """
        if index_type == "seed":
            self._build_seed_index(barcode_sequences)
        elif index_type == "kmer":
            self._build_compact_kmer_index(barcode_sequences)
        else:
            raise ValueError(f"Unknown index_type: {index_type}")
        self.index_type = index_type  # track the built index type


    def _load_trie(self, pickle_file: str, mismatches: int):
        """
        Loads a pre-built Aho-Corasick trie and index from a pickle file.
        """
        self.logger.info(f"Importing pickle from '{pickle_file}'")
        try:
            # Load automata and k-mer index data
            with gzip.open(pickle_file, "rb") as f:
                loaded_data = pickle.load(f)

            # Set mismatch budget
            self.mismatches = mismatches

            # Load trie
            automaton = loaded_data.get("automata", {})
            for whitelist_key, trie in automaton.items():
                self.automata[whitelist_key] = trie
            self.logger.info(f"Loaded trie for whitelists: {list(self.automata.keys())}")

            # Load k-mer index (legacy or new SeedIndex)
            kmer_index_data = loaded_data.get("kmer_index", {})
            self.kmer_length = loaded_data.get("kmer_length", None)
            version = loaded_data.get("version", "unknown")
            self.logger.info(f"Pickle version: {version}")

            for whitelist_key, ki_data in kmer_index_data.items():
                index_type = ki_data.get("index_type", "seed")
                self.index_type = index_type

                if index_type == "seed":
                    idx_obj = SeedIndex(max_mismatches=self.mismatches)
                    idx_obj.build_from_barcodes(ki_data.get("barcodes", []), logger=self.logger)
                elif index_type == "kmer":
                    idx_obj = CompactKmerIndex(k=self.kmer_length)
                    idx_obj.barcodes = ki_data.get("barcodes", [])
                    arr = ki_data.get("kmer_array", [])
                    if self.kmer_length > 8:
                        idx_obj.kmer_array = np.array(arr, dtype=np.uint32)
                    else:
                        idx_obj.kmer_array = np.array(arr, dtype=np.uint16)
                    idx_obj.barcode_indices = ki_data.get("barcode_indices", [])
                else:
                    raise ValueError(f"Unknown index_type in pickle: {index_type}")

                self.kmer_index[whitelist_key] = idx_obj
                self.logger.info(f"Loaded index for whitelists: {list(self.kmer_index.keys())}")

        except (ImportError, FileNotFoundError) as e:
            self.logger.error(f"Error loading trie and index from {pickle_file}: {e}")
            raise ImportError

    def _save_trie(self, pickle_file: str):
        """
        Save the Aho-Corasick trie and k-mer index to a pickle file,
        """
        self.logger.info(
            f"Pickling and compressing Aho-Corasick trie and index to '{pickle_file}'"
        )

        kmer_index_serializable = {}
        for whitelist_key, idx_obj in self.kmer_index.items():

            if isinstance(idx_obj, SeedIndex):
                index_type = "seed"
                barcodes = idx_obj.barcodes
                kmer_array = []
                barcode_indices = []
            elif isinstance(idx_obj, CompactKmerIndex):
                index_type = "kmer"
                barcodes = idx_obj.barcodes
                kmer_array = idx_obj.kmer_array.tolist() if idx_obj.kmer_array is not None else []
                barcode_indices = [list(v) for v in idx_obj.barcode_indices] if idx_obj.barcode_indices is not None else []
            else:
                raise TypeError(f"Unknown index object type: {type(idx_obj)}")

            kmer_index_serializable[whitelist_key] = {
                "index_type": index_type,
                "barcodes": barcodes,
                "kmer_array": kmer_array,
                "barcode_indices": barcode_indices,
            }


        data_to_save = {
            "version": __version__,
            "automata": self.automata,
            "kmer_index": kmer_index_serializable,
            "kmer_length": self.kmer_length,
        }

        with gzip.open(pickle_file, "wb") as f:
            pickle.dump(data_to_save, f, protocol=pickle.HIGHEST_PROTOCOL)

        self.logger.info("Pickling complete")


    def output_kmer_index_details(self, num_keys=5, num_values=2):
        """
        Output details about the seed index for the first whitelist.
        """
        if not self.kmer_index:
            self.logger.info("Seed index is empty.")
            return

        whitelist_key = next(iter(self.kmer_index))
        idx = self.kmer_index[whitelist_key]

        self.logger.info(f"\nIndex details for whitelist: '{whitelist_key}'")


        if isinstance(idx, SeedIndex):
            total_barcodes = len(idx.barcodes)
            lengths = sorted(idx.seed_index.keys())

            total_seeds = 0
            for L in lengths:
                for seed_id, seed_dict in idx.seed_index[L].items():
                    total_seeds += len(seed_dict)

            self.logger.info("Index type: SeedIndex")
            self.logger.info(f"  Total barcodes: {total_barcodes}")
            self.logger.info(f"  Barcode lengths: {lengths}")
            self.logger.info(f"  Total unique seeds: {total_seeds}")

            # Show a few example seeds
            shown = 0
            for L in lengths:
                for seed_id, seed_dict in idx.seed_index[L].items():
                    for seed, barcode_ids in seed_dict.items():
                        barcodes = [idx.barcodes[i] for i in list(barcode_ids)[:num_values]]
                        self.logger.info(f"  Seed (L={L}, id={seed_id}): {seed}")
                        self.logger.info(f"    Example barcodes: {barcodes}")
                        shown += 1
                        if shown >= num_keys:
                            return

        elif isinstance(idx, CompactKmerIndex):
            self.logger.info("Index type: CompactKmerIndex")
            self.logger.info(f"  k-mer length: {idx.k}")
            self.logger.info(f"  Total k-mers: {len(idx.kmer_array)}")
            self.logger.info(f"  Total barcodes: {len(idx.barcodes)}")

            for i in range(min(num_keys, len(idx.kmer_array))):
                kmer_int = idx.kmer_array[i]
                barcode_ids = idx.barcode_indices[i][:num_values]
                barcodes = [idx.barcodes[j] for j in barcode_ids]
                kmer = int_to_kmer(kmer_int, idx.k)

                self.logger.info(f"  k-mer: {kmer}")
                self.logger.info(f"    Example barcodes: {barcodes}")

        else:
            self.logger.warning(f"Unknown index type: {type(idx)}")


    def _find_approximate_matches(
        self, sequence: str, whitelist_key: str, k=None, max_mismatches=None
    ) -> list[str]:
        """
        Find approximate matches using the SeedIndex while ensuring mismatch counts
        are consistent with the per-barcode max_mismatches logic.

        Parameters:
            sequence: the read sequence segment to match
            whitelist_key: which whitelist to query
            k, max_mismatches: ignored for backward compatibility

        Returns:
            List of barcodes from the whitelist matching within max_mismatches
        """
        if whitelist_key not in self.kmer_index:
            self.logger.warning(f"No seed index found for whitelist: {whitelist_key}")
            return []

        seed_index = self.kmer_index[whitelist_key]
        candidates = seed_index.query(sequence)  # get candidate barcodes using seeds
        matched_barcodes = []

        for candidate in candidates:
            # compute full Hamming distance between read segment and candidate barcode
            dist = hamming_distance(sequence, candidate)
            if dist <= self.mismatches:  # only accept if within user-specified mismatch budget
                matched_barcodes.append(candidate)

        return matched_barcodes

    def find_matches(
        self,
        sequence: List[tuple[str, int]],
        whitelist_key: str,
        orientation: str,
        original_start: int,
    ) -> List[MatchResult]:
        """
        Find all matching barcodes in a sequence using Aho-Corasick.
        If no exact match is found, searches for approximate matches using a k-mer index.
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
            if orientation == "reverse":
                seq = reverse_complement(seq)

            # Encode the sequence into a compact string
            encoded_seq = encode_dna_to_compact_string(seq)

            # Check for exact matches in the automaton
            for end_index, (wl_key, original_seq, mismatches) in automaton.iter(
                encoded_seq
            ):
                pat_len = len(original_seq)
                seq_len = len(seq)

                if orientation == "reverse":
                    r_end = end_index
                    r_start = r_end - pat_len + 1
                    f_start = seq_len - r_end - 1
                    f_end   = seq_len - r_start - 1
                    start = start_pos + f_start
                    end   = start_pos + f_end
                else:
                    match_start = start_pos + (end_index - pat_len + 1)
                    start = match_start
                    end   = match_start + pat_len - 1

                match_dist = abs(start - original_start)
                #self.logger.info(f"{seq} - {orientation} {original_seq} @ match_start: {match_start} expected_start: {original_start} match_dist: {match_dist}")
                matches.append(
                    MatchResult(
                        barcode     = original_seq,
                        whitelist   = wl_key,
                        orientation = orientation,
                        start       = start,
                        end         = end,
                        mismatches  = mismatches,
                        distance    = match_dist,
                    )
                )
                exact_match = True
            # if matches:
            # self.logger.info(f"{matches}")

            # If no exact matches found, search with the k-mer index
            if not exact_match and self.mismatches > 0:
                # self.logger.info(f"Check for approx match for {seq}")
                approximate_matches = self._find_approximate_matches(
                    seq,
                    whitelist_key
                )
                # if approximate_matches:
                #self.logger.info(f"Approximate matches (m={self.mismatches} for {seq}): {approximate_matches}")
                for match in approximate_matches:
                    match_start = start_pos
                    match_dist = abs(match_start - original_start)
                    if(match_start < 1):
                        match_dist -= 1
                    match_mismatches = hamming_distance(seq, match)
                    matches.append(
                        MatchResult(
                            barcode=match,
                            whitelist=whitelist_key,
                            orientation=orientation,
                            start=match_start,
                            end=match_start + len(match),
                            mismatches=match_mismatches,
                            distance=match_dist,
                        )
                    )
                    #if seq[:2] == "NN":
                    #    self.logger.info(f"{seq} match:{match} start:{match_start} end:{match_start + len(match)} distance:{match_dist} mismatches:{match_mismatches}")

        return sorted(matches, key=lambda x: x.start)


def _encode_dna_to_bits(seq: str) -> int:
    """
    Encode a DNA sequence into an integer using 2 bits per base.
    A=00, C=01, G=10, T=11; non-ACGT treated as A.
    """
    encoding = {"A": 0b00, "C": 0b01, "G": 0b10, "T": 0b11}
    val = 0
    for base in seq:
        val = (val << 2) | encoding.get(base, 0)
    return val


def _hamming_distance_bits(a: int, b: int) -> int:
    """
    Hamming distance for bit-packed DNA sequences.
    Each mismatch flips exactly 2 bits.
    """
    x = a ^ b
    return bin(x).count("1") // 2

def parser_encode(parser):
    subparser = parser.add_parser(
        "encode",
        description="""
Generate Aho-Corasick Trie.

Example:

scarecrow encode --barcodes whitelist.txt --pickle
---
""",
        epilog="If the out file exists then this will be loaded rather than re-generated.",
        help="Generate Aho-Corasick trie and k-mer index, pickle to compressed file",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument(
        "-b",
        "--barcodes",
        metavar="<file>",
        help="Barcode whitelist files in format <barcode_name>:<whitelist_name>:<whitelist_file>\n\t(e.g. BC1:v1:barcodes1.txt)",
        type=str,
        default=None,
    )
    out_format = subparser.add_mutually_exclusive_group(required=True)
    out_format.add_argument(
        "-p",
        "--pickle",
        action="store_true",
        help="Pickle whitelist as an Aho-Corasick trie and k-mer index [true]",
    )
    subparser.add_argument(
        "-i",
        "--index",
        metavar="<str>",
        help = ("Indexing method, can be either seed or kmer [seed]"),
        type = str,
        default = "seed",
    )
    subparser.add_argument(
        "-m",
        "--mismatches",
        metavar="<int>",
        type=int,
        default=1,
        help="Number of allowed mismatches in barcode [1]",
    )
    subparser.add_argument(
        "-f",
        "--force_overwrite",
        action="store_true",
        help="Force overwrite of existing trie file if it exists [false]",
    )
    return subparser


def validate_encode_args(parser, args):
    """
    Validate arguments
    """
    # Global logger setup
    logfile = "{}_{}.{}".format(
        "./scarecrow_sam2fastq", generate_random_string(), "log"
    )
    logger = setup_logger(logfile)
    logger.info(f"scarecrow version {__version__}")
    logger.info(f"logfile: '{logfile}'")

    run_encode(
        barcodes = args.barcodes,
        pickle = args.pickle,
        force_overwrite = args.force_overwrite,
        index = args.index,
        mismatches = args.mismatches,
    )


@log_errors
def run_encode(
    barcodes: str = None,
    pickle: bool = True,
    force_overwrite: bool = False,
    index: str = "seed",
    mismatches: int = 1,
) -> None:
    """
    Function to encode whitelist in a format suitable for efficient barcode matching
    """
    logger = logging.getLogger("scarecrow")

    # Parse barcode whitelists
    key, label, file_path = barcodes.split(":")
    logger.info(f"Parsing '{key}' '{label}' in '{file_path}'")
    if os.path.exists(file_path):
        expected_barcodes = parse_seed_arguments([barcodes])

        # Generate Aho-Corasick Trie
        pickle_file = f"{file_path}.{index}.pkl.gz"
        if pickle:
            if force_overwrite and os.path.exists(pickle_file):
                logger.info(f"Removing existing file '{pickle_file}'")
                os.remove(pickle_file)
            matcher = BarcodeMatcherAhoCorasick(
                barcode_sequences = {k: set(v) for k, v in expected_barcodes.items()},
                pickle_file = pickle_file,
                index_type = index,
                mismatches = mismatches,
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
    encoding = {0: "A", 1: "C", 2: "G", 3: "T"}
    kmer = []
    for _ in range(k):
        kmer.append(encoding[kmer_int & 0b11])  # Extract the last 2 bits
        kmer_int >>= 2  # Shift right by 2 bits
    return "".join(reversed(kmer))


def hamming_distance(s1, s2):
    """
    Calculate the Hamming distance between two sequences.
    """
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def encode_dna_to_compact_string(sequence):
    """
    Encodes a DNA sequence into a compact string representation using 2 bits per nucleotide.
    """
    nucleotide_to_char = {"A": "0", "C": "1", "G": "2", "T": "3", "N": "4"}
    encoded = []
    for nucleotide in sequence:
        encoded.append(
            nucleotide_to_char.get(nucleotide, "0")
        )  # Default to '0' for unknown bases
    return "".join(encoded)


def setup_worker_logger(log_file: str = None):
    """Configure logger for worker processes with file output"""
    logger = logging.getLogger("scarecrow")
    if not logger.handlers:  # Only add handlers if none exist
        # Create formatters
        formatter = logging.Formatter(
            "%(asctime)s - %(processName)s - %(name)s - %(levelname)s - %(message)s"
        )

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
