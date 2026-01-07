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
        self.barcodes = []

        # build-time only
        self.kmer_to_barcode_indices = defaultdict(list)

        # finalized arrays (NO dtype=object)
        self.kmer_array = None              # uint16 / uint32 [N]
        self.barcode_offsets = None         # int64 [N+1]
        self.barcode_indices = None         # int32 [M]

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
        Convert the temporary dictionary to compact numeric arrays.
        Eliminates dtype=object completely.
        """
        items = sorted(self.kmer_to_barcode_indices.items())
        n_kmers = len(items)
        total_hits = sum(len(v) for _, v in items)

        if self.k > 8:
            self.kmer_array = np.empty(n_kmers, dtype=np.uint32)
        else:
            self.kmer_array = np.empty(n_kmers, dtype=np.uint16)

        self.barcode_offsets = np.empty(n_kmers + 1, dtype=np.int64)
        self.barcode_indices = np.empty(total_hits, dtype=np.int32)

        offset = 0
        for i, (kmer_int, barcode_list) in enumerate(items):
            self.kmer_array[i] = kmer_int
            self.barcode_offsets[i] = offset
            self.barcode_indices[offset : offset + len(barcode_list)] = barcode_list
            offset += len(barcode_list)

        self.barcode_offsets[n_kmers] = offset

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
        if max_mismatches is None:
            max_mismatches = 1

        matches = set()

        if (
            self.kmer_array is None
            or self.kmer_array.size == 0
            or len(self.barcodes) == 0
        ):
            return matches

        k = self.k
        n_seeds = max_mismatches + 1
        seq_len = len(sequence)
        if k is None or k <= 0:
            k = seq_len // n_seeds

        seq_bits = _encode_dna_to_bits(sequence)
        seen = set()

        for i in range(seq_len - k + 1):
            kmer = sequence[i : i + k]
            kmer_int = self._kmer_to_int(kmer)

            idx_matches = np.where(self.kmer_array == kmer_int)[0]
            for idx in idx_matches:
                start = self.barcode_offsets[idx]
                end = self.barcode_offsets[idx + 1]

                for barcode_index in self.barcode_indices[start:end]:
                    if barcode_index in seen:
                        continue
                    seen.add(barcode_index)

                    candidate = self.barcodes[barcode_index]
                    if (
                        _hamming_distance_bits(
                            seq_bits, _encode_dna_to_bits(candidate)
                        )
                        <= max_mismatches
                    ):
                        matches.add(candidate)

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

        self.index_type = None
        self.kmer_index = {}
        self.kmer_length = kmer_length
        self.barcode_length = None
        self.jitter = jitter
        self.mismatches = mismatches
        self.logger = setup_worker_logger()

        if pickle_file is not None and os.path.exists(pickle_file):
            self._load_trie(pickle_file, mismatches = mismatches)
            if index_type is None:
                self.index_type = self.index_type  # loaded from pickle
        else:
            self._build_trie(barcode_sequences)
            self._build_index(barcode_sequences, index_type=index_type or "seed")
            self._save_trie(pickle_file)

        self.report_memory_usage()

    def summarize(self, limit: int = 3):
        """
        Log a summary of the loaded automata and indexes, matching inspect_pkl output.
        """
        logger = self.logger

        logger.info("=== BarcodeMatcherAhoCorasick Summary ===")

        # Top-level metadata
        logger.info(f"Index type: {self.index_type}")
        logger.info(f"Mismatches: {self.mismatches}")
        logger.info(f"Barcode length: {self.barcode_length}")
        logger.info(f"kmer length: {self.kmer_length}")

        # Version is only known if loaded from pickle
        version = getattr(self, "version", None)
        logger.info(f"Pickle version: {version}")

        # Automata
        logger.info(f"Whitelists in automata: {list(self.automata.keys())}")

        # Indexes
        for whitelist_key, idx_obj in self.kmer_index.items():
            logger.info(f"Whitelist: {whitelist_key}")

            # Barcodes
            barcodes = getattr(idx_obj, "barcodes", [])
            logger.info(f"  Number of barcodes: {len(barcodes)}")

            # Index type
            if isinstance(idx_obj, SeedIndex):
                index_type = "seed"
                has_seed = bool(getattr(idx_obj, "seed_index", None))
                kmer_array = []
            elif isinstance(idx_obj, CompactKmerIndex):
                index_type = "kmer"
                has_seed = False
                kmer_array = (
                    idx_obj.kmer_array.tolist()
                    if idx_obj.kmer_array is not None
                    else []
                )
            else:
                index_type = "unknown"
                has_seed = False
                kmer_array = []

            logger.info(f"  Index type stored: {index_type}")
            logger.info(f"  Has seed_index?: {'Yes' if has_seed else 'No'}")

            # Automaton preview
            automaton = self.automata.get(whitelist_key)
            if automaton is None:
                trie_preview = []
            elif isinstance(automaton, dict):
                trie_preview = list(automaton.keys())[:limit]
            else:
                # e.g. pyahocorasick Automaton
                try:
                    trie_preview = [k for k, _ in list(automaton.items())[:limit]]
                except Exception:
                    trie_preview = []

            logger.info(
                f"  Automaton (first {limit} keys/values): "
                f"{trie_preview}{'...' if automaton and len(trie_preview) == limit else ''}"
            )

            # Index-specific preview
            if isinstance(idx_obj, CompactKmerIndex):
                kmer_array = (
                    idx_obj.kmer_array.tolist()
                    if idx_obj.kmer_array is not None
                    else []
                )
                preview = kmer_array[:limit]
                logger.info(
                    f"  k-mer array (first {limit} values): "
                    f"{preview}{'...' if len(kmer_array) > limit else ''}"
                )

            elif isinstance(idx_obj, SeedIndex):
                seed_index = idx_obj.seed_index or {}
                seed_lengths = sorted(seed_index.keys())
                preview_lengths = seed_lengths[:limit]

                logger.info(
                    f"  Seed index lengths (first {limit}): "
                    f"{preview_lengths}{'...' if len(seed_lengths) > limit else ''}"
                )

                if preview_lengths:
                    L = preview_lengths[0]
                    seeds = list(seed_index[L].keys())[:limit]
                    logger.info(
                        f"  Example seeds for length {L}: "
                        f"{seeds}{'...' if len(seed_index[L]) > limit else ''}"
                    )

        logger.info("=== End Summary ===")

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
                # Memory usage of barcode_offsets (if finalized)
                if hasattr(kmer_index, "barcode_offsets"):
                    kmer_index_memory += asizeof.asizeof(kmer_index.barcode_offsets)
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
            first_seq = next(iter(sequences))
            self.barcode_length = len(first_seq)  # <-- set here
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
            self.barcode_length = L
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
        Load a pre-built Aho-Corasick trie and k-mer/seed indexes from a pickle file.
        Supports multiple pickle files without overwriting existing data.
        Rebuilds indexes if the mismatch budget differs from the stored value.
        """
        self.logger.info(f"Importing pickle from '{pickle_file}'")

        if not os.path.exists(pickle_file):
            raise FileNotFoundError(f"Pickle file not found: {pickle_file}")

        try:
            with gzip.open(pickle_file, "rb") as f:
                loaded_data = pickle.load(f)

            # Stored metadata
            stored_index_type = loaded_data.get("index_type")
            stored_mismatches = loaded_data.get("index_mismatches", 0)
            stored_barcode_length = loaded_data.get("barcode_length", None)
            top_level_kmer_length = loaded_data.get("kmer_length", None)
            version = loaded_data.get("version", "unknown")

            self.logger.info(f"Pickle version: {version}")
            self.logger.info(f"Stored index type: {stored_index_type}")
            self.logger.info(f"Stored mismatches: {stored_mismatches}")
            self.logger.info(f"Stored barcode length: {stored_barcode_length}")

            # Set mismatch budget (requested by caller)
            self.mismatches = mismatches

            # Load automata (merge with existing)
            automaton = loaded_data.get("automata", {})
            for whitelist_key, trie in automaton.items():
                self.automata[whitelist_key] = trie
            self.logger.info(f"Loaded automata whitelists: {list(automaton.keys())}")

            # Load k-mer / seed indexes
            kmer_index_data = loaded_data.get("kmer_index", {})
            for whitelist_key, ki_data in kmer_index_data.items():
                barcodes = ki_data.get("barcodes", [])

                # Determine if rebuild is needed
                need_rebuild = self.mismatches != stored_mismatches

                # Handle SeedIndex
                if stored_index_type == "seed":
                    if need_rebuild:
                        idx_obj = SeedIndex(max_mismatches=self.mismatches)
                        idx_obj.build_from_barcodes(barcodes, logger=self.logger)
                        self.index_type = "seed"
                        if barcodes:
                            self.barcode_length = len(barcodes[0])
                        self.kmer_length = None
                    else:
                        idx_obj = SeedIndex(max_mismatches=self.mismatches)
                        idx_obj.barcodes = barcodes
                        idx_obj.seed_index = ki_data.get("seed_index", None)

                # Handle CompactKmerIndex
                elif stored_index_type == "kmer":
                    # Determine k-mer length for this whitelist
                    k = top_level_kmer_length
                    if k is None:
                        if barcodes:
                            k = max(1, len(barcodes[0]) // (self.mismatches + 1))
                        else:
                            raise ValueError(
                                f"No barcodes found for whitelist '{whitelist_key}' to derive k-mer length"
                            )

                    if need_rebuild:
                        idx_obj = CompactKmerIndex(k=k)
                        for bc in barcodes:
                            idx_obj.add_barcode(bc)
                        idx_obj.finalize()
                        self.index_type = "kmer"
                        if barcodes:
                            self.barcode_length = len(barcodes[0])
                        self.kmer_length = k
                    else:
                        idx_obj = CompactKmerIndex(k=k)
                        idx_obj.barcodes = barcodes
                        arr = ki_data.get("kmer_array", [])
                        idx_obj.kmer_array = np.array(
                            arr, dtype=np.uint32 if k > 8 else np.uint16
                        )
                        idx_obj.barcode_indices = np.array(
                            ki_data.get("barcode_indices", []), dtype=np.int32
                        )
                        idx_obj.barcode_offsets = np.array(
                            ki_data.get("barcode_offsets", []), dtype=np.int64
                        )
                else:
                    raise ValueError(f"Unknown index_type in pickle: {stored_index_type}")

                # Store index against its whitelist
                self.kmer_index[whitelist_key] = idx_obj

            # Optionally save rebuilt indexes once
            if need_rebuild:
                self.logger.info(f"Rebuilt indexes with mismatches={self.mismatches}, saving back to pickle")
                self._save_trie(pickle_file)

            self.logger.info(f"Loaded k-mer/seed indexes for whitelists: {list(self.kmer_index.keys())}")

        except (ImportError, FileNotFoundError, ValueError) as e:
            self.logger.error(f"Error loading trie and index from {pickle_file}: {e}")
            raise

    def _save_trie(self, pickle_file: str):
        """
        Save the Aho-Corasick automata and k-mer/seed indexes to a pickle file.
        Preserves all whitelist keys and automata.
        """
        self.logger.info(f"Pickling and compressing automata and indexes to '{pickle_file}'")

        kmer_index_serializable = {}

        for whitelist_key, idx_obj in self.kmer_index.items():

            if isinstance(idx_obj, SeedIndex):
                index_type = "seed"
                barcodes = idx_obj.barcodes
                seed_index = idx_obj.seed_index  # store seed index if already built
                kmer_array = []
                barcode_indices = []
                barcode_offsets = []

            elif isinstance(idx_obj, CompactKmerIndex):
                index_type = "kmer"
                barcodes = idx_obj.barcodes
                seed_index = None
                kmer_array = idx_obj.kmer_array.tolist() if idx_obj.kmer_array is not None else []
                barcode_indices = idx_obj.barcode_indices.tolist() if idx_obj.barcode_indices is not None else []
                barcode_offsets = idx_obj.barcode_offsets.tolist() if idx_obj.barcode_offsets is not None else []

            else:
                raise TypeError(f"Unknown index object type: {type(idx_obj)}")

            kmer_index_serializable[whitelist_key] = {
                "index_type": index_type,
                "barcodes": barcodes,
                "seed_index": seed_index,
                "kmer_array": kmer_array,
                "barcode_indices": barcode_indices,
                "barcode_offsets": barcode_offsets,
            }

        data_to_save = {
            "version": __version__,
            "automata": self.automata,               # preserve all automata
            "kmer_index": kmer_index_serializable,
            "kmer_length": self.kmer_length,
            "index_type": self.index_type,           # "seed" or "kmer"
            "index_mismatches": self.mismatches,     # mismatches used
            "barcode_length": self.barcode_length,
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
                start = idx.barcode_offsets[i]
                end = idx.barcode_offsets[i + 1]
                barcode_ids = idx.barcode_indices[start:end][:num_values]
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
            self.logger.warning(
                f"No index found for whitelist: '{whitelist_key}'. "
                f"Available keys: {list(self.kmer_index.keys())}"
            )
            return []

        idx_obj = self.kmer_index[whitelist_key]
        if isinstance(idx_obj, SeedIndex):
            candidates = idx_obj.query(sequence)
        elif isinstance(idx_obj, CompactKmerIndex):
            candidates = idx_obj.query(sequence, max_mismatches=self.mismatches)
        else:
            self.logger.warning(f"Unknown index type for whitelist: {whitelist_key} (object type: {type(idx_obj)})")
            return []

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

def parser_pickle(parser):
    subparser = parser.add_parser(
        "pickle",
        description="""
Generate and pickle to file an Aho-Corasick Trie and seed or k-mer index.

Example:

scarecrow pickle --barcodes BC1:v1:barcodes1.txt --index seed --mismatches 2
---
""",
        epilog="If the out file exists then this will be loaded rather than re-generated.",
        help="Generate and pickle an Aho-Corasick trie and seed or k-mer index",
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


def validate_pickle_args(parser, args):
    """
    Validate arguments
    """
    # Global logger setup
    logfile = "{}_{}.{}".format(
        "./scarecrow_pickle", generate_random_string(), "log"
    )
    logger = setup_logger(logfile)
    logger.info(f"scarecrow version {__version__}")
    logger.info(f"logfile: '{logfile}'")

    run_pickle(
        barcodes = args.barcodes,
        force_overwrite = args.force_overwrite,
        index = args.index,
        mismatches = args.mismatches,
    )


@log_errors
def run_pickle(
    barcodes: str = None,
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
        pickle_file = f"{file_path}.{index}.m{mismatches}.pkl.gz"
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
