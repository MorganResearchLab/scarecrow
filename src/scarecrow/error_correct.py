#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Cell barcode error correction from a whitelist
"""

import numpy as np
import re
import logging

# inputs: data container of whitelist
# FASTQ read sequence
# correct barcode position
def make_barcode(a_mat, t_mat, c_mat, g_mat, bc_idx):
    '''
    Given a set of binary numpy matrices that map barcode sequence 
    positions to nucleotides, reconstruct a specific barcode
    '''

    bcode = np.array(["N" for i in range(a_mat.shape[1])])
    bcode[a_mat[best_match, :] == 1] = "A"
    bcode[t_mat[best_match, :] == 1] = "T"
    bcode[c_mat[best_match, :] == 1] = "C"
    bcode[g_mat[best_match, :] == 1] = "G"

    return("".join(bcode))


def precompute_matrices(whitelist_path, bc_len):
    '''
    Given a path to a text file containing a cell barcode sequence per row,
    parse the file into 5 matrices - one per nucleotide + one for N.

    input:
        whitelist_path: str - absolute path to file of cell barcode whitelist

    output:
        mat_dict: dict - dictionary with 5 elements corresponding to A, T, C, G, N
        each of which contains a binary matrix of whitelist X poition
    '''

    # get number of BCs in whitelist quickly
    with open(whitelist_path, "r") as f:
        nbc = sum(1 for _ in f)
    
    # storage matrices
    a_mat = np.zeros((nbc, bc_len))
    t_mat = np.zeros((nbc, bc_len))
    c_mat = np.zeros((nbc, bc_len))
    g_mat = np.zeros((nbc, bc_len))

    # nt regex
    a_re = re.compile("A")
    t_re = re.compile("T")
    c_re = re.compile("C")
    g_re = re.compile("G")

    ibc = 0
    with open(whitelist_path, "r") as wfile:
        for bcline in wfile.read_lines():
            bc = bc.line.rstrip("\n")
            
            for match in a_re.finditer(bc):
                a_mat[ibc, match.span()[0]] = 1

            for match in t_re.finditer(bc):
                t_mat[ibc, match.span()[0]] = 1

            for match in c_re.finditer(bc):
                c_mat[ibc, match.span()[0]] = 1

            for match in g_re.finditer(bc):
                g_mat[ibc, match.span()[0]] = 1

            ibc += 1

    mat_dict = {"A": a_mat, "T": t_mat, "C": c_mat,
                "G": g_mat}

    return(mat_dict)


def correct_barcode(sequence, whitelist_dict, edit):
    '''
    Given a dictionary of matrices that map nucleotide 
    barcode positions to nucleotides, find the closest 
    matching barcode(s) within the edit distance.

    '''

    ## Preamble
    # nt regex
    a_re = re.compile("A")
    t_re = re.compile("T")
    c_re = re.compile("C")
    g_re = re.compile("G")
    n_re = re.compile("N")

    bc_len = len(sequence)
    exp_len = whitelist_dict["A"].shape[1]

    if bc_len != exp_len:
        IOError("Wrong barcode sequence length. Expected {} nts, observed {}.".format(exp_len, bc_len))
    
    bc_a_mat = np.zeros((1, bc_len))
    bc_t_mat = np.zeros((1, bc_len))
    bc_c_mat = np.zeros((1, bc_len))
    bc_g_mat = np.zeros((1, bc_len))
    bc_n_mat = np.zeros((1, bc_len))
    
    for match in a_re.finditer(bc):
        bc_a_mat[0, match.span()[0]] = 1

    for match in t_re.finditer(bc):
        bc_t_mat[0, match.span()[0]] = 1

    for match in c_re.finditer(bc):
        bc_c_mat[0, match.span()[0]] = 1

    for match in g_re.finditer(bc):
        bc_g_mat[0, match.span()[0]] = 1

    for match in n_re.finditer(bc):
        bc_n_mat[0, match.span()[0]] = 1

    # if the nucleotides match then they score > 1
    a_edit = (bc_a_mat + a_mat > 1).sum(axis=1)
    t_edit = (bc_t_mat + t_mat > 1).sum(axis=1)
    c_edit = (bc_c_mat + c_mat > 1).sum(axis=1)
    g_edit = (bc_g_mat + g_mat > 1).sum(axis=1)

    # how do we handle N nucleotides? Assume this increases the edit distance?
    tot_edit = bc_len - (a_edit + t_edit + c_edit + g_edit)
    best_match = np.argmin(tot_edit)

    # check if within edit distance threshold
    if len(best_match) > 1:
        if any(tot_edit[best_match] < edit):
            choices = tot_edit[best_match] < edit
            bc_correct = random.sample(best_match[choices], 1)
        else:
            logging.warning("Ambiguous barcode assignment - skipping")
    else:
        bc_correct = best_match

    # compute the barcode from the matrices
    new_bc = make_barcode(a_mat, t_mat, c_mat, g_mat, bc_correct)
        
    return(new_bc)

    
