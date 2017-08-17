#!/usr/bin/python3

import pickle
import argparse

from functions import loadCounts

from Bio import SeqIO
from Bio.Alphabet import IUPAC

def set_up_total(reference):
    valid_counts = loadCounts(reference, '.counts.p')
    total = {}
    for t in valid_counts.keys():
        total[t] = {}
        for e in valid_counts[t].keys():
            total[t][e] = {'counts': {}, 'depth': {}}
    return total

def calculate_depth(depth_1, depth_2):
    """
    Collect samtools depth output into a list. 2nd column = 1-based position, 3rd column = coverage.
    Samtools gives two separate files for assembled and unassembled reads
    :param depth_file: output of samtools depth, tab delimited
    :return: a list with coverage per position as ints
    """

    depth = []
    with open(depth_1, 'r') as f:
        for line in f.readlines():
            l = line.split()
            depth.append(int(l[2]))
    with open(depth_2, 'r') as f:
            l = line.split()
            i = int(l[1]) - 1
            depth[i] += int(l[2])
    return depth

def average_depth(errors, depth):
    """
    Use errors description to find which part of reference it is located in, average coverage over that sequence.
    :param errors: a tuple describing a set of mutations
    :param depth: a list giving coverage per position
    :return: a float giving average coverage for this mutation
    """
    # to allow non-consecutive mutations
    e_depth = []
    nt = len(errors) - 1
    for i in range(1,nt,3):
        start = int(errors[i])
        e_depth += depth[start:start+3]
    avg = sum(e_depth) / nt
    return avg

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reference', help='Reference fasta file', required=True)
    parser.add_argument('-f', '--files', help='File listing counts to be combined')
    parser.add_argument('-o', '--output', help='Output file name', required=True)
    args = parser.parse_args()

    reference = SeqIO.read(args.reference,'fasta', alphabet=IUPAC.ambiguous_dna)
    total = set_up_total(reference)

    # input file need: path to counts file, fraction name
    with open(args.files, 'r') as f:
        for line in  f.readlines():
            name, count_loc, depth1_loc, depth2_loc = line.rstrip().split(',')
            depth = calculate_depth(depth1_loc, depth2_loc)
            with open(count_loc, 'rb') as p:
                counts = pickle.load(p)

            for t in counts.keys():
                for errors, number in counts[t].items():
                    if len(errors) == 1:
                        continue
                    total[t][errors]['depth'][name] = average_depth(errors, depth)
                    total[t][errors]['counts'][name] = number

    with open(args.output, 'wb') as o:
        pickle.dump(total, o)