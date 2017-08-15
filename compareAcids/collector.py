#!/usr/bin/python3

import pickle
import argparse

from functions import loadCounts

from Bio import SeqIO
from Bio.Alphabet import IUPAC

def set_up_total(reference):
    valid_counts = loadCounts(reference, '.del_counts.p')
    total = {}
    for t in valid_counts.keys():
        total[t] = {}
        for e in valid_counts[t].keys():
            total[t][e] = {'counts': {}, 'depth': {}}
    return total


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
            count_loc, name = line.rstrip().split(',')
            with open(count_loc, 'rb') as p:
                counts = pickle.load(p)
            for t in counts.keys():
                for e, number in counts[t].items():
                    total[t][e]['counts'][name] = number

    with open(args.output, 'wb') as o:
        pickle.dump(total, o)