#!/usr/bin/python3

import pickle
import argparse
import pprint

from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.Data import CodonTable
from Bio.Alphabet import IUPAC

def del_by_position(protein):
    summary = {'d': {}, 'dd': {}, 'ddd': {}}
    for k in ('d', 'dd', 'ddd'):
        for error in protein[k].keys():
            prot_pos = int(int(error[1]) / 3 + 1)
            summary[k][prot_pos] = protein[k][error]['pred_activity']
    return summary


def compare_ns(n_short, s_short):
    dif = {'d': {'resi': [], 'effect': []}, 'dd': {'resi': [], 'effect': []}, 'ddd': {'resi': [], 'effect': []}}

    for k in n_short.keys():
        for pos in n_short[k].keys():
            # look at same position in both normal and short
            if pos not in s_short[k]:
                continue
            elif n_short[k][pos][1] < 3 or s_short[k][pos][1] < 3:
                continue
            elif n_short[k][pos][0] != s_short[k][pos][0]:
                effect = (pos, n_short[k][pos], s_short[k][pos])
                dif[k]['resi'].append(str(pos))
                dif[k]['effect'].append(effect)
        dif[k]['resi'] = '+'.join(dif[k]['resi'])
    return dif


def print_dif(dif):
    for k in dif.keys():
        print(k, 'eGFP', 'GFP8', sep=':')
        for e in dif[k]['effect']:
            print(':'.join(str(element) for element in e))
        print()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--normal', help='Low activity', required=True)
    parser.add_argument('-s', '--stable', help='High activity', required=True)
    args = parser.parse_args()


    with open(args.normal, 'rb') as f:
        normal = pickle.load(f)

    with open(args.stable, 'rb') as f:
        stable = pickle.load(f)

    n_short = del_by_position(normal['protein'])
    s_short = del_by_position(stable['protein'])

    dif = compare_ns(n_short, s_short)

    print_dif(dif)
