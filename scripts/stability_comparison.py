#!/usr/bin/python3

import pickle
import argparse


def extract_deletions(total):
    """
    Input dictionary is keyed by DNA mutations. Same deletion coming from different protein background will appear as
    two different deletion, therefore appearing unrecognizable. Therefore, for comparison, restructure deletions dict
    to work as normal[dna position] = all other information
    :param total:
    :return: short_deletions
    """
    short_deletions = {'d': {}, 'dd': {}, 'ddd': {}}
    # type, position, activity
    for e in total.keys():
        protein_type = total[e]['protein_type']
        if protein_type not in ('d', 'dd', 'ddd'):
            continue
        protein_position = int(int(e[1]) / 3 + 1)
        short_deletions[protein_type][protein_position] = total[e]['pred_activity']

    return short_deletions


def compare_background(n_short, s_short, cutoff):
    """
    Take output of del_by_position for two background and look for difference
    :param n_short: deletions dictionary for normal background
    :param s_short: deletions dictionary for stable background
    :param cutoff: minimum quality score
    :return:
    """
    dif = {'d': {'resi': [], 'effect': []}, 'dd': {'resi': [], 'effect': []}, 'ddd': {'resi': [], 'effect': []}}

    for k in n_short.keys():
        for pos in n_short[k].keys():
            # look at same position in both normal and short
            if pos not in s_short[k]:
                continue
            elif n_short[k][pos][1] < cutoff or s_short[k][pos][1] < cutoff:
                continue
            # have data in both backgrounds with good confidence, look for differences
            elif n_short[k][pos][0] != s_short[k][pos][0]:
                effect = (pos, n_short[k][pos], s_short[k][pos])
                dif[k]['resi'].append(str(pos))
                dif[k]['effect'].append(effect)
        dif[k]['resi'] = '+'.join(dif[k]['resi'])
    return dif


def print_differences(dif):
    for k in dif.keys():
        print(k, 'eGFP', 'GFP8', sep=':')
        for e in dif[k]['effect']:
            print(':'.join(str(element) for element in e))
        print()


if __name__ == "__main__":
    """
    1. Import two 'total' dictionaries, output of multiple_activity_fractions.py, for normal and stable background
    2. Restructure dictionaries to be coded as 'deletion at position X'
    3. Compare effect of mutation
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--normal', help='Low activity', required=True)
    parser.add_argument('-s', '--stable', help='High activity', required=True)
    args = parser.parse_args()

    # Get two 'total' counts dictionaries
    with open(args.normal, 'rb') as f:
        normal = pickle.load(f)

    with open(args.stable, 'rb') as f:
        stable = pickle.load(f)

    # Go through all mutations and look for mutations where predicted activity differs
    n_short = extract_deletions(normal)
    s_short = extract_deletions(stable)

    dif = compare_background(n_short, s_short, 3)

    print_differences(dif)
