#!/usr/bin/python3

import argparse
import random
import sys
import pickle
import csv

from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC

import numpy as np
import matplotlib.pyplot as plt
from sklearn import decomposition

sys.path.insert(0, "/home/maya/Install/Acids")
from indels.ind import classify_mutation, find_dna_mutations, get_total, classify_point_protein


def get_experimental_data(args, reference):
    """
    If experimental data is provided, collect that data into a dictionary
    The file is in CSV format with header, columns give start & end of deletion, activity (float) and protein effect
    If no data is provided, return empty dictionary.
    :param args: arguments processed with argparse, must include args.experimentsl (a csv file)
    :param reference: SeqRecord of reference fasta file
    :return: experimental[errors]{'activity': float, 'protein': string}
    """

    exp_data = {}
    rejected = defaultdict(int)

    if not args.experimental:
        return exp_data

    with open(args.experimental, 'r') as f:
        lit = csv.DictReader(f, dialect='excel')
        ref = reference.seq.upper().tomutable()
        r = str(ref)
        start, end, activity, protein = lit.fieldnames

        for line in lit:
            s = int(line.get(start))
            e = int(line.get(end))
            a = float(line.get(activity))
            length = e - s
            read = MutableSeq(r[:s] + ("-" * length) + r[e:], ref.alphabet)
            errors = find_dna_mutations(read, ref, rejected, MAX_ERROR_INDEX=720)
            exp_data[errors] = {'activity': a, 'protein': line.get(protein)}
    print("Loaded experimental data")

    return exp_data


def predict_activity(total, style='depth', position_counts=None):
    """
    Use entries in total [errors] to predict activity
    - pred_activity: H/M/L, whichever has the highest frequency of observing this mutation
    - confidence: 0 for N<5, 1 for 6 <= N < 10, 2 for 11 <= N < 20, 3 for 21 <= N < 50, 4 for <= N < 100, 5 for 100 <= N
    :param total: total = full DNA dictionary, must contain all information at this point
    :return: tuple (H/M/L, confidence) where confidence is an int
    """

    for e in total.keys():
        if e == ():
            continue
        # make an inverted dictionary: enrichment - fraction, eg. 1.80:H, 0.1:MM, 0.0: L
        enrichment = calculate_enrichment(e, total, position_counts=position_counts, style=style)
        inverted = {}
        if enrichment == {}:
            pred_activity = ''
        else:
            for f, value in enrichment.items():
                inverted[value] = f
            max_activity = max(inverted.keys())
            if max_activity == 0 :
                pred_activity = ''
            else:
                pred_activity = inverted[max_activity]


        # Add a confidence indicator to this prediction
        N = total[e]['counts']['N']

        if N > 100:
            pred_confidence = 5
        elif N > 50:
            pred_confidence = 4
        elif N > 20:
            pred_confidence = 3
        elif N > 5:
            pred_confidence = 2
        elif N > 0:
            pred_confidence = 1
        else:
            pred_confidence = 0

        if pred_confidence == 0:
            pred_activity = ''

        total[e]['pred_activity'] = (pred_activity, pred_confidence)

    return total


def generate_experimental_points(experimental, total):
    """
    Convert experimental data into points for pyplot.
    :return: points, dictionary containing three lists: x, y coordinates & colour
    """

    # set color maps
    green_map = plt.cm.get_cmap('Greens')
    yellow_map = plt.cm.get_cmap('Blues')
    red_map = plt.cm.get_cmap('Reds')
    # confidence is an integer in range [0, 1, 2, 3, 4, 5]
    colors = {'H': [green_map(i/10) for i in range(3, 9, 1)],
              'MM': [yellow_map(i/10) for i in range(3, 9, 1)],
              'L': [red_map(i/10) for i in range(3, 9, 1)]}
    # make a list of points with exp & pred activity
    points = {'x': [], 'y': [], 'c': []}

    for e in experimental.keys():
        pred_activity = total[e]['pred_activity']
        if pred_activity[0] != '':
            if pred_activity[0] == 'L':
                points['x'].append(random.random())
            else:
                points['x'].append(random.uniform(0.4, 0.6))

            points['y'].append(total[e]['exp_activity'])
            points['c'].append(colors[pred_activity[0]][pred_activity[1]])

    return points


def plot_known_mutations(points):
    """
    Scatterplot mutations with known activity
    x - random number (0,1)
    y - experimental fluorescence, log scale
    c - colour, red = low, blue = medium, green = high. Intensity indicates confidence.
    """
    plt.figure()
    plt.scatter(points['x'], points['y'], c=points['c'], marker='.')
    plt.yscale('symlog', linthreshy=0.02)
    plt.axhline(y=0.05, color='k', linewidth=1)
    plt.axhline(y=0.005, color='k', linewidth=1)
    plt.ylabel('Log % eGFP fluorescence')
    plt.xticks([])
    return


def pymol_helper(errors):
    """
    Which residues in a mutation actually affect the protein. Deliberately ignores synonymous mutations.
    Reports 1st amino acid for d/dd/ddd and 2nd amino acid for sd/sdd/sddd
    :param errors: tuple describing codon changes
    :return: pymol-select compatible description of residues
    """
    prot_effect = classify_mutation(errors, style='protein')
    if prot_effect in ('d', 'dd', 'ddd'):
        # take first reported deletion
        prot_pos = str(int(int(errors[0])/3) + 1)
        return prot_pos
    elif prot_effect in ('sd', 'sdd', 'sddd'):
        # take first reported deletion, that is second mutation
        prot_pos = str(int(int(errors[3])/3) + 1)
        return prot_pos
    elif prot_effect in ('s', 'ss'):
        # report substitution residues
        residues = []
        for i in range(0, len(errors), 3):
            if classify_point_protein(errors[1], errors[2]) == 's':
                dna_pos = int(errors[i])
                prot_pos = str(int(dna_pos / 3 + 1))
                residues.append(prot_pos)
                return '+'.join(residues)


def pymol_resi_list(fraction, total):
    """
    Return a dictionary with a list of residues for Pymol plotting.
    Fraction describes the mutation type on protein level, eg. 'sd'
    """
    fractions = {'H': [], 'MM': [], 'L': []}

    for e in total.keys():
        if e == ():
            continue
        if total[e]['protein_type'] == fraction:
            # add location to f, a[0] is prediction of activity, a[1] is confidence
            a = total[e]['pred_activity']
            if a[0] != '' and a[1] > 1:
                if pymol_helper(e) != '':
                    fractions[a[0]].append(pymol_helper(e))

    for frac, resi in fractions.items():
        fractions[frac] = sorted(set(resi))

    return fractions


def find_epistatic_combinations(total, cutoff, output):
    """
    Look for substitutions that change the effect of deletions.
    Go through all errors that on protein level are classified as 'sd', 'sdd' or 'sddd' and compare their effect with
    corresponding pure deletion. If there is a difference, report:
    - which deletion & which combined substitution-deletion
    - predicted activity (incl. confidence) for both
    :param total: dictionary listing all mutations. total[e]['protein_type'/'pred_activity']
    :param cutoff: minimum level of confidence
    :param
    """
    with open(output + '_min_confidence_' + str(cutoff) + '.csv', 'w', newline='') as f:
        epiwriter = csv.writer(f, delimiter='\t')
        i = 1
        # Look at all mutations that combine substitutions and deletions
        for e_sd in total.keys():
            if e_sd == ():
                continue
            if total[e_sd]['protein_type'] not in ('sd', 'sdd', 'sddd'):
                continue

            # generate related substitution and deletion
            e_d = e_sd[3:]
            e_s = e_sd[:3]

            # retrieve predicted activities and require high confidence
            a_sd = total[e_sd]['pred_activity']
            a_d = total[e_d]['pred_activity']
            a_s = total[e_s]['pred_activity']

            for a in (a_sd, a_d, a_s):
                if a[1] < cutoff:
                    break
            else:
                if len(set(k[0] for k in (a_sd, a_d, a_s))) > 1:
                    min_confidence = min(k[1] for k in (a_sd, a_d, a_s))
                    epiwriter.writerow([i, min_confidence, 'SUB+DEL', total[e_sd]['protein_mutation'], a_sd, e_sd, total[e_sd]['counts']])
                    epiwriter.writerow([i, min_confidence, 'DELETION', total[e_d]['protein_mutation'], a_d, e_d, total[e_d]['counts']])
                    epiwriter.writerow([i, min_confidence, 'SUBSTITUTION', total[e_s]['protein_mutation'], a_s, e_s, total[e_s]['counts']])
                    i += 1
    return


def count_mutations_per_position(total):
    """
    How many mutations are observed per position?
    - s or d: add 1 to first nucleotide of errors. So (3, AAA, TTT) -> counts['dna'][3] +=1
    - longer: (i, AAA, TTT, i+3, TTT, ---) -> counts[dna][i] += 1/2, same for i+3
    Same idea for protein effect
    :param total: Contains all observed DNA mutations, incl. dna effect
    :return:
    """
    [e] = random.sample(total.keys(), 1)
    fractions = total[e]['counts'].keys()
    position_counts = {'protein': {f: {} for f in fractions}, 'dna': {f: {} for f in fractions}}

    for e in total.keys():
        if e == ():
            continue

        dna_len = len(total[e]['dna_type']) # number of codons affected
        for f in fractions:
            for i in range(0, dna_len):
                try:
                    position_counts['dna'][f][int(e[3*i])] += total[e]['counts'][f] / dna_len
                except KeyError:
                    position_counts['dna'][f][int(e[3*i])] = total[e]['counts'][f] / dna_len

        # some substitutions have no effect on the protein, ignore them
        if total[e]['protein_type'] == 'wt':
            continue
        prot_len = len(total[e]['protein_type'])
        if prot_len == dna_len:
            start = 0
        else:
            start = 1
        for f in fractions:
            for i in range(start, prot_len):
                try:
                    position_counts['protein'][f][int(float(e[3*i])/3) + 1] += total[e]['counts'][f] / prot_len
                except KeyError:
                    position_counts['protein'][f][int(float(e[3*i])/3) + 1] = total[e]['counts'][f] / prot_len

    return position_counts

def calculate_enrichment(e, total, position_counts=None, style='protein_mutations'):
    """
    Helper function to calculate enrichments for all fractions relative to baseline
    :param e: which mutation
    :param total: dictionary with all sequencing data
    :param style: divide by depth or total number of mutations observed at this position
    :return: {'H': float, 'L': float etc.}
    """
    assert style in ('depth', 'protein_mutations', 'dna_mutations')

    if style == 'depth':
        N_freq = total[e]['counts']['N'] / total[e]['depth']['N']
        if N_freq == 0:
            return {}
        enrichments = {}
        for f in total[e]['counts'].keys():
            if f == 'N':
                continue
            frac_freq = total[e]['counts'][f] / total[e]['depth'][f]
            enrichments[f] = frac_freq / N_freq

    else:
        assert position_counts is not None
        if style == 'protein_mutations':
            level = 'protein'
            if total[e]['protein_type'] == 'wt': # eliminate wt single condon substitutions
                return {}
        else:
            level = 'dna'

        # enrichement = (fraction count / fraction total ) / (N count / N total) = frac_freq / N_freq
        # first calculate all single counts

        dna_len = len(total[e]['dna_type'])  # number of codons affected
        prot_len = len(total[e]['protein_type'])
        start = 0
        if (level == 'protein') and (dna_len != prot_len): # the substitution is wt; for 'ss' assign to second codon
            start = 1
        fraction_frequency = {}
        for f in total[e]['counts'].keys():
            frac_count = total[e]['counts'][f]
            frac_total = 0
            for i in range(start, prot_len):
                frac_total += position_counts[level][f][int(float(e[3 * i]) / 3) + 1]
            try:
                fraction_frequency[f] = frac_count / frac_total
            except ZeroDivisionError:
                fraction_frequency[f] = 0

        # convert individual frequencies to enrichment - X compared to naive library
        enrichments = {}
        if 'N' not in fraction_frequency.keys():
            return {}
        if fraction_frequency['N'] == 0:
            return {}

        for f in total[e]['counts'].keys():
            if f == 'N':
                continue
            enrichments[f] = fraction_frequency[f] / fraction_frequency['N']

    return enrichments


def generate_training(experimental, total, position_counts):
    """
    Prepare np arrays for all mutants for which activity is known.
    :param experimental:
    :param total:
    :return:
    """
    no_points = len(experimental.keys())
    # start with lists
    X_training = []
    y_training = []

    classes = {'L': 0, 'MM': 1, 'H': 2}
    order = ['H', 'MM', 'L']

    # add H/L/MM enrichements to X and measured activity to y
    for e in experimental.keys():
        enrich = calculate_enrichment(e, total, position_counts)
        if enrich == {}:
            continue
        entry = [enrich[fr] for fr in order]
        X_training.append(entry)

        a = experimental[e]['activity']
        if a >= 0.5:
            y_training.append(classes['H'])
        elif a > 0:
            y_training.append(classes['MM'])
        elif a == 0:
            y_training.append(classes['L'])

    X_training = np.array(X_training)
    y_training = np.array(y_training)

    return X_training, y_training, order


def pca_experimental(X_training, y_training, order, total):

    pca = decomposition.PCA(n_components=2)  # keep 2 components, they explain 96% of variance
    pca.fit(X_training)  # learn the coordinates
    X_pca = pca.transform(X_training)

    print(order)
    print(pca.components_)
    print()
    print(pca.explained_variance_ratio_)

    plt.figure()
    colors = ['green', 'turquoise', 'red']
    for color, i, target_name in zip(colors, [2, 1, 0], order):
        plt.scatter(X_pca[y_training == i, 0], X_pca[y_training == i, 1], color=color, label=target_name, marker=".")
    plt.legend(loc='best', shadow=False, scatterpoints=1)
    plt.title('Unscaled PCA of experimental data')
    plt.xlabel('PCA 1')
    plt.ylabel('PCA 2')

    return X_pca, y_training, pca, order


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reference', help='Reference fasta file', required=True)
    parser.add_argument('-f', '--files', help='File listing counts to be combined')
    parser.add_argument('-e', '--experimental', help='Experimental data on mutations', required=False)
    parser.add_argument('-o', '--output', help='Name of output file', required=True)
    args = parser.parse_args()

    # Import counts for individual fractions, of which one is 'N', the naive (unsorted) library
    reference = SeqIO.read(args.reference, 'fasta', alphabet=IUPAC.ambiguous_dna)
    experimental = get_experimental_data(args, reference)
    total = get_total(args.files, experimental, reference)

    position_counts = count_mutations_per_position(total)

    X_training, y_training, order = generate_training(experimental, total, position_counts)
    X_exp, y_exp, pca, order = pca_experimental(X_training, y_training, order, total)

    # Crude predictions
    total = predict_activity(total, style='protein_mutations', position_counts=position_counts)
    find_epistatic_combinations(total, 2, args.output)

    # Save results
    with open(args.output + '.total_with_predictions.p', 'wb') as f:
        pickle.dump(total, f)

    # VALIDATION: PLOT PREDICTIONS FOR EXPERIMENTALLY TESTED MUTATIONS
    points = generate_experimental_points(experimental, total)
    plot_known_mutations(points)
    plt.show()

    # PYMOL HELPER SEQUENCE
    pymol_dict = {}
    for k in ('d', 'dd', 'ddd', 'sd', 'sdd', 'sddd'):
        pymol_dict[k] = pymol_resi_list(k, total)
        for frac, resi in pymol_dict[k].items():
            print(k, frac, '+'.join(resi))
