#!/usr/bin/python3

import argparse
import random
import sys
import csv

from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC

import matplotlib.pyplot as pyplot

sys.path.insert(0, "/home/maya/Install/Acids")
from indels.ind import classify_mutation, find_dna_mutations, get_total, classify_point_protein


def get_experimental_data(args, reference):
    """
    If experimental data is provided, collect that data into a dictionary
    The file is in CSV format with header, columns give start & end of deletion, activity (float) and protein effect
    If no data is provided, return empty dictionary.
    :param args: arguments processed with argparse, must include args.experimentsl (a csv file)
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


def predict_activity(total):
    """
    Use entries in total [errors] to predict activity
    - pred_activity: H/M/L, whichever has the highest frequency of observing this mutation
    - confidence: 0 for N<5, 1 for 6 <= N < 10, 2 for 11 <= N < 20, 3 for 21 <= N < 50, 4 for <= N < 100, 5 for 100 <= N
    :param errors: which mutation to look up
    :param total: total = full DNA dictionary, must contain all information at this point
    :return: tuple (H/M/L, confidence) where confidence is an int
    """

    for e in total.keys():
        # make an inverted dictionary: enrichment - fraction, eg. 1.80:H, 0.1:MM, 0.0: L
        enrich = {}
        if e == ():
            continue
        N_freq = total[e]['counts']['N'] / total[e]['depth']['N']

        if N_freq == 0:
            pred_activity = ''
        else:
            for f in total[e]['counts'].keys():
                if f == 'N':
                    continue
                frac_freq = total[e]['counts'][f] / total[e]['depth'][f]
                enr_f = frac_freq / N_freq
                enrich[enr_f] = f
            max_activity = max(enrich.keys())
            if max_activity == 0 :
                pred_activity = ''
            else:
                pred_activity = enrich[max_activity]
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
    green_map = pyplot.cm.get_cmap('Greens')
    yellow_map = pyplot.cm.get_cmap('Blues')
    red_map = pyplot.cm.get_cmap('Reds')
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
    pyplot.figure(1)
    pyplot.scatter(points['x'], points['y'], c=points['c'])
    pyplot.yscale('symlog', linthreshy=0.02)
    pyplot.axhline(y=0.05, color='k', linewidth=1)
    pyplot.axhline(y=0.005, color='k', linewidth=1)
    pyplot.ylabel('Log % eGFP fluorescence')
    pyplot.xticks([])
    return


def pymol_residues(errors):
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
    f = {'H': [], 'MM': [], 'L': []}

    for e in total.keys():
        if e == ():
            continue
        if total[e]['protein_type'] == fraction:
            # add location to f, a[0] is prediction of activity, a[1] is confidence
            a = total[e]['pred_activity']
            if a[0] != '' and a[1] > 1:
                if pymol_residues(e) != '':
                    f[a[0]].append(pymol_residues(e))

    for frac, resi in f.items():
        f[frac] = sorted(set(resi))

    return f


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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reference', help='Reference fasta file', required=True)
    parser.add_argument('-f', '--files', help='File listing counts to be combined')
    parser.add_argument('-e', '--experimental', help='Experimental data on mutations', required=False)
    parser.add_argument('-o', '--output', help='Name of output file', required=True)
    args = parser.parse_args()

    reference = SeqIO.read(args.reference, 'fasta', alphabet=IUPAC.ambiguous_dna)
    experimental = get_experimental_data(args, reference)
    total = get_total(args.files, experimental, reference)
    total = predict_activity(total)

    find_epistatic_combinations(total, 2, args.output)

    # # VALIDATION: PLOT PREDICTIONS FOR EXPERIMENTALLY TESTED MUTATIONS
    # points = generate_experimental_points(experimental, total)
    # plot_known_mutations(points)
    # pyplot.show()

    # # PYMOL HELPER SEQUENCE
    # pymol_dict = {}
    # for k in ('d', 'dd', 'ddd', 'sd', 'sdd', 'sddd'):
    #     pymol_dict[k] = pymol_resi_list(k, total)
    #     for frac, resi in pymol_dict[k].items():
    #         print(k, frac, '+'.join(resi))
