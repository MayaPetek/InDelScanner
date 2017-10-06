#!/usr/bin/python3

import pickle
import argparse
import random
import sys
import csv

from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.Data import CodonTable
from Bio.Alphabet import IUPAC

import matplotlib.pyplot as pyplot

sys.path.insert(0, "/home/maya/Install/Acids")
from indels.ind import set_up_total, mutation_to_protein_notation, classify_mutation, find_dna_mutations, \
    calculate_depth, average_depth

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


def get_sequencing_data(files, total, experimental):
    """
    The input file in args.file lists where counts & depth data is located:
    line = name,counts file loc,depth file loc, depth 2 file loc, activity (H/N/MM/L)
    Read in all mutations, predict activity and place them in total.
    :param total: at first empty dictionary with keys for all expected mutations
    :return: added all known information about expected mutations
    """
    # input file need: path to counts file, fraction name
    with open(files, 'r') as f:
        for line in  f.readlines():
            # name = H/N/L/MM, denotes activity fraction
            name, count_loc, depth1_loc, depth2_loc = line.rstrip().split(',')
            depth = calculate_depth(depth1_loc, depth2_loc)
            # load individual counts
            with open(count_loc, 'rb') as p:
                counts = pickle.load(p)

            # read all data into t
            for errors in counts.keys():
                if len(errors) == 0:
                    continue
                total[errors]['depth'][name] = average_depth(errors, depth)
                total[errors]['counts'][name] = counts[errors]
                total[errors]['protein_mutation'] = mutation_to_protein_notation(errors)
                total[errors]['dna_type'] = classify_mutation(errors, style='dna')
                total[errors]['protein_type'] = classify_mutation(errors, style='protein')
                if errors in experimental:
                    total[errors]['exp_activity'] = experimental[errors]['activity']
                else:
                    total[errors]['exp_activity'] = None
    return total


def predict_activity(total):
    """
    Use entries in total [errors] to predict activity
    - pred_activity: H/M/L, whichever has the highest frequency of observing this mutation
    - confidence: 0 for N<5, 1 for 6 <= N < 10, 2 for 11 <= N < 20, 3 for 21 <= N < 50, 4 for <= N < 100, 5 for 100 <= N
    :param errors: which mutation to look up
    :param total: total = full DNA dictionary, must contain all information at this point
    :return: tuple (H/M/L, confidence) where confidence is an int
    """
    # once all fractions are in total, compare H/MM/L/N counts to get an overall prediction
    for m in total.keys():
        for errors in total[m].keys():
            total[m][errors]['pred_activity'] = predict_activity(errors, total)

    # make a dictionary: enrichment - fraction, eg. 1.80:H, 0.1:MM, 0.0: L
    enrichment = {}
    N_freq = total[errors[0]][errors]['counts']['N'] / total[errors[0]][errors]['depth']['N']
    if N_freq == 0:
        pred_activity = ''
    else:
        for f in ('H', 'MM', 'L'):
            frac_freq = total[errors[0]][errors]['counts'][f] / total[errors[0]][errors]['depth'][f]
            e = frac_freq / N_freq
            enrichment[e] = f
        max_activity = max(enrichment.keys())
        if max_activity == 0 :
            pred_activity = ''
        else:
            pred_activity = enrichment[max_activity]

    N = total[errors[0]][errors]['counts']['N']
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

    return (pred_activity, pred_confidence)

def generate_experimental_points():
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
        pred_activity = total[e[0]][e]['pred_activity']
        if pred_activity[0] != '':
            if pred_activity[0] == 'L':
                points['x'].append(random.random())
            else:
                points['x'].append(random.uniform(0.4, 0.6))

            points['y'].append(total[e[0]][e]['exp_activity'])
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

    if errors[0] in ('d', 'dd', 'ddd'):
        # take first reported deletion
        assert errors[3] == '---'
        prot_pos = str(int(int(errors[1])/3) + 1)
        return prot_pos
    elif errors[0] in ('sd', 'sdd', 'sddd'):
        # take first reported deletion, that is second mutation
        assert errors[6] == '---'
        prot_pos = str(int(int(errors[4])/3) + 1)
        return prot_pos
    elif errors[0] in ('s', 'ss'):
        # report substitution residues
        residues = []
        for i in range(1, len(errors), 3):
            dna_pos = int(errors[i])
            dna_ref = Seq(errors[i+1], alphabet=IUPAC.ambiguous_dna)
            dna_mut = Seq(errors[i+2], alphabet=IUPAC.ambiguous_dna)

            prot_pos = str(int(dna_pos / 3 + 1))
            prot_ref = str(dna_ref.translate())
            if dna_mut == '---':
                prot_mut = 'Δ'
            else:
                prot_mut = str(dna_mut.translate())

            if prot_ref == prot_mut:
                continue

            residues.append(prot_pos)
        return '+'.join(residues)

def print_pymol(fraction, dictionary):
    """
    Return a dictionary with a list of residues for Pymol plotting.
    """
    f = {'H': [], 'MM':[], 'L': []}

    for e in dictionary[fraction].keys():
        a = dictionary[fraction][e]['pred_activity']
        if a[0] != '' and a[1] > 1:
            if pymol_residues(e) != '':
                f[a[0]].append(pymol_residues(e))

    for frac, resi in f.items():
        f[frac] = sorted(set(resi))

    return f


def sub_that_change_del(protein):
    """
    Look for substitutions that change the effect of deletions.
    Go through all errors that on protein level are classified as 'sd', 'sdd' or 'sddd' and compare their effect with
    corresponding pure deletion. If there is a difference, report:
    - which deletion & which combined substitution-deletion
    - predicted activity (incl. confidence) for both
    :param protein: dictinary listing all mutations by effect on protein. protein['d'][errors]['protein'/'pred_activity']
    :return: interesting mutations
    """

    sd_dif = {}

    for prot_key in ('sd', 'sdd', 'sddd'):
        for error in protein[prot_key].keys():
            del_error = list(error)
            sub_error = tuple(['s'] + del_error[1:4])
            del_error = tuple([error[0][1:]] + del_error[4:])

            act_sd = protein[prot_key][error]['pred_activity']
            act_d = protein[del_error[0]][del_error]['pred_activity']
            act_s = protein['s'][sub_error]['pred_activity']

            if act_sd[1] < 3 or act_d[1] < 3:
                continue
            elif act_sd[0] != act_d[0]:

                print('protein', total[error[0]][error]['protein'], act_sd)
                print('counts', total[error[0]][error]['counts'])

                print('protein', total[del_error[0]][del_error]['protein'], act_d)
                print('counts', total[del_error[0]][del_error]['counts'])

                print('protein', total['s'][sub_error]['protein'], act_s)
                print('counts', total['s'][sub_error]['counts'])

                print(error, del_error, sub_error)
                print()

    return sd_dif

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reference', help='Reference fasta file', required=True)
    parser.add_argument('-f', '--files', help='File listing counts to be combined')
    parser.add_argument('-e', '--experimental', help='Experimental data on mutations', required=False)
    args = parser.parse_args()

    reference = SeqIO.read(args.reference, 'fasta', alphabet=IUPAC.ambiguous_dna)
    total = set_up_total(reference)

    experimental = get_experimental_data(args, reference)
    total = get_sequencing_data(args.files, total, experimental)

    n = 0
    for errors in total.keys():
        if total[errors]['counts'] != {}:
            n += 1
    print('Found', n, 'mutations')

    # points = generate_experimental_points()
    # plot_known_mutations(points)
    #
    # for k in ('d', 'dd', 'ddd', 'sd', 'sdd', 'sddd'):
    #     protein[k]['pymol'] = print_pymol(k, protein)
    #     for frac, resi in protein[k]['pymol'].items():
    #         print(k, frac, '+'.join(resi))

    #sub_that_change_del(protein)