#!/usr/bin/python3

import pickle
import argparse
import csv

from functions import loadCounts, findErrors
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.Data import CodonTable
from Bio.Alphabet import IUPAC

import random
import matplotlib.pyplot as pyplot

def set_up_total(reference):
    valid_counts = loadCounts(reference, '.counts.p')
    total = {}
    for t in valid_counts.keys():
        total[t] = {}
        for e in valid_counts[t].keys():
            total[t][e] = {'counts': {}, 'depth': {}, 'protein': '', 'exp_activity': '', 'pred_activity': ''}
    return total

def get_sequencing_data(total):
    # input file need: path to counts file, fraction name
    with open(args.files, 'r') as f:
        for line in  f.readlines():
            # name = H/N/L/MM, denotes activity fraction
            name, count_loc, depth1_loc, depth2_loc = line.rstrip().split(',')
            depth = calculate_depth(depth1_loc, depth2_loc)
            # load individual counts
            with open(count_loc, 'rb') as p:
                counts = pickle.load(p)

            # read all data into t
            for t in counts.keys():
                for errors, number in counts[t].items():
                    if len(errors) == 1:
                        continue
                    total[t][errors]['depth'][name] = average_depth(errors, depth)
                    total[t][errors]['counts'][name] = number
                    total[t][errors]['protein'] = error_to_protein(errors)

                    if errors in experimental:
                        total[t][errors]['exp_activity'] = experimental[errors]['activity']
                    else:
                        total[t][errors]['exp_activity'] = None

    # once all fractions are in total, compare H/MM/L/N counts to get an overall prediction
    for m in total.keys():
        for errors in total[m].keys():
            total[m][errors]['pred_activity'] = predict_activity(errors, total)

    return total

def get_experimental_data():
    # If experimental data is provided, collect that data into a dictionary
    # experimental[errors]{'activity': float, 'protein': string}
    # The file is CSV format with header, columns give start & end of deletion, activity (float) and protein effect
    experimental = {}
    rejected = defaultdict(int)
    codons = list(CodonTable.unambiguous_dna_by_name["Standard"].forward_table.keys())
    codons += CodonTable.unambiguous_dna_by_name["Standard"].stop_codons

    if not args.experimental:
        return experimental

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
            errors = findErrors(read, ref, rejected, codons, MAX_ERROR_INDEX=720)
            experimental[errors] = {'activity': a, 'protein': line.get(protein)}

    return experimental

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

def error_to_protein(errors):
    protein = []
    for i in range(1, len(errors), 3):
        dna_pos = int(errors[i])
        dna_ref = Seq(errors[i+1], alphabet=IUPAC.ambiguous_dna)
        dna_mut = Seq(errors[i+2], alphabet=IUPAC.ambiguous_dna)
        try:
            prot_pos = str(int(dna_pos / 3 + 1))
            prot_ref = str(dna_ref.translate())
            prot_mut = str(dna_mut.translate())
        except CodonTable.TranslationError:
            if dna_mut == '---':
                prot_mut = 'Δ'
        if prot_ref == prot_mut:
            continue
        protein.append(prot_ref + prot_pos + prot_mut)

    return '/'.join(protein)

def predict_activity(errors, total):
    """use entries in total[errors[0]][errors] to predict activity
    - pred_activity: H/M/L, whichever has the highest frequency of observing this mutation
    - confidence: 0 for N<5, 1 for 6 <= N < 10, 2 for 11 <= N < 20, 3 for 21 <= N < 50, 4 for <= N < 100, 5 for 100 <= N
    """
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


def classify_point_protein(expected, actual):
    """
    Assume errors is a well-behaved mutation, that is s/d only
    :param expected:
    :param actual:
    :return: string, 's' or 'd'
    """
    dna_ref = Seq(expected, alphabet=IUPAC.ambiguous_dna)
    dna_mut = Seq(actual, alphabet=IUPAC.ambiguous_dna)
    prot_ref = str(dna_ref.translate())

    if str(dna_mut) == '---':
        return 'd'
    else:
        prot_mut = str(dna_mut.translate())

    if prot_ref != prot_mut:
        return 's'
    else:
        return ''


def classify_protein(errors):
    protein_errors = []
    for i in range(1, len(errors), 3):
        protein_errors.append(classify_point_protein(errors[i+1], errors[i+2]))
    prot = ''.join(protein_errors)
    if prot == '':
        return 'wt'
    else:
        return prot

def total_protein(total):
    """
    Total is a large dictionary that contains all sequencing data on DNA level. Some substitutions are between synonymous
    codons and therefore have no effect on amino acids. An 'sd' mutation on DNA might be a 'd' when translated.
    This function goes through the entire DNA dictionary and redoes it based on protein effect.
    :param total: contains all mutations in the format total['d'][error tuple]['pred_activity']
    :return:
    """
    protein  = {k:{} for k in total.keys()}
    protein['wt'] = {}
    for k in total.keys():
        for errors in total[k]:
            prot_key = classify_protein(errors)
            protein[prot_key][errors] = {'protein': total[k][errors]['protein'],
                                         'pred_activity': total[k][errors]['pred_activity']}
    return protein

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reference', help='Reference fasta file', required=True)
    parser.add_argument('-f', '--files', help='File listing counts to be combined')
    parser.add_argument('-e', '--experimental', help='Experimental data on mutations', required=False)
    args = parser.parse_args()

    reference = SeqIO.read(args.reference, 'fasta', alphabet=IUPAC.ambiguous_dna)
    total = set_up_total(reference)

    experimental = get_experimental_data()
    total = get_sequencing_data(total)
    protein = total_protein(total)

    points = generate_experimental_points()

    plot_known_mutations(points)

    for k in ('d', 'dd', 'ddd', 'sd', 'sdd', 'sddd'):
        protein[k]['pymol'] = print_pymol(k, protein)
        for frac, resi in protein[k]['pymol'].items():
            print(k, frac, '+'.join(resi))

