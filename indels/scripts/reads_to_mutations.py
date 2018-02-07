#!/usr/bin/python3

import sys
import pickle
import re
import os
import csv
import argparse

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from sklearn import decomposition
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

from collections import defaultdict
from functools import partial

from Bio import AlignIO
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq, translate

from indels.ind import trim_read, findEnds, endMatch, findGap, gapAlign, ab1_to_fastq, needle_align
from indels.output import print_coloured_diff

# Demand Python 3.
if sys.version_info[0] < 3:
    print("Python 3 is required, but you are using Python %i.%i.%i") % (
        sys.version_info[0], sys.version_info[1], sys.version_info[2])
    sys.exit(1)

def indel_len(sequence, start):
    l = 0
    while sequence[start + l] == '-':
        l += 1
    return l


def find_DNA_diff(read, ref, verbose):
    """
    @ read, ref: MutableSeq objects
    :return errors - tuple (position, expected triplet, actual triplet, ) / none if broken read

    Letter by letter report mutations in NGS read, all counts 1- based in result (code in 0-count).
    - substitution: 78C = nt 78 in reference is changed to C
    - deletions: 78d6 = 6 nt deleted after 78: 1-78, d6, 85-end
    - insertion: 78iATC = after nt 78 inserted seq. ATC
    """

    if read is None:
        if verbose:
            print('no read provided')
        return

    # No gap realignment at this point

    # quality control that there are no mutations at ends of reads
    ends = findEnds(read, ref)
    if not endMatch(read, ref, ends):
        if verbose:
            print('ends do not match')
        return

    # scan read & reference letter by letter, counting position in reference
    # reads have been trimmed so that reference starts @ 0
    dna_errors = []
    ref_index = ends.get('start')
    i = ends.get('start')

    while i < ends.get('end'):
        # check for differences
        if read[i] == ref[i]:
            ref_index += 1
            i += 1

        elif read[i] == '-':
            # start of a deletion
            l = indel_len(read, i)
            # now we know the length of a deletion, check for frameshifts
            if l % 3 == 0:
                dna_errors += [(str(ref_index), 'd', str(l))]  # deletion length l starting at ref_index in 0-count
                i += l
                ref_index += l
            else:
                dna_errors += [(str(ref_index), 'f')]
                break

        elif ref[i] == '-':
            # start of an insertion
            l = indel_len(ref, i)
            # check for frameshifts
            if l % 3 == 0:
                dna_errors += [(str(ref_index), 'i', str(read[i:i+l]) )]
                i += l
            else:
                dna_errors += [(str(ref_index), 'f')]
                break

        else:
            # substitution
            dna_errors += [(str(ref_index + 1), 's', str(read[i]) )]
            i += 1
            ref_index += 1

    return tuple(dna_errors)


def find_protein_diff(read, ref, verbose):

    # quality control
    if read is None:
        return
    ends = findEnds(read, ref)
    if not endMatch(read, ref, ends):
        return

    newread = read
    newref = ref

    # scan reference triplet by triplet
    # move letters when encountering an indel
    prot_errors = []
    i = ends.get('aligned')
    ref_index = int(ends.get('aligned') / 3) + 1  # reference amino acid index

    while i <= ends.get('end'):
        if newread is None:
            break
        ref_codon = newref[i:i+3]
        read_codon = newread[i:i+3]

        if '-' in read_codon:  # found a deletion
            # Check if this is the last acid, and it's incomplete, ignore it.
            if re.search('[ATGC]', str(newread[i + 3:])) is None:
                break

            if '-' in ref_codon:  # something very broken
                prot_errors.append((ref_index,'f'))
                return tuple(prot_errors)
            if read_codon == '---':  # single codon deletion
                prot_errors += [(ref_index, 'd')]
                i += 3
                ref_index += 1

            else:  # check it's not a frame shift
                l = indel_len(newread, i)
                if l % 3 != 0:
                    prot_errors.append((ref_index, 'f'))
                    return tuple(prot_errors)
                # realign gap and repeat loop at same position to compare the codons
                gap = findGap(newread[i - 1:])
                gap = (gap[0] + i - 1, gap[1] + i - 1)
                newread = gapAlign(newread, gap)
                continue

        elif '-' in ref_codon:  # found an insertion
            l = indel_len(newref, i)
            if l % 3 != 0:
                prot_errors.append((ref_index, 'f'))
                return tuple(prot_errors)
            gap = findGap(newref[i-1:])
            if gap[0] == 1:  # insertion after codon
                insertion = newread[gap[0] + i - 1:gap[1] + i - 1]
                if '-' in insertion:
                    prot_errors.append((ref_index, 'f'))
                    return tuple(prot_errors)
                prot_errors.append((ref_index, 'i', str(translate(insertion)) ))
                i += l
                ref_index += 1
            else:  # realign gap and repeat loop at same position to compare the codons
                gap = (gap[0] + i - 1, gap[1] + i - 1)
                newref = gapAlign(newref, gap)
                continue

        elif translate(read_codon) != translate(ref_codon):  # must be a substitution
            prot_errors.append((ref_index, 's', str(translate(read_codon))))
            i += 3
            ref_index += 1

        else:
            i += 3
            ref_index += 1

    if verbose:
        print(prot_errors)

    return tuple(prot_errors)


def depth_by_nt_position(depth_1, depth_2):
    """
    Collect samtools depth output into a list. 2nd column = 1-based position, 3rd column = coverage.
    Samtools gives two separate files for assembled and unassembled reads
    :param depth_1: output of samtools depth, tab delimited
    :param depth_2: same for other set of reads
    :return: a list of ints with coverage per position
    """

    nt_depth = []
    # depth of assembled and unassembled reads is in two separate files
    with open(depth_1, 'r') as f:
        for line in f.readlines():
            l = line.split()
            # depth in third column, 'samtools depth' output
            nt_depth.append(int(l[2]))
    # open second file and add the count to same position in depth list
    with open(depth_2, 'r') as f:
        for line in f.readlines():
            l = line.split()
            i = int(l[1]) - 1
            nt_depth[i] += int(l[2])
    return nt_depth


def depth_by_aa_position(nt_depth):
    """
    Convert nucleotide depth into average depth at amino acid position
    :param nt_depth:
    :return:
    """
    i = 0
    aa_depth = {}
    while i < len(nt_depth):
        pos = 1 + (i / 3)
        d = sum(nt_depth[i:i+3]) / 3
        aa_depth[pos] = d
        i += 3
    return aa_depth


def aa_depth_for_mutation(prot_error, aa_depth):

    if prot_error:
        e_depth = [aa_depth[point[0]] for point in prot_error]
    else:
        return

    avg_depth = sum(e_depth) / len(e_depth)

    return avg_depth


def count_one_fraction(alignment, aa_depth, debug):
    """
    Don't bother with expected/allowed mutations, just find everything and filter later
    Final format: {DNA error: [(protein error), fraction,
    1. Read reference file
    2. Scan over reference sequence to generate all possible mutations
    3. For each ref & read in multiple alignment:
        - verify the read is good quality
        - call the mutation
        - add to count table
    4. Print counts
    """
    # use a regular dictionary
    # when a protein mutation is first encountered, create an entry including depth
    one_lane_counts = {}

    # reading & looping over read/reference sequence in multiple sequence alignment
    # use AlignIO parser and keep sequence only, allowing it to change (important for gap shifts) 
    for pair in AlignIO.parse(alignment, "fasta", alphabet=IUPAC.ambiguous_dna, seq_count=2):
        # both read and ref are MutableSeq
        ref = pair[0].seq.tomutable()
        read = pair[1].seq.tomutable()
        read = MutableSeq(str(read).replace('N', '.'), read.alphabet)
        readname = pair[1].id

        # trim sequencing read to reference
        ref, read = trim_read(ref, read)

        try:
            dna_errors = find_DNA_diff(read, ref, debug)  # errors = a tuple
            prot_errors = find_protein_diff(read, ref, debug)
        except:
            if not dna_errors:
                print(dna_errors)
            print_coloured_diff(readname, read, ref, debug)
            raise

        try:
            one_lane_counts[prot_errors]['total'] += 1
            one_lane_counts[prot_errors]['dna'][dna_errors] += 1
        except KeyError:
            one_lane_counts[prot_errors] = {'dna': defaultdict(int),
                                            'depth': aa_depth_for_mutation(prot_errors, aa_depth),
                                            'total': 1}
            one_lane_counts[prot_errors]['dna'][dna_errors] += 1

    # count the mutations
    n = 0
    threshold = 10
    for error in one_lane_counts.keys():
        if one_lane_counts[error]['total'] > threshold:
            n += 1

    print('Fount {0} total protein mutations, of which {1} have more than {2} counts'
          .format(len(one_lane_counts), n, threshold))

    return one_lane_counts


def count_multiple_fractions(folder, baseline, debug):
    """
    Process all reference.fraction.aln files in given  folder in combination with
    reference.fraction.assembled.depth.txt and reference.fraction.unassembled.depth.txt
    :param folder: contains all *.aln and depth files
    :param baseline: string containing name of baseline
    :return:
    """
    all_references = {}

    print(os.listdir(folder))
    for f in os.listdir(folder):
        if f.endswith('.aln'):
            aln_path = os.path.join(folder, f)
            ref, fraction, suffix = f.rsplit(".", 2)
            print('Counting alignment {0} in background {1} and activity fraction {2}'
                  .format(f, ref, fraction))
            # prepare sequencing coverage / depth
            assembled = os.path.join(folder, ref + '.' + fraction + '.assembled.depth.txt')
            unassembled = os.path.join(folder, ref + '.' + fraction + '.unassembled.depth.txt')
            aa_depth = depth_by_aa_position(depth_by_nt_position(assembled, unassembled))

            if ref not in all_references.keys():
                all_references[ref] = {}

            if fraction == baseline:
                fraction = 'baseline'
            all_references[ref][fraction] = count_one_fraction(aln_path, aa_depth, debug)

    return all_references


def calculate_enrichements(all_references):

    enrichements = {}
    for ref in all_references.keys():
        if 'baseline' not in all_references[ref].keys():
            return all_references
        enrichements[ref] = {}
        for fraction in all_references[ref].keys():
            if fraction == 'baseline':
                continue
            enrichements[ref][fraction] = {}
            for e in all_references[ref][fraction].keys():
                try:
                    baseline_count = all_references[ref]['baseline'][e]['total']
                    baseline_depth = all_references[ref]['baseline'][e]['depth']
                    error_count = all_references[ref][fraction][e]['total']
                    error_depth = all_references[ref][fraction][e]['depth']
                except KeyError:
                    continue

                try:
                    baseline_freq = baseline_count / baseline_depth
                    error_freq = error_count / error_depth
                    # all_references[ref][fraction][e]['enriched'] = error_freq / baseline_freq
                    enrichements[ref][fraction][e] = error_freq / baseline_freq
                except TypeError:
                    enrichements[ref][fraction][e] = 0

    return enrichements


def combine_totals_same_reference(counts):
    """
    Start with a dictionary
    {prot_errors: {'H': 17, 'N': 12, ...}
    :param counts: dict containing {fraction: one_lane_counts} pairings
    :return:
    """

    one_reference_counts = defaultdict(partial(defaultdict, int))

    for fraction, one_lane_counts in counts.items():
        for prot_errors in one_lane_counts.keys():
            one_reference_counts[prot_errors][fraction] = one_lane_counts[prot_errors]['total']

    # df_protein = pd.DataFrame.from_dict(one_reference_counts, orient='index')

    return one_reference_counts


def find_shared_entries(counts1, counts2, cutoff=20):
    """
    Find mutations observed in both backgrounds
    :param counts1: {mutation: {H: 13, N: 0, M: 88}}
    :param counts2:
    :return:
    """
    c1 = {}
    c2 = {}
    for entry, counts in counts1.items():
        if max(counts.values()) > cutoff:
            try:
                if max(counts2[entry].values()) > cutoff:
                    c2[entry] = counts2[entry]
                    c1[entry] = counts1[entry]
            except ValueError:
                continue
    return c1, c2


def process_sanger_plate(summary, experimental, rc=False, q_cutoff=40, verbose=False):

    with open(summary) as f:
        reader = csv.DictReader(f, fieldnames=('ab1', 'activity', 'use', 'protein'))
        next(reader)
        for row in reader:
            if row['use'] == 'no':
                continue
            # read ab1, convert to trimmed fastq
            fqname = ab1_to_fastq(row['ab1'], rc=rc, q_cutoff=q_cutoff)
            try:
                ref, read, id = needle_align(fqname, args.reference)
                dna_errors = find_DNA_diff(read, ref, verbose)  # errors = a tuple
                protein_errors = find_protein_diff(read, ref, verbose)
                if verbose:
                    print(row['ab1'], dna_errors, protein_errors, row['protein'])
                if dna_errors:
                    try:
                        experimental[protein_errors]['dna'][dna_errors] = float(row['activity'])
                    except KeyError:
                        experimental[protein_errors] = {'dna': {}}
                        experimental[protein_errors]['dna'][dna_errors] = float(row['activity'])
            except ValueError:
                prefix, suffix = os.path.splitext(fqname)
                outname = prefix + '.aln'
                if os.stat(outname).st_size == 0:
                    print(outname, 'is empty.')
                    continue

    return experimental


def import_sanger(sanger_folder, verbose=False):
    """
    Expect all sequencing files to be in a folder
    :param sanger_folder:
    :return:
    """

    # set up experimental dictionary in same format as one lane counts
    # experimental[protein_mutation][dna_mutation] = activity
    # experimental[protei_mutation]['average'] = activity
    experimental = {}


    # add all folders - one folder is one sequencing plate
    cwd = os.getcwd()

    for d in os.listdir(sanger_folder):
        plate = cwd + os.sep + sanger_folder + os.sep + d
        print(plate)
        summary = ''
        for file in os.listdir(plate):
            if file.endswith(d + '.csv'):
                summary = d + '.csv'
                break
        os.chdir(plate)
        if summary == '':
            print('Cannot find activity CSV in {0} folder.'.format(plate))
        else:
            process_sanger_plate(summary, experimental, verbose=verbose, rc=True)

    # average synonymous mutations
    for p in experimental.keys():
        experimental[p]['average'] = sum(experimental[p]['dna'].values()) / len(experimental[p]['dna'])

    os.chdir(cwd)
    return experimental


def prepare_exp_for_pca(experimental, all_references, position_counts, style):
    """
    Prepare np arrays for all mutants for which activity is known.
    :param experimental:
    :param total:
    :return:
    """

    # get everything into one Data Frame
    seq_data = []  # enrichements in NGS
    act_data = []  # experimentally measured activity, a number
    index = []  # protein mutations
    columns = ['H', 'MM', 'L']

    n = 0  # counter for mutations with poor NGS data
    # experimental[protein_mutation]['average'] = activity

    for e in experimental.keys():
        enrich = enrichements['eGFP'][e]
        # discard mutations for which there is no sequencing data
        if all_references['eGFP']['baseline'][e]['total'] < 10:
            n += 1
            continue
        seq = {c: enrich[c] for c in columns}
        act = {'activity': experimental[e]['average']}

        if act['activity'] >= 0.6:
            act['fraction'] = 'H'
        elif act['activity'] > 0:
            act['fraction'] = 'MM'
        elif act['activity'] == 0:
            act['fraction'] = 'L'

        seq_data.append(seq)
        act_data.append(act)
        index.append(e)

    exp_seq = pd.DataFrame(seq_data, index=index)
    exp_act = pd.DataFrame(act_data, index=index)

    act = {'H': 2, 'MM': 1, 'L': 0}
    exp_act['fraction_num'] = exp_act.apply(lambda x: act[x['fraction']], axis=1)

    print('No enrichement:', n)
    return exp_seq, exp_act


def pca_experimental(exp_seq, exp_act, graph=False):
    """
    PCA on sequencing enrichment.
    :param exp_seq: data frame with sequencing enrichment indexed by error name
    :param exp_act: data frame with activity data & classification indexed by error name
    :return:
    """
    # keep original data safe
    exp_std = exp_seq.copy()
    # scale data to unit variance on each component
    exp_std[exp_std.columns] = StandardScaler().fit_transform(exp_std[exp_std.columns])
    pca = decomposition.PCA(3)

    exp_pca = pd.DataFrame(pca.fit_transform(exp_std), columns=['PCA%i' % i for i in range(1, 4)], index=exp_seq.index)

    print('Explained variance by PCA', pca.explained_variance_ratio_)
    print('Using {} experimentally measured mutations'.format(len(exp_pca)))

    if graph:
        # make plots
        max_act = exp_act.activity.max()
        exp_act['col'] = exp_act['activity'] / max_act
        exp_act['col'] = exp_act['col'].astype(str)

        fig, ax = plt.subplots()
        cax = ax.scatter(exp_pca.PCA1, exp_pca.PCA2, c=exp_act.col, cmap='Greens', marker='o', s=20, edgecolors='grey',
                   norm=colors.SymLogNorm(linthresh=0.1, linscale=1, vmin=0, vmax=max_act))
        ax.set_title('Log fluorescence')
        ax.set_xlabel('n PCA 1')
        ax.set_ylabel('n PCA 2')

        cbar = fig.colorbar(cax)
        cbar.ax.set_yticklabels(['0'] + ['']*8 + ['{0:.1f}'.format(0.1*k) for k in range(1,11)])

    return pca, exp_pca


def pca_kmeans(exp_pca, nclusters=2):

    kmeans = KMeans(n_clusters=nclusters, n_init=50).fit(exp_pca)

    plt.figure()
    for i in range(77):
        if kmeans.labels_[i] == 0:
            plt.scatter(exp_pca.ix[i, 'PCA1'], exp_pca.ix[i, 'PCA2'], c='r', marker='o')
        elif kmeans.labels_[i] == 1:
            plt.scatter(exp_pca.ix[i, 'PCA1'], exp_pca.ix[i, 'PCA2'], c='g', marker='o')
        elif kmeans.labels_[i] == 2:
            plt.scatter(exp_pca.ix[i, 'PCA1'], exp_pca.ix[i, 'PCA2'], c='b', marker='o')
    for i in range(0, nclusters):
        plt.scatter(kmeans.cluster_centers_[i][0], kmeans.cluster_centers_[i][1], marker='+')

    return kmeans

if __name__ == "__main__":
    """
    Arguments:
    1.   Alignment of sequencing reads to a reference sequence in FASTA format. The
         reference comes before read sequence. (>Ref, seq, >Read, seq).
    (2.) Optional: DEBUG print a coloured representation of mismatches
    """
    parser = argparse.ArgumentParser(description='Finds all in-frame mutations in a gene')
    parser.add_argument('-f', '--folder', help='Folder containing multiple sequence alignments', required=True)
    parser.add_argument('-r', '--reference', required=False)
    parser.add_argument('-d', '--debug', help='Turn on debugging', required=False, action="store_true")  # Visual
    parser.add_argument('-b', '--baseline', help='Name of baseline fraction', required=False)
    args = parser.parse_args()


    """
    PROCESS NGS DATA
    """

    # all_references = count_multiple_fractions(args.folder, args.baseline, args.debug)
    #
    # with open('everything.p', 'wb') as f:
    #     pickle.dump(all_references, f)

    with open('everything.p', 'rb') as f:
        all_references = pickle.load(f)

    enrichements = calculate_enrichements(all_references)

    # """
    # PROCESS SANGER SEQUENCING DATA
    # """
    # exp = import_sanger(os.path.join(args.folder, 'sanger'), verbose=True)

    #
    # eGFP = combine_totals_same_reference(all_references['eGFP'])
    # GFP8 = combine_totals_same_reference(all_references['GFP8'])
    #
    # for e in experimental.keys():
    #     print('Error: {} \t Activity: {} \n eGFP:{} \n GFP8:{} \n'.format(e, experimental[e], eGFP[e], GFP8[e]))
    #
    # # eGFP_shared, GFP8_shared = find_shared_entries(eGFP, GFP8)
    # #
    # # df_eGFP = pd.DataFrame.from_dict(eGFP_shared, orient='index')
    # # df_GFP8 = pd.DataFrame.from_dict(GFP8_shared, orient='index')
    #
    #
    #
    # # # both = pd.merge(df_eGFP, df_GFP8, how='inner', sort=False)
    #
