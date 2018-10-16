#!/usr/bin/python3

import sys
import pickle
import re
import os
import csv
import argparse
import pprint
import random
import time

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from sklearn import decomposition
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

from collections import defaultdict

from Bio import AlignIO, SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq, translate

from indels.ind import trim_read, findEnds, endMatch, findGap, gapAlign, ab1_to_fastq, needle_align
from indels.output import print_coloured_diff

# Demand Python 3.
if sys.version_info[0] < 3:
    print("Python 3 is required, but you are using Python %i.%i.%i") % (
        sys.version_info[0], sys.version_info[1], sys.version_info[2])
    sys.exit(1)

# Identifying mutations from fasta alignment


def indel_len(sequence, start):
    l = 0
    while sequence[start + l] == '-':
        l += 1
    return l


def find_DNA_hgvs(read, ref, refname, verbose=False, start_offset=3, end_trail=3):
    """@ read, ref: MutableSeq objects
    :return errors - tuple (position, expected triplet, actual triplet, ) / none if broken read

    The assumption is that the reference includes an offset of 3 nt either side of the gene of interest. The starting triplet is
    reported as 'amino acid 0'. If the offset is less or more, it needs to be set explicitly. end_trail specifies the
    number of nt after end of gene and is ignored
    """
    if read is None:
        if verbose:
            print('no read provided')
        return

    # No gap realignment at this point
    prefix = str(refname) + ':c.'

    # quality control that there are no mutations at ends of reads
    ends = findEnds(read, ref, start_offset)
    if not endMatch(read, ref, ends):
        if verbose:
            print('ends do not match')
        return

    # scan read & reference letter by letter, counting position in reference
    # reads have been trimmed so that reference starts @ 3 (0,1,2 is the extra triplet)
    # in the general case, reference starts @ offset in 0-count
    # This is equal to the number of nt before ATG
    # ref_index denotes HGVS DNA position labeling, i is used for accessing sequence
    dna_errors = []
    ref_index = ends.get('start') - start_offset + 1  # if the read starts at 3, this becomes nt 1 (1-based as is HGVS)
    i = ends.get('start')

    if ref_index < 1: # don't wish to consider mutations in the vector sequence before start of gene
        skip = 1 - ref_index
        i += skip
        ref_index += skip

    while i < ends.get('end') - end_trail: # the trailing nt are ignored when reading mutations
        # check for differences
        if read[i] == ref[i]:
            ref_index += 1
            i += 1

        elif read[i] == '-':
            # start of a deletion, format depends on length
            l = indel_len(read, i)
            if l == 1:  # format is POSdel
                dna_errors.append(str(ref_index) + 'del')
            else:
                # format is FIRST_LASTdel
                dna_errors.append(str(ref_index) + '_' + str(ref_index + l - 1) + 'del')
            i += l
            ref_index += l

        elif ref[i] == '-':
            # start of an insertion, format is FLANK_FLANKinsSEQ
            l = indel_len(ref, i)
            dna_errors.append(str(ref_index -1) + '_' + str(ref_index) + 'ins' + str(read[i:i+l]) )
            i += l

        else:
            # substitution: need to include ref. sequence in format 8A>G
            dna_errors.append(str(ref_index) + str(ref[i]) + '>' + str(read[i]))
            i += 1
            ref_index += 1

    # format the result including name of sequence
    if len(dna_errors) == 1:
        dna_hgvs = prefix + dna_errors[0]
    else:
        dna_hgvs = prefix + '[' + (';').join(dna_errors) + ']'

    return dna_hgvs


def find_DNA_diff(read, ref, verbose=False, start_offset=3, end_trail=3):
    """
    @ read, ref: MutableSeq objects
    :return errors - tuple (position, expected triplet, actual triplet, ) / none if broken read

    The assumption is that the reference includes 3 nt either side of the gene of interest. The starting triplet is
    reported as 'amino acid 0'.
    As for HGVS, the starting offset and number of trailing nt are variable
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
    ends = findEnds(read, ref, start_offset)
    if not endMatch(read, ref, ends):
        if verbose:
            print('ends do not match')
        return

    # scan read & reference letter by letter, counting position in reference
    # reads have been trimmed so that reference starts @ offset=3 by default (0,1,2 is the extra triplet)
    dna_errors = []
    ref_index = ends.get('start') - start_offset + 1
    i = ends.get('start')

    if ref_index < 1: # don't wish to consider mutations in the vector sequence before start of gene
        skip = 1 - ref_index
        i += skip
        ref_index += skip

    while i < ends.get('end') - end_trail:
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


def find_protein_diff(read, ref, verbose=False, start_offset=3, end_trail=3):

    # quality control
    if read is None:
        return
    ends = findEnds(read, ref, start_offset)
    if not endMatch(read, ref, ends):
        return

    newread = read
    newref = ref

    # scan reference triplet by triplet
    # move letters when encountering an indel
    prot_errors = []
    i = ends.get('aligned')
    ref_index = int((ends.get('aligned') - start_offset )/ 3) + 1  # reference amino acid index

    if ref_index < 1: # don't wish to consider mutations in the vector sequence before start of gene
        skip = 1 - ref_index
        i += 3*skip
        ref_index += skip

    while i <= ends.get('end') - end_trail:
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
            if str(translate(read_codon)) == '*':
                return tuple(prot_errors)
            i += 3
            ref_index += 1

        else:
            i += 3
            ref_index += 1

    if verbose:
        print(prot_errors)

    return tuple(prot_errors)


# Add sequencing depth and convert to enrichments

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


# Raw processing of all alignments, get composition

def count_one_fraction(alignment, aa_depth, refname, debug):
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
        dna_errors, dna_hgvs, prot_errors = None, None, None

        try:
            dna_errors = find_DNA_diff(read, ref, debug)  # errors = a tuple
            dna_hgvs = find_DNA_hgvs(read, ref, refname, debug)  # string according to HGVS format (ish)
            prot_errors = find_protein_diff(read, ref, debug)
        except:
            if not dna_errors:
                print(dna_errors)
            print_coloured_diff(readname, read, ref, debug)
            raise

        try:
            one_lane_counts[prot_errors]['total'] += 1
            one_lane_counts[prot_errors]['dna'][dna_errors] += 1
            one_lane_counts[prot_errors]['dna_hgvs'][dna_hgvs] += 1
        except KeyError:
            one_lane_counts[prot_errors] = {'dna': defaultdict(int), 'dna_hgvs': defaultdict(int),
                                            'depth': aa_depth_for_mutation(prot_errors, aa_depth),
                                            'total': 1}
            one_lane_counts[prot_errors]['dna'][dna_errors] += 1
            one_lane_counts[prot_errors]['dna_hgvs'][dna_hgvs] += 1


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
            refname, fraction, suffix = f.rsplit(".", 2)
            print('Counting alignment {0} in background {1} and activity fraction {2}'
                  .format(f, refname, fraction))
            # prepare sequencing coverage / depth
            assembled = os.path.join(folder, refname + '.' + fraction + '.assembled.depth.txt')
            unassembled = os.path.join(folder, refname + '.' + fraction + '.unassembled.depth.txt')
            # assembled = os.path.join(folder, ref + '.assembled.depth.txt')
            # unassembled = os.path.join(folder, ref + '.unassembled.depth.txt')
            aa_depth = depth_by_aa_position(depth_by_nt_position(assembled, unassembled))

            if refname not in all_references.keys():
                all_references[refname] = {}

            if fraction == baseline:
                fraction = 'baseline'
            all_references[refname][fraction] = count_one_fraction(aln_path, aa_depth, refname, debug)

    return all_references


def classify_dna(dna_error):
    if dna_error is None:  # empty or broken reads
        return 'b'
    elif len(dna_error) > 1:
        # expect substitutions
        if dna_error[-1][1] == 'f':  # frameshifts are always the last mutation
            return 'f'
        else:
            for k in range(len(dna_error)):
                if dna_error[k][1] == 'i':
                    return 'si'
                elif dna_error[k][1] == 'd':
                    return 'sd'
            return 's'
    elif len(dna_error) == 1:
        if dna_error[0][1] == 'f':
            return 'f'
        elif dna_error[0][1] == 'd':
            return 'd' + dna_error[0][2]  # length of deletion
        elif dna_error[0][1] == 'i':
            return 'i' + str(len(dna_error[0][2]))
        elif dna_error[0][1] == 's':
            return 's'
    else:
        return 'b'


def get_dna_composition(all_references, cutoff=10):
    dna_count = {}
    dna_reads = {}
    for background in all_references.keys():
        for fraction in all_references[background].keys():
            print('Analysing background {0} and fraction {1}. Unusal mutations: '.format(background, fraction))
            distinct_mutations = 0
            total_count = 0
            dna_count[background + fraction] = {'s': 0, 'd3': 0, 'd6': 0, 'd9': 0, 'i3': 0, 'i6': 0, 'i9': 0, 'f': 0,
                                                'b': 0, 'other': 0, 'sd': 0, 'si': 0}
            dna_reads[background + fraction] = {'s': 0, 'd3': 0, 'd6': 0, 'd9': 0, 'i3': 0, 'i6': 0, 'i9': 0, 'f': 0,
                                                'b': 0, 'other': 0, 'sd': 0, 'si': 0}
            for mutation in all_references[background][fraction].keys():
                mut_total = all_references[background][fraction][mutation]['total']
                if mut_total >= cutoff:
                    # find all DNA entries with high enough counts
                    for dna_error, c in all_references[background][fraction][mutation]['dna'].items():
                        # dna_error = max(all_references[background][fraction][mutation]['dna'],
                        #  key=lambda key: all_references[background][fraction][mutation]['dna'][key])
                        if c >= cutoff:
                            try:
                                dna_type = classify_dna(dna_error)
                                dna_count[background + fraction][dna_type] += 1
                                dna_reads[background + fraction][dna_type] += c
                                distinct_mutations += 1
                                total_count += c
                            except KeyError:
                                # print(dna_type, dna_error, mut_total)
                                dna_count[background + fraction]['other'] += 1
                                dna_reads[background + fraction]['other'] += c
            print('In background {0} and fraction {1} found {2} distinct mutations with total read count {3}'.format(
                background, fraction, distinct_mutations, total_count
            ))

    return pd.DataFrame.from_dict(dna_count), pd.DataFrame.from_dict(dna_reads)


def is_mutation_consecutive(mutation):
    """
    If only consecutve amino acids are affected, return 'c', else return 'nc'
    :param mutation:
    :return:
    """
    for pos in range(1, len(mutation)):
        if mutation[pos][0] != (mutation[pos - 1][0] + 1) :
            return 'nc'
    return 'c'


def classify_protein(mutation):
    if mutation is None: # came from empty or broken reads
        return 'b'
    else:
        m = []
        for pos in range(len(mutation)):
            t = mutation[pos][1]
            if t != 'i':
                m.append(t)
            else:
                m.append(t + str(len(mutation[pos][2])) )
        # need to distinguish between consecutive mutations and likely sequencing errors
        if 'f' in m:
            return 'f'
        elif len(m) == 1:
            return ''.join(m)
        elif len(m) >= 1:
            c = is_mutation_consecutive(mutation) + '-' + ''.join(m)
            return c


def get_protein_composition(all_references, cutoff=10):
    protein_count = {}
    protein_reads = {}
    for background in all_references.keys():
        for fraction in all_references[background].keys():
            print('Analysing protein composition in background {0} and fraction {1}. Unusal mutations: '.format(background, fraction))
            distinct_mutations = 0
            total_count = 0
            protein_count[background + '.' + fraction] = {'s': 0, 'd': 0, 'c-dd': 0, 'c-ddd': 0, 'i1': 0, 'i2': 0,
                            'i3': 0, 'c-ss': 0, 'c-sd': 0, 'c-sdd': 0, 'c-sddd': 0, 'c-si1': 0, 'c-si2': 0, 'c-si3': 0,
                            'f': 0, 'other': 0, 'b': 0}
            protein_reads[background + '.' + fraction] = {'s': 0, 'd': 0, 'c-dd': 0, 'c-ddd': 0, 'i1': 0, 'i2': 0,
                            'i3': 0, 'c-ss': 0, 'c-sd': 0, 'c-sdd': 0, 'c-sddd': 0, 'c-si1': 0, 'c-si2': 0, 'c-si3': 0,
                            'f': 0, 'other': 0, 'b': 0}
            for mutation in all_references[background][fraction].keys():
                mut_total = all_references[background][fraction][mutation]['total']
                if mut_total >= cutoff:
                    # find all DNA entries with high enough counts
                    try:
                        prot_type = classify_protein(mutation)
                        protein_count[background + '.' + fraction][prot_type] += 1
                        protein_reads[background + '.' + fraction][prot_type] += mut_total
                        distinct_mutations += 1
                        total_count += mut_total
                    except KeyError:
                        # print(prot_type, mutation, mut_total)
                        protein_count[background + '.' + fraction]['other'] += 1
                        protein_reads[background + '.' + fraction]['other'] += mut_total
            print('In background {0} and fraction {1} found {2} distinct mutations with total read count {3}'.format(
                background, fraction, distinct_mutations, total_count
            ))

    return pd.DataFrame.from_dict(protein_count), pd.DataFrame.from_dict(protein_reads)


def what_appears_n_times(all_references, n, style='protein', fractions=None, effect=None):
    for background in all_references.keys():
        if fractions is None:
            fractions = all_references[background].keys()
        for fraction in fractions:
            print('In {0} background and fraction {1} these mutations appear {2} times:'.format(background, fraction, n))
            for mutation in all_references[background][fraction].keys():
                if style == 'protein':
                    if all_references[background][fraction][mutation]['total'] == n\
                            and classify_protein(mutation) in effect:
                        print(pretty_mutation(mutation))
                elif style == 'dna':
                    for dna, c in all_references[background][fraction][mutation]['dna'].items():
                        if c == n and classify_dna(dna) in effect:
                            print(pretty_mutation(dna))
                else:
                    print('Please specify either dna or protein as style.')
                    return


def insertion_composition(all_references, cutoff=2, l=(3,6,9)):

    comp = {length: {k: {'A':0, 'C':0, 'T':0, 'G':0} for k in range(length)} for length in l}
    for background in all_references.keys():
        for fraction in all_references[background].keys():
            for mutation in all_references[background][fraction].keys():
                for dna_error, c in all_references[background][fraction][mutation]['dna'].items():
                    if dna_error is None:
                        continue
                    if c >= cutoff and len(dna_error) == 1 and classify_dna(dna_error) in ('i3', 'i6', 'i9'):
                        ins = dna_error[0][2]
                        ins_len = len(ins)
                        if ins_len in l:
                            for pos in range(ins_len):
                                comp[ins_len][pos][ins[pos]] += 1
    return comp


def find_transposon_histogram(all_references, background, baseline='baseline', transposon='d3'):
    """
    Find all mutation of a certain type and count where they are in DNA
    :return: dict
    """
    hist = defaultdict(int)

    for prot_mutation in all_references[background][baseline]:
        if prot_mutation is None:
            continue
        elif len(prot_mutation) <= 2:
            for dna_mutation, count in all_references[background][baseline][prot_mutation]['dna'].items():
                if classify_dna(dna_mutation) == transposon:
                    hist[int(dna_mutation[0][0])] += count

    return hist


def rev_comp(nt):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return complement[nt]

def transposon_consensus_seq(all_references, reference, fraction='d3', transposon='d3'):
    """
    Determine the consensus sequence for transposon insertion, by analysing the position of d3 or i3 mutations
    :param all_references: dictionary containing all mutation data
    :param reference: BioSeq fasta reference
    :param fraction: name of the activity fraction / library to be analysed
    :param transposon: which type of mutations are we counting
    :return: dict with composition by position
    """
    consensus = {pos: {'A': 0, 'T': 0, 'C': 0, 'G': 0} for pos in range(5)}
    baseline = {pos: {'A': 0, 'T': 0, 'C': 0, 'G': 0} for pos in range(5)}
    ref = SeqIO.read(reference, 'fasta')

    for i in range(2, len(ref)-7):
        trans_seq = str(ref[i:i+5].seq)
        for pos in range(5):
            baseline[pos][trans_seq[pos]] += 1
            baseline[pos][rev_comp(trans_seq[4 - pos])] += 1

    background = str(ref.name)

    deletions = all_references[background][fraction]
    for prot_mutation in deletions.keys():
        for dna_mutation, count in deletions[prot_mutation]['dna'].items():
            if classify_dna(dna_mutation) == transposon:  # found the simple mutations: d3 or i3
                start = int(dna_mutation[0][0])
                # this position refers to the nt BEFORE the indel in 1-count: hence this is the 1st nt of 5 nt site
                # say the reference is nnn ATG CTG AAC: for start = 1, we want to retrieve ATGCT -> ref[3:8]
                if transposon == 'd3':
                    trans_seq = str(ref[start+2:start+7].seq)
                # for insertions we want 4 nt before and 1 after the insertion
                elif transposon == 'i3':
                    trans_seq = str(ref[start - 1:start + 4].seq)
                if len(trans_seq) != 5:
                    continue
                for pos in range(5):
                    consensus[pos][trans_seq[pos]] += count
                    consensus[pos][rev_comp(trans_seq[4 - pos])] += count

    return baseline, consensus


def dna_mutation_frequencies(all_references):
    """
    Generate data for a histogram of DNA mutation frequencies - how many occur once, twice, ...
    :param all_references:
    :return:
    """
    freq = {}  # for each library ('d3', 'i6', etc.) give a dictionary {1: 5 times, 2: 3 times...}
    for background in all_references.keys():
        for fraction in all_references[background].keys():
            freq[background + fraction] = defaultdict(int)
            for prot_mutation in all_references[background][fraction].keys():
                for dna_mutation, c in all_references[background][fraction][prot_mutation]['dna'].items():
                    if classify_dna(dna_mutation) == fraction:
                        freq[background + fraction][c] += 1
    return freq


def insertion_frequencies(all_references):
    """
    Generate data for number of different insertions per position.
    :param all_references:
    :return:
    """
    ins_freq = {}  # for each library ('d3', 'i6', etc.) give a dictionary {1: 5 times, 2: 3 times...}
    for background in all_references.keys():
        for fraction in all_references[background].keys():
            ins_freq[fraction] = defaultdict(int)
            for prot_mutation in all_references[background][fraction].keys():
                for dna_mutation, c in all_references[background][fraction][prot_mutation]['dna'].items():
                    if classify_dna(dna_mutation) == fraction:
                        ins_freq[fraction][dna_mutation[0][0]] += 1
    return ins_freq


def export_hgvs(all_references, output):
    """
    Export DNA mutations and their counts in HGVS format, indended for use with Enrich2.
    :param all_references:
    :param output: prefix for CSV files
    :return:
    """
    for background in all_references.keys():
        for fraction in all_references[background].keys():
            with open('.'.join((output, background, fraction, 'csv')), 'w') as f:
                hgvs_writer = csv.writer(f, delimiter=',')
                for prot in all_ref[background][fraction].keys():
                    for hgvs, count in all_ref[background][fraction][prot]['dna_hgvs'].items():
                        if hgvs is None:
                            continue
                        hgvs_writer.writerow([hgvs, count])


# Find activity

def calculate_enrichments(all_references):

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


def combine_same_reference(dictionary, style='enrich'):
    """
    Start with a dictionary in the form reference[fraction][prot_errors]('total)=number and turn it into
    reference[prot_errors] = {'H': 0.8, 'L': 12, ...}
    :param dictionary: dict containing {fraction: number} pairings
    :return:
    """
    assert style in ['enrich', 'counts']

    combined = {ref: {} for ref in dictionary.keys()}

    # invert the order

    for ref in dictionary.keys():
        for fraction in dictionary[ref].keys():
            for prot_errors in dictionary[ref][fraction].keys():
                if style == 'enrich':
                    value = dictionary[ref][fraction][prot_errors]
                elif style == 'counts':
                    value = dictionary[ref][fraction][prot_errors]['total']
                try:
                    combined[ref][prot_errors][fraction] = value
                except KeyError:
                    combined[ref][prot_errors] = {fraction: value}

    return combined


def max_activity_fraction(enrich):
    """
    Find the activity fraction with highest enrichment in a dict of format {H: 0.1, M: 17, L:0}
    :param enrich:
    :return:
    """
    inverted = {value: frac for frac, value in enrich.items()}
    max_enrichment = max(inverted.keys())
    max_fraction = inverted[max_enrichment]
    return max_fraction


def find_shared_entries(counts1, counts2, cutoff=10):
    """
    Find mutations observed in both backgrounds
    :param counts1: {mutation: {H: 13, N: 0, M: 88}}
    :param counts2: {mutation: {H: 13, N: 0, M: 88}}
    :return: list of mutations with good counts in both
    """
    shared = []

    for prot_mutation, counts in counts1.items():
        if max(counts.values()) > cutoff:
            try:
                if max(counts2[prot_mutation].values()) > cutoff:
                    shared.append(prot_mutation)
            except KeyError:
                continue

    return shared


# Collect all counts and enrichments into a sequence similarity network

def are_mutations_similar(mutation1, mutation2):
    """
    Similar if they share at least one type & position of mutation, eg. 88d.
    :param mutation1:
    :param mutation2:
    :return:
    """
    positions1 = set(p[0:2] for p in mutation1)
    positions2 = set(p[0:2] for p in mutation2)

    if positions1 == positions2:
        return False

    max_length = max([len(positions1), len(positions2)])
    cutoff = max(1, max_length - 1)

    shared = positions1 & positions2

    if len(shared) >= cutoff:
        return True
    else:
        return False


def pretty_mutation(prot_error):
    if prot_error is None:
        return ''
    else:
        flat = []
        for point in prot_error:
            s = ''.join(str(e) for e in point)
            flat.append(s)
        return ' '.join(flat)


def generate_positions(ref):
    """
    Make a list of (int, s) for all positions in reference fasta
    :param ref: name of fasta file
    :return: list of tuples
    """
    reference = SeqIO.read(ref, "fasta")
    positions = [((i, 's'), (i, 'd')) for i in range(1, len(reference) + 1)]
    return positions


def generate_mutations_for_network(shared, combined_counts={}, combined_enrich={}):

    mutations = {e: {'type': 'NGS'} for e in shared if e}

    for e in shared:
        if e is None:
            continue

        for ref in combined_enrich.keys():
            try:
                mutations[e]['.'.join([str(ref), 'max_enrich'])] = max_activity_fraction(combined_enrich[ref][e])
            except KeyError:
                continue
            for frac in combined_counts[ref][e].keys():
                mutations[e]['.'.join([str(ref), str(frac), 'NGS_count'])] = combined_counts[ref][e][frac]

    return mutations


def make_similarity_network(mutations, ref):
    """
    Build a graph connecting mutations that
    :param mutations: dictionary with mutation tuples as keys and a dictionary of attributes as values
    :return: network
    """

    network = nx.Graph()

    # generate the protein backbone using generic 1s, 2s, 3s... mutations and connect them
    positions = generate_positions(ref)
    for e in positions:
        network.add_node(pretty_mutation(e), type='backbone', position=e[0][0])
    for i in range(1, len(positions)):
        network.add_edge(pretty_mutation(positions[i-1]), pretty_mutation(positions[i]), type='backbone')

    backbone = network.copy()

    # add in all mutations, their attributes and connect them to backbone
    for name, att in mutations.items():
        if name is None:
            continue
        network.add_node(pretty_mutation(name), **att)
        for e in positions:
            if are_mutations_similar(e, name):
                network.add_edge(pretty_mutation(e), pretty_mutation(name), type='NGS_to_backbone')

    # now add all other connections
    for mut1 in mutations.keys():
        if mut1 is None:
            continue
        for mut2 in mutations.keys():
            if mut2 is None:
                continue
            if are_mutations_similar(mut1, mut2):
                network.add_edge(pretty_mutation(mut1), pretty_mutation(mut2), type='NGS_to_NGS')

    return backbone, network


# Add activity values for Sanger sequenced clones - currently broken

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

    # """
    # PROCESS SANGER SEQUENCING DATA
    # """
    # exp = import_sanger(os.path.join(args.folder, 'sanger'), verbose=True)

    """
    PROCESS NGS DATA
    1. Count all mutations
    2. Calculate enrichments
    3. Find list of good quality mutations
    4. Call activity scores
    5. Make sequence similarity network
    """

    all_ref = count_multiple_fractions(args.folder, args.baseline, args.debug)
    export_hgvs(all_ref, 'new_script')

    # with open('GFP8_with_hgvs.p', 'wb') as f:
    #     pickle.dump(all_ref, f)

    # with open('gfp_with_hgvs_20180827.p', 'rb') as f:
    #     all_ref = pickle.load(f)
    #
    # d3_baseline, d3_cons = transposon_consensus_seq(all_ref, args.reference, fraction='baseline', transposon='d3')

    # with open('S6.p', 'rb') as f:
    #     all_ref = pickle.load(f)
    #
    # d3_baseline, d3_cons = transposon_consensus_seq(all_ref, args.reference, fraction='d3', transposon='d3')
    #
    # pprint.pprint(d3_baseline)
    # pprint.pprint(d3_cons)
    #
    # i3_cons = transposon_consesus_seq(all_ref, args.reference, fraction = 'i3', transposon = 'i3')

    # for cutoff in [1,2,5,10]:
    #     dna_count, dna_reads = get_dna_composition(all_ref, cutoff)
    #     protein_count, protein_reads = get_protein_composition(all_ref, cutoff)
    #
    #     print(cutoff)
    #     print(dna_count)
    #     print(dna_reads)
    #     print(protein_count)
    #     print(protein_reads)

    # what_appears_n_times(all_ref, 20, style='protein', effect=['d'], fractions=['d3'])

    # comp = insertion_composition(S6)
    #

    # hist_i3 = find_transposon_histogram(all_ref, 'S6', 'i3', transposon = 'i3')
    # print(hist_i3)
    # hist8 = find_transposon_histogram(all_ref, 'GFP8', 'baseline')
    # freq = dna_mutation_frequencies(all_ref)
    # ins_freq = insertion_frequencies(all_ref)

    # combined_enrich = combine_same_reference(calculate_enrichments(all_references))
    # combined_counts = combine_same_reference(all_references, style='counts')
    #
    # hist = find_transposon_histogram(all_references, 'eGFP')


    # shared = find_shared_entries(combined_counts['eGFP'], combined_counts['GFP8'])
    # mutations = generate_mutations_for_network(shared, combined_counts, combined_enrich)
    # backbone, net = make_similarity_network(mutations, args.reference)
    # # nx.write_graphml(backbone, 'GFP_v7_backbone.graphml')
    # nx.write_graphml(net, 'GFP_v8.graphml')

    #
    # for e in experimental.keys():
    #     print('Error: {} \t Activity: {} \n eGFP:{} \n GFP8:{} \n'.format(e, experimental[e], eGFP[e], GFP8[e]))
    #
    # #
    # #
    # # df_eGFP = pd.DataFrame.from_dict(eGFP_shared, orient='index')
    # # df_GFP8 = pd.DataFrame.from_dict(GFP8_shared, orient='index')
    #
    #
    #
    # # # both = pd.merge(df_eGFP, df_GFP8, how='inner', sort=False)
    #
