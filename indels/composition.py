#!/usr/bin/env python3

import sys
import pickle
import re
import os
import csv
import argparse
import pprint

import pandas as pd
from collections import defaultdict

from Bio import AlignIO, SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq, translate

from ind import trim_read, findEnds, endMatch, findGap, gapAlign
from output import print_coloured_diff, printErrors

# Demand Python 3.
if sys.version_info[0] < 3:
    print("Python 3 is required, but you are using Python %i.%i.%i") % (
        sys.version_info[0], sys.version_info[1], sys.version_info[2])
    sys.exit(1)

# Identifying mutations from fasta alignment


def read_is_wt(read, ref):
    """
    Check whether read sequence is exactly equal to wt.
    :param read:
    :param ref:
    :return:
    """
    trimmed_read = re.search(r'^-+([AGCTN][ACGTN-]+[ACGTN])-+$', str(read))
    if trimmed_read is None:
        return True
    else:
        return re.search(str(trimmed_read.group(1)), str(ref)) is not None


def indel_len(sequence, start):
    l = 0
    while sequence[start + l] == '-':
        l += 1
    return l


def find_DNA_hgvs(read, ref, refname, verbose=False, start_offset=3, end_trail=3):
    """@ read, ref: MutableSeq objects
    :return errors - tuple (position, expected triplet, actual triplet, ) / none if broken read

    The assumption is that the reference includes an offset of 3 nt either side of the gene of interest, such that the
    starting triplet is reported as 'amino acid 0'. If the offset is different, it needs to be set explicitly.
    end_trail specifies the  number of nt after end of gene and is ignored
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
    max_i = len(ref) - end_trail

    while i < ends.get('end'):
        if i > max_i:  # the trailing nt are ignored when reading mutations
            break
        # check for differences
        if read[i] == ref[i]:
            ref_index += 1
            i += 1

        elif read[i] == '-':
            # start of a deletion, format depends on length
            l = indel_len(read, i)
            if ref_index > 0:
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
            if ref_index > 0:
                dna_errors.append(str(ref_index - 1) + '_' + str(ref_index) + 'ins' + str(read[i:i+l]))
            i += l

        else:
            # substitution: need to include ref. sequence in format 8A>G
            if ref_index > 0:
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
    - deletions: 78d6 = 6 nt deleted starting with 78: 1-77, d6, 84-end
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
    max_i = len(ref) - end_trail

    while i < ends.get('end'):
        if i > max_i:
            break
        # check for differences
        if read[i] == ref[i]:
            ref_index += 1
            i += 1

        elif read[i] == '-':
            # start of a deletion
            l = indel_len(read, i)
            # now we know the length of a deletion, check for frameshifts
            if l % 3 == 0:
                if ref_index > 0:
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
                if ref_index > 0:
                    dna_errors += [(str(ref_index), 'i', str(read[i:i+l]))]
                i += l
            else:
                dna_errors += [(str(ref_index), 'f')]
                break

        else:
            # substitution
            if ref_index > 0:
                dna_errors += [(str(ref_index + 1), 's', str(read[i]))]
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
    ref_index = int((ends.get('aligned') - start_offset)/3) + 1  # reference amino acid index
    max_i = len(ref) - end_trail

    while i <= ends.get('end'):
        if i > max_i:
            break

        if newread is None:
            break
        ref_codon = newref[i:i+3]
        read_codon = newread[i:i+3]

        if '-' in read_codon:  # found a deletion
            # Check if this is the last acid, and it's incomplete, ignore it.
            if re.search('[ATGC]', str(newread[i + 3:])) is None:
                break

            if '-' in ref_codon:  # something very broken
                prot_errors.append((ref_index, 'f'))
                return tuple(prot_errors)
            elif read_codon == '---':  # single codon deletion
                if ref_index > 0:
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
                newread = gapAlign(newread, gap, start_offset)
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
                if ref_index > 0:
                    prot_errors.append((ref_index, 'i', str(translate(insertion))))
                i += l
                ref_index += 1
            else:  # realign gap and repeat loop at same position to compare the codons
                gap = (gap[0] + i - 1, gap[1] + i - 1)
                newref = gapAlign(newref, gap, start_offset)
                continue

        elif translate(read_codon) != translate(ref_codon):  # must be a substitution
            if ref_index > 0:
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


# Raw processing of all alignments, get composition

def count_one_fraction(alignment, refname, debug, start_offset, end_trail):
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
    # when a protein mutation is first encountered, create an entry
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

        if read_is_wt(read, ref):
            if debug:
                trimmed_read = re.search(r'^-+([AGCTN][ACGTN-]+[ACGTN])-+$', str(read))
                print()
                print(trimmed_read.group(1))
                printErrors("WT", read, ref, True)
            continue

        dna_errors, dna_hgvs, prot_erros = None, None, None

        try:
            dna_errors = find_DNA_diff(read, ref, debug, start_offset, end_trail)  # errors = a tuple
            dna_hgvs = find_DNA_hgvs(read, ref, refname, debug, start_offset, end_trail)  # string in HGVS format (ish)
            prot_errors = find_protein_diff(read, ref, debug, start_offset, end_trail)
            # print()
            # print(readname)
            # print(dna_hgvs, prot_errors)
            # printErrors(dna_errors, read, ref, True)

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
            one_lane_counts[prot_errors] = {'dna': defaultdict(int), 'dna_hgvs': defaultdict(int), 'total': 1}
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


def count_multiple_fractions(folder, baseline, debug, start_offset, end_trail):
    """
    Process all reference.fraction.aln files in given  folder
    :param folder: contains all *.aln
    :param baseline: string containing name of baseline - needed for DNA HGVS format
    :param debug: if True, print reads with mutations
    :param start_offset: number of nt in fasta file before the gene starts
    :param end_trail: number of nt to ignore at the end
    :return:
    """
    all_references = {}

    print(os.listdir(folder))
    for alnfile in os.listdir(folder):
        if alnfile.endswith('.aln'):
            aln_path = os.path.join(folder, alnfile)
            refname, fraction, suffix = alnfile.rsplit(".", 2)
            print('Counting alignment {0} in background {1} and activity fraction {2}'
                  .format(alnfile, refname, fraction))

            if refname not in all_references.keys():
                all_references[refname] = {}

            if fraction == baseline:
                fraction = 'baseline'
            all_references[refname][fraction] = count_one_fraction(aln_path, refname, debug, start_offset, end_trail)

    return all_references


def classify_dna(dna_error):
    if dna_error is None:  # empty or broken reads
        return 'b'
    elif dna_error == ():
        return 'wt'
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


def classify_protein(mutation):
    """
    Detect whether a protein mutation is an insertion/deletion/substituion and if more than 1, (non)consecutive
    :param mutation:
    :return:
    """
    if mutation is None:  # came from empty or broken reads
        return 'b'
    else:
        m = []
        for pos in range(len(mutation)):
            t = mutation[pos][1]
            if t != 'i':
                m.append(t)
            else:
                m.append(t + str(len(mutation[pos][2])))
        # need to distinguish between consecutive mutations and likely sequencing errors
        if 'f' in m:
            return 'f'
        elif len(m) == 1:
            return ''.join(m)
        elif len(m) >= 1:
            c = is_mutation_consecutive(mutation) + '-' + ''.join(m)
            return c


# Start composition statistics

def get_dna_composition(all_references, cutoff=10):
    dna_count = {}
    dna_reads = {}
    for background in all_references.keys():
        for fraction in all_references[background].keys():
            # print('Analysing background {0} and fraction {1}. Unusal mutations: '.format(background, fraction))
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
                            dna_type = classify_dna(dna_error)
                            if dna_type == 'wt':
                                continue
                            elif dna_type == 'b':
                                print(dna_type, dna_error, mut_total)
                            try:
                                dna_count[background + fraction][dna_type] += 1
                                dna_reads[background + fraction][dna_type] += c
                                distinct_mutations += 1
                                total_count += c
                            except KeyError:
                                print(dna_type, dna_error, mut_total)
                                dna_count[background + fraction]['other'] += 1
                                dna_reads[background + fraction]['other'] += c

    return pd.DataFrame.from_dict(dna_count), pd.DataFrame.from_dict(dna_reads)


def is_mutation_consecutive(mutation):
    """
    If only consecutve amino acids are affected, return 'c', else return 'nc'
    :param mutation:
    :return:
    """
    for pos in range(1, len(mutation)):
        if mutation[pos][0] != (mutation[pos - 1][0] + 1):
            return 'nc'
    return 'c'


def get_protein_composition(all_references, cutoff=10):
    protein_count = {}
    protein_reads = {}
    for background in all_references.keys():
        for fraction in all_references[background].keys():
            # print('Analysing protein composition in background {0} and fraction {1}. Unusal mutations: '.format(background, fraction))
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
            # print('In background {0} and fraction {1} found {2} distinct mutations with total read count {3}'.format(
            #     background, fraction, distinct_mutations, total_count
            # ))

    return pd.DataFrame.from_dict(protein_count), pd.DataFrame.from_dict(protein_reads)


def insertion_composition(all_references, cutoff=2, l=(3, 6, 9)):

    comp = {length: {k: {'A': 0, 'C': 0, 'T': 0, 'G': 0} for k in range(length)} for length in l}
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


def transposon_consensus_seq(all_references, reference, fraction='d3', transposon='d3', start_offset=3):
    """
    Determine the consensus sequence for transposon insertion, by analysing the position of d3 or i3 mutations
    :param all_references: dictionary containing all mutation data
    :param reference: BioSeq fasta reference
    :param fraction: name of the activity fraction / library to be analysed
    :param transposon: which type of mutations are we counting
    :param start_offset
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
                n1 = int(dna_mutation[0][0]) + start_offset - 2
                # this position refers to the nt in the first position of  the indel in 1-count:
                # hence this is the 2nd nt of 5 nt site -> need to retrieve start-1:start+4
                # Correct for offset: if start = 1: this is ref[start_offset]
                # Therefore N1 = start_offset + start -2
                if transposon == 'd3':
                    trans_seq = str(ref[n1:n1+5].seq)
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
    parser.add_argument('-s', '--start_offset', help='Number of nt before starting ATG (integer)', required=False,
                        default=3, type=int)
    parser.add_argument('-e', '--end_trail', help='Number of nt after end of gene (integer)', required=False, default=0,
                        type=int)
    parser.add_argument('-o', '--output', help='Filename')
    args = parser.parse_args()

    # On the first run, analyse all *.aln files in the target folder and create a dictionary of errors
    # Structure: all_ref[background][fraction][protein mutation] = all data about the mutations
    # all_ref = count_multiple_fractions(args.folder, args.baseline, args.debug, args.start_offset, args.end_trail)
    # export_hgvs(all_ref, 'new')
    # with open(args.output +'.p', 'wb') as f:
    #     pickle.dump(all_ref, f)

    # Later, load in the results for analysis
    with open(args.output + '.p', 'rb') as f:
        all_ref = pickle.load(f)

    # Now start with statistics
    # 1. Generate CSV files that give overall composition of libraries on protein and DNA level
    for cutoff in [1, 10]:
        dna_count, dna_reads = get_dna_composition(all_ref, cutoff)
        protein_count, protein_reads = get_protein_composition(all_ref, cutoff)
        print(cutoff)
        print(dna_count)
        print(dna_reads)
        print(protein_count)
        print(protein_reads)
    print()

    # 2. Generate data for a histogram of tranposon insertion sites: best used for -3 bp library, it shows how many
    #    times a mutation is detected at each DNA position. Spikes corresponds to sites close to tranposon preferred
    #    insertion sequence (GC rich). Main Fig. 3A
    histogram = find_transposon_histogram(all_ref, 'S6', baseline='d3', transposon='d3')
    print('Histogram of transposon bias')
    pprint.pprint(histogram)
    print()

    # 3. Find TransDel consensus sequence in -3 bp library - gives number of observations for the NNNNN target site
    # This is used for WebLogo in Fig. 3A
    d3_baseline, d3_cons = transposon_consensus_seq(all_ref, args.reference, fraction='d3',
                                                    start_offset=args.start_offset)
    print('TransDel consensus sequence: baseline counts (reflect GC composition)')
    pprint.pprint(d3_baseline)
    print('TransDel consensus counts')
    pprint.pprint(d3_cons)
    print()

    # # 4. Get position by position ACGT composition of insertions - SI Figure
    ins_composition = insertion_composition(all_ref)
    print('ACGT composition of insertions, for each position separately')
    pprint.pprint(ins_composition)
    print()
    #
    # 5. Get data for histogram of how often mutations are detected on average - Fig. 3B Poisson-like distribution
    detection_histogram = dna_mutation_frequencies(all_ref)
    print('How many mutations occur once, twice, etc. per library')
    pprint.pprint(detection_histogram)
    print()
    #
    # # 6. How diverse are insertions at each position? Maximum of 64 from 64 possible triplets in NNN, more in i6/i9
    ins_freq = insertion_frequencies(all_ref)
    print('Diversity of insertions per position')
    pprint.pprint(ins_freq)
    print()
