#!/usr/bin/python3

import sys
import pickle
import re
import os
import csv
import argparse

import pandas as pd
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


def find_DNA_diff(read, ref):
    """
    @ read, ref: MutableSeq objects
    :return errors - tuple (position, expected triplet, actual triplet, ) / none if broken read

    Letter by letter report mutations in NGS read, all counts 1- based in result (code in 0-count).
    - substitution: 78C = nt 78 in reference is changed to C
    - deletions: 78d6 = 6 nt deleted after 78: 1-78, d6, 85-end
    - insertion: 78iATC = after nt 78 inserted seq. ATC
    """

    if read is None:
        print('no read provided')
        return

    # No gap realignment at this point

    # quality control that there are no mutations at ends of reads
    ends = findEnds(read, ref)
    if not endMatch(read, ref, ends):
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
                dna_errors += [str(ref_index) + 'd' + str(l)]  # deletion length l starting at ref_index in 0-count
                i += l
                ref_index += l
            else:
                dna_errors += [str(ref_index) + 'f']
                break

        elif ref[i] == '-':
            # start of an insertion
            l = indel_len(ref, i)
            # check for frameshifts
            if l % 3 == 0:
                dna_errors += [str(ref_index) + 'i' + str(read[i:i+l])]
                i += l
            else:
                dna_errors += [str(ref_index) + 'f']
                break

        else:
            # substitution
            dna_errors += [str(ref_index + 1) + str(read[i])]
            i += 1
            ref_index += 1

    return tuple(dna_errors)


def find_protein_diff(read, ref):

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
    ref_index = int(ends.get('aligned') / 3) + 1 # reference amino acid index

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
                prot_errors.append(str(ref_index) + 'f')
                return tuple(prot_errors)
            if read_codon == '---':  # single codon deletion
                prot_errors += [str(ref_index) + 'd']
                i += 3
                ref_index += 1

            else:  # check it's not a frame shift
                l = indel_len(newread, i)
                if l % 3 != 0:
                    prot_errors.append(str(ref_index) + 'f')
                    return tuple(prot_errors)
                # realign gap and repeat loop at same position to compare the codons
                gap = findGap(newread[i - 1:])
                gap = (gap[0] + i - 1, gap[1] + i - 1)
                newread = gapAlign(newread, gap)
                continue

        elif '-' in ref_codon:  # found an insertion
            l = indel_len(newref, i)
            if l % 3 != 0:
                prot_errors.append(str(ref_index) + 'f')
                return tuple(prot_errors)
            gap = findGap(newref[i-1:])
            if gap[0] == 1:  # insertion after codon
                insertion = newread[gap[0] + i - 1:gap[1] + i - 1]
                if '-' in insertion:
                    prot_errors.append(str(ref_index) + 'f')
                    return tuple(prot_errors)
                prot_errors.append(str(ref_index) + 'i' + str(translate(insertion)))
                i += l
                ref_index += 1
            else:  # realign gap and repeat loop at same position to compare the codons
                gap = (gap[0] + i - 1, gap[1] + i - 1)
                newref = gapAlign(newref, gap)
                continue

        elif translate(read_codon) != translate(ref_codon):  # must be a substitution
            prot_errors.append(str(ref_index) + str(translate(read_codon)))
            i += 3
            ref_index += 1

        else:
            i += 3
            ref_index += 1

    return tuple(prot_errors)


def countOneAlignment(alignment, debug):
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
    one_lane_counts = defaultdict(partial(defaultdict, int))

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
            dna_errors = find_DNA_diff(read, ref)  # errors = a tuple
            prot_errors = find_protein_diff(read, ref)
        except:
            if not dna_errors:
                print(dna_errors)
            print_coloured_diff(readname, read, ref, debug)
            raise

        one_lane_counts[prot_errors]['total'] += 1
        one_lane_counts[prot_errors][dna_errors] += 1

    n = 0
    threshold = 10
    for error in one_lane_counts.keys():
        if one_lane_counts[error]['total'] > threshold:
            n += 1

    print('Fount {0} total protein mutations, of which {1} have more than {2} counts'
          .format(len(one_lane_counts), n, threshold))

    return one_lane_counts


def count_multiple_fractions(folder, debug):
    """
    Process all reference.fraction.aln files in given  folder
    :param folder:
    :return:
    """
    all_references = {}

    for aln in os.listdir(folder):
        if aln.endswith('.aln'):
            ref, fraction, suffix = aln.rsplit(".", 2)
            print('Counting alignment {0} in background {1} and activity fraction {2}'
                  .format(aln, ref, fraction))
            if ref not in all_references.keys():
                all_references[ref] = {}
            all_references[ref][fraction] = countOneAlignment(aln, debug)

    return all_references


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
                dna_errors = find_DNA_diff(read, ref)  # errors = a tuple
                protein_errors = find_protein_diff(read, ref)
                if verbose:
                    print(row['ab1'], dna_errors, protein_errors, row['protein'])
                if dna_errors:
                    experimental[protein_errors][dna_errors] = float(row['activity'])
            except ValueError:
                prefix, suffix = os.path.splitext(fqname)
                outname = prefix + '.aln'
                if os.stat(outname).st_size == 0:
                    print(outname, 'is empty.')
                    continue


def import_sanger(sanger_folder):
    """
    Expect all sequencing files to be in a folder
    :param sanger_folder:
    :return:
    """

    # set up experimental dictionary in same format as one lane counts
    # experimental[protein_mutation][dna_mutation] = activity
    experimental = defaultdict(partial(defaultdict, int))


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
            process_sanger_plate(summary, experimental, rc=True)

    # convert to data frame

    os.chdir(cwd)
    return experimental



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
    args = parser.parse_args()

    """
    PROCESS SANGER SEQUENCING DATA
    """

    experimental = import_sanger(os.path.join(args.folder, 'sanger'))

    """
    PROCESS NGS DATA
    """

    # all_references = count_multiple_fractions(args.folder, args.debug)
    # with open('everything.p', 'wb') as f:
    #     pickle.dump(all_references, f)

    with open('everything.p', 'rb') as f:
        all_references = pickle.load(f)

    eGFP = combine_totals_same_reference(all_references['eGFP'])
    GFP8 = combine_totals_same_reference(all_references['GFP8'])

    for e in experimental.keys():
        print('Error: {} \t Activity: {} \n eGFP:{} \n GFP8:{} \n'.format(e, experimental[e], eGFP[e], GFP8[e]))

    # eGFP_shared, GFP8_shared = find_shared_entries(eGFP, GFP8)
    #
    # df_eGFP = pd.DataFrame.from_dict(eGFP_shared, orient='index')
    # df_GFP8 = pd.DataFrame.from_dict(GFP8_shared, orient='index')



    # # both = pd.merge(df_eGFP, df_GFP8, how='inner', sort=False)

