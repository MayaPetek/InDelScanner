#!/usr/bin/python3

import sys
import pickle
import re
import os
import argparse

import pandas as pd
from collections import defaultdict
from functools import partial

from Bio import AlignIO
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq, translate

from indels.ind import trim_read, findEnds, endMatch, findGap, gapAlign
from indels.output import print_coloured_diff


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
        return

    # No gap realignment at this point

    # quality control that there are no mutations at ends of reads
    ends = findEnds(read, ref)
    if not endMatch(read, ref, ends):
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
    ref_index = int(ends.get('aligned') / 3)  # reference amino acid index

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
            print_coloured_diff(readname, read, ref, verbose)
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

    df_protein = pd.DataFrame.from_dict(one_reference_counts, orient='index')

    return df_protein

if __name__ == "__main__":
    """
    Arguments:
    1.   Alignment of sequencing reads to a reference sequence in FASTA format. The
         reference comes before read sequence. (>Ref, seq, >Read, seq).
    (2.) Optional: DEBUG print a coloured representation of mismatches
    """
    # Demand Python 3.
    if sys.version_info[0] < 3:
        print("Python 3 is required, but you are using Python %i.%i.%i") % (
            sys.version_info[0], sys.version_info[1], sys.version_info[2])
        sys.exit(1)

    parser = argparse.ArgumentParser(
        description='Finds all in-frame mutations in a gene')
    parser.add_argument('-f', '--folder', help='Folder containing multiple sequence alignments', required=True)
    parser.add_argument('-d', '--debug', help='Turn on debugging', required=False, action="store_true")  # Visual
    args = parser.parse_args()

    all_references = count_multiple_fractions(args.folder, args.debug)
    with open('everything.p', 'wb') as f:
        pickle.dump(all_references, f)

    total = {}
    for i in all_references.keys():
        total[i] = combine_totals_same_reference(all_references[i])

    print(total.keys())
