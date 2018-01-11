#!/usr/bin/python3

import sys
import pickle
import argparse
from collections import defaultdict

from Bio import SeqIO, AlignIO
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq

from indels.ind import trim_read, findEnds, endMatch, find_dna_mutations
from indels.output import print_coloured_diff


def find_DNA_diff(read, ref):
    """
    @ read, ref: MutableSeq objects
    :return errors - tuple (position, expected triplet, actual triplet, ) / none if broken read

    We've got to find and report several types of difference: insertion, deletion, and substitution.
    We'll begin by running through the strings until we find the first 3-letter block that contains a "-"
    in either string, or which has letters in both strings, but they differ.
    After that, we will regard any "-" in read as being part of a deletion error, any "-" in ref as
    being part of an insertion error, until we reach the point where all remaining symbols in read are
    # "-" (at which point we are finished).
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
            l = 0
            while read[i+l] == '-':
                l += 1
            # now we know the length of a deletion, check for frameshifts
            if l % 3 == 0:
                dna_errors += [str(ref_index) + 'd' + str(l)] # deletion length l starting at ref_index in 0-count
                i += l
                ref_index += l
            else:
                dna_errors += [str(ref_index) + 'f']
                break

        elif ref[i] == '-':
            # start of an insertion
            l = 0
            while ref[i+l] == '-':
                l += 1
            # check for frameshifts
            if l % 3 == 0:
                dna_errors += [str(ref_index) + 'i' + str(read[i:i+l])]
                i += l
            else:
                dna_errors += [str(ref_index) + 'f']
                break

        else:
            # substitution
            dna_errors += [str(ref_index) + str(read[i])]
            i += 1
            ref_index += 1

    return tuple(dna_errors)

def find_protein_diff():

    return

def countOneAlignment(args):
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
    # Demand Python 3.
    if sys.version_info[0] < 3:
        print("Python 3 is required, but you are using Python %i.%i.%i") % (
            sys.version_info[0], sys.version_info[1], sys.version_info[2])
        sys.exit(1)

    reference = SeqIO.read(args.reference, 'fasta', alphabet=IUPAC.ambiguous_dna)

    ctr = 0 # counter for DEBUG printing
    one_lane_counts = defaultdict(int)
    rejected = defaultdict(int)

    
    # reading & looping over read/reference sequence in multiple sequence alignment
    # use AlignIO parser and keep sequence only, allowing it to change (important for gap shifts) 
    for pair in AlignIO.parse(args.alignment, "fasta", alphabet=IUPAC.ambiguous_dna, seq_count=2):

        ctr += 1

        # both read and ref are MutableSeq
        ref = pair[0].seq.tomutable()
        read = pair[1].seq.tomutable()
        read = MutableSeq(str(read).replace('N', '.'), read.alphabet)
        readname = pair[1].id

        # trim sequencing read to reference
        ref, read = trim_read(ref, read)

        dna_errors = find_DNA_diff(read, ref) # errors = a tuple
        print(dna_errors)
        one_lane_counts[dna_errors] += 1

    return one_lane_counts


if __name__ == "__main__":
    """
    Arguments:
    1.   Alignment of sequencing reads to a reference sequence in FASTA format. The
         reference comes before read sequence. (>Ref, seq, >Read, seq).
    2.   Reference sequence in FASTA format. Same file that was used to create the 
         alignment, in particular the reference name here needs to match reference
         name in the alignment.
    3.   Name of output file
    (4.) Optional: DEBUG print a coloured representation of mismatches
    """
    parser = argparse.ArgumentParser(
        description='Finds all 3 bp substitutions and 3/6/9 bp deletions in a gene')
    parser.add_argument('-a', '--alignment', help='Multiple sequence alignment', required=True)
    parser.add_argument('-r', '--reference', help='Reference fasta file with with the alignment was constructed',
                        required=True)
    parser.add_argument('-o', '--output', help='Output file name', required=False)
    parser.add_argument('-d', '--debug', help='Turn on debugging', required=False, action="store_true") # Visual
    args = parser.parse_args()

    one_lane_counts = countOneAlignment(args)

    #  Libraries contain up to range from 1 to 4 mutations per gene
    #  shortest = single deletion, longest = 1 substitution + 3 deletions
    # so total of 4 errors

    # How many nucleotides need to match at the end of a read for a valid alignment:
    # - matching 2 should correct alignment errors, 3 avoids problems with InDel
    #   repositioning
    # - 3 also simplifies handling the first codon: it's either complete or it's OK
    #   to move 1 or 2 bases over to the next triplet

    # valid_counts, bad_counts, rejected = countOneAlignment(MAX_ERRORS=4, MATCH_N_END=3)
    #
    # with open(args.output + '.valid_counts.p', 'wb') as f:
    #     pickle.dump(valid_counts, f)
    # with open(args.output + '.bad_counts.p', 'wb') as f:
    #     pickle.dump(bad_counts, f)

    n = 0
    for error, count in one_lane_counts.items():
        if count > 0:
            n +=1
            print(error, count)
    print('Found', n, 'mutations')