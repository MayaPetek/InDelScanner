#!/usr/bin/python3

import sys
import pickle
import argparse

from collections import defaultdict

from functions import loadCounts, saveCounts, prepareCounts, findErrors, verifyRead
from outputUtils import printDiff

from Bio import SeqIO, AlignIO
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio.Seq import MutableSeq


def countOneAlignment():
    """
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

    codons = list(CodonTable.unambiguous_dna_by_name["Standard"].forward_table.keys())
    codons += CodonTable.unambiguous_dna_by_name["Standard"].stop_codons

    reference = SeqIO.read(args.reference, 'fasta', alphabet=IUPAC.ambiguous_dna)

    MAX_ERROR_INDEX = len(reference) - 1

    try:
        valid_counts = loadCounts(reference, '.counts.p')
        print("Imported counts")
    except IOError:
        print("Making counts")
        valid_counts = saveCounts(reference, '.counts.p')

    all_counts = defaultdict(int)
    rejected = defaultdict(int)

    ctr=0 # counter for DEBUG printing
    
    # reading & looping over read/reference sequence in multiple sequence alignment
    # use AlignIO parser and keep sequence only, allowing it to change (important for gap shifts) 
    for pair in AlignIO.parse(args.alignment,"fasta", alphabet=IUPAC.ambiguous_dna, seq_count=2):

        ctr += 1

        ref = pair[0].seq.tomutable()
        read = pair[1].seq
        read = MutableSeq(str(read).replace('N', '.'), read.alphabet)
        id = pair[1].id

        errors = findErrors(read, ref, rejected, codons, MAX_ERROR_INDEX) # errors = a tuple

        # is the read broken?
        if not verifyRead(read, ref, rejected, MATCH_N_END):
            if args.debug:
                printDiff(errors, read, ref, id, ctr, args.debug)
            continue

        if (len(errors) // 3) > MAX_ERRORS: # len(errors) = 1, 4, 7,...
            rejected['too many errors'] +=1
            if args.debug:
                printDiff(errors, read, ref, id, ctr, args.debug)
            continue

        # count up all errors in one massive dictionary
        elif len(errors) > 1:
            all_counts[errors] += 1
            if errors[0] in valid_counts.keys():
                    valid_counts[errors[0]].insert(errors)

    return valid_counts, all_counts

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
        description='A script to find and counts all 3 bp substitutions and 3/6/9 bp deletions in a gene from a multiple alignment')
    parser.add_argument('-a', '--alignment', help='Multiple sequence alignment', required=True)
    parser.add_argument('-r', '--reference', help='Reference fasta file with with the alignment was constructed',
                        required=True)
    parser.add_argument('-o', '--output', help='Output file name', required=True)
    parser.add_argument('-d', '--debug', help='Turn on debugging', required=False, action="store_true") # Visual representation of how reads match reference
    args = parser.parse_args()

    # max number of mutations called in a read
    MAX_ERRORS = 4

    # How many nucleotides need to match at the end of a read for a valid alignment:
    # - matching 2 should correct alignment errors, 3 avoids problems with InDel
    #   repositioning
    # - 3 also simplifies handling the first codon: it's either complete or it's OK
    #   to move 1 or 2 bases over to the next triplet
    MATCH_N_END = 3

    valid_counts, all_counts = countOneAlignment()

    with open(args.output + '.valid_counts.p', 'wb') as f:
        pickle.dump(valid_counts, f)
    with open(args.output + '.all_counts.p', 'wb') as f:
        pickle.dump(all_counts, f)

    for k, v in valid_counts['d'].items():
        if v > 0:
            print(k,v)
