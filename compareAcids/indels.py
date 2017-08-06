#!/usr/bin/python3

import sys

import time
import argparse

from collections import defaultdict

from functions import loadCounts

from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, AlignIO
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable


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


    reference = SeqIO.read(args.reference,'fasta', alphabet = IUPAC.ambiguous_dna)

    try:
        valid_counts = loadCounts(reference)
    except IOError:
        print("Making counts")
        valid_counts = prepareCounts(reference)

    all_counts = defaultdict(int)
    rejected = defaultdict()

    read = None
    ref = None
    
    ctr=0 # counter for DEBUG printing
    
    # reading & looping over read/reference sequence in multiple sequence alignment
    # use AlignIO parser and keep sequence only, allowing it to change (important for gap shifts) 
    for pair in AlignIO.parse(args.alignment,"fasta", alphabet = IUPAC.ambiguous_dna, seq_count = 2):

        ctr += 1

        ref = pair[0].seq.tomutable()
        read = pair[1].seq.tomutable()
       
        
        # is the read broken?
        if verifyRead(read, ref, rejected) == None:
            if PRINT_COLOURED_DIFF:
                printDiff(errors, read, ref, ctr, PRINT_COLOURED_DIFF)
            continue

        errors = findErrors(read, ref, rejected)
        assert errors

        if (len(errors) // 3) > MAX_ERRORS: # len(errors) = 1, 4, 7,...
            rejected['too many errors'] +=1
            if PRINT_COLOURED_DIFF:
                printDiff(errors, read, ref, ctr, PRINT_COLOURED_DIFF)
            continue
        # count up all errors in one massive dictionary
        elif len(errors) > 1:
            all_counts[errors] += 1
            if errors[0] in valid_counts.keys():
                valid_counts[errors[0]].insert(errors)

    return valid_counts, all_counts

    # # Print everything interesting to output
    #
    # with open(args.output + ".csv", "w") as output:
    #     print("Rejected", file = output)
    #     for (key, value) in rejected.items():
    #         print(key, value, sep=",", file = output)
    #     print("Substitutions", file = output)
    #     for (key, value) in subcounts.items():
    #         print(' '.join(key), value, sep=",", file = output)
    #     print("Deletions", file = output)
    #     for (key, value) in delcounts.items():
    #         print(' '.join(key), value, sep=",", file = output)
    #
    # with open(args.output + ".all.csv", "w") as output:
    #     for (key, value) in allcounts.items():
    #         print(' '.join(key), value, sep=",", file = output)
    #
    # with open(args.output + ".d3transposon.csv", "w") as output:
    #     for (key, value) in delcounts.items():
    #         if key.count('---') == 1:
    #             print(key[0], value, sep=",", file = output)
    #
    # time_4 = time.time()
    # print("Printed results in", (time_4 - time_3) // 60, "min and", (time_4 - time_3) % 60, "s")


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
    parser.add_argument('-d', '--debug', help='Turn on debugging', required=False, action="store_true")
    args = parser.parse_args()

    # Visual representation of how reads match reference
    PRINT_COLOURED_DIFF = args.debug

    valid_counts, all_counts = countOneAlignment()