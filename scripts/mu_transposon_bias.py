#!/usr/bin/python3

import sys
import csv
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

# Demand Python 3.
if sys.version_info[0] < 3:
    print ("Python 3 is required, but you are using Python %i.%i.%i") % (
           sys.version_info[0], sys.version_info[1], sys.version_info[2])
    sys.exit(1)

"""
1. Read in reference sequence
2. Read in list of most frequent deletions from the counts file. This gives pair of position & 3nt deletion counts for the most common positions
   These are selected to be the most common and unambigious.
3. The relevant transposon insertion site is 1 nt before start of deletion, 3 bp of deletion and 1 extra nt.
   Eg. for CSV row 223, 37, we want bases reference[222:227], that is for X -> [x-1:x+4]
4. Print these N nt to fasta as many times as the original deletion was observed
(5. Later - create consensus sequence weighed by number of occurences)
"""


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='A script to find and counts all 3 bp substitutions and 3/6/9 bp deletions in a gene from a multiple alignment')
    parser.add_argument('-c', '--counts', help='CSV file with 1 column giving position and 2nd counts of deletions',
                        required=True)
    parser.add_argument('-r', '--reference', help='Reference fasta file with with the alignment was constructed',
                        required=True)
    parser.add_argument('-o', '--output', help='Output file name', required=True)
    args = parser.parse_args()

    """ Read in a list of counts per position in the format:
    0-based position, count of d3 at this position in eGFP, count in GFP8
    All counts have been corrected for coverage
    """
    ref = SeqIO.read(args.reference,'fasta', alphabet = IUPAC.ambiguous_dna) #SeqRecord

    with open(args.counts) as f:
        counts = csv.DictReader(f, dialect = 'excel')
        pos, num = counts.fieldnames
        mu = {}
        for line in counts:
            mu[int(line.get(pos))] = int(line.get(num))

        
    with open(args.output, 'w') as out:
        for position, count in mu.items():
            if count > 0:
                for i in range(count):
                    # interesting is [pos-1:pos+4], add 5 nt each side to be sure
                    SeqIO.write(ref[position-6:position+9], out, 'fasta')
            



