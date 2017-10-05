#!/usr/bin/env python3
"""
1. Check current directory for *.ab1 files
Then for each file
2. Read in file
3. Convert to FASTQ with quality scores
4. Trim poor quality: http://bi1x.caltech.edu/2015/tutorials/bi1x_bioinformatics.html
5. Align forward & reverse to reference with Bio.align
6. Give alignment to count_one_lane.py
7. Record name of sequencing file & errors present with location.

"""
import argparse
import glob
import os

from collections import defaultdict

import numpy as np

from Bio import SeqIO
from Bio import AlignIO
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable


# from Bio.Seq import MutableSeq


from indels.ind import findErrors, error_to_protein, convert_ab1, needle_align, get_codons

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reference', help='FASTA file with reference sequence', required=True)
    args = parser.parse_args()

    reference = SeqIO.read(args.reference, 'fasta', alphabet=IUPAC.ambiguous_dna)
    codons = get_codons(withdeletions=False)
    codons_with_deletions = get_codons(withdeletions=True)

    # Specify FASTQ directory
    ab1_dir = os.getcwd()
    convert_ab1(ab1_dir)

    listing = listing = glob.glob(os.path.join('fastq', '*_trimmed.fastq'))
    rejected = defaultdict(int)



    for fqname in listing:
        ref, read, id = needle_align(fqname, args.reference)
        errors = findErrors(read, ref, rejected, codons_with_deletions, 720)

        protein = error_to_protein(errors)
        if protein == '':
            protein = "WT"

        print(id, protein, ' '.join(errors), sep='\t')
