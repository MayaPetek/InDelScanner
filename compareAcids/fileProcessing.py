#!/usr/bin/python3

import sys

from Bio import SeqIO
from Bio.Alphabet import IUPAC

from functions import prepareCounts, saveCounts, loadCounts, allowedDict, findErrors, classifyPoint, verifyRead

if __name__ == "__main__":
    ref = SeqIO.read(sys.argv[1],'fasta', alphabet = IUPAC.ambiguous_dna)
    saveCounts(ref, '.counts.p')
    c = loadCounts(ref, '.counts.p')