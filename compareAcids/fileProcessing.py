#!/usr/bin/python3

import sys


from Bio import SeqIO
from Bio.Alphabet import IUPAC

from functions import prepareCounts, saveCounts, loadCounts, allowedDict, findErrors, classifyPoint, verifyRead

# max number of mutations called in a read
MAX_ERRORS = 4

# How many nucleotides need to match at the end of a read for a valid alignment:
# - matching 2 should correct alignment errors, 3 avoids problems with InDel
#   repositioning
# - 3 also simplifies handling the first codon: it's either complete or it's OK
#   to move 1 or 2 bases over to the next triplet
MATCH_N_END = 3

if __name__ == "__main__":
    ref = SeqIO.read(sys.argv[1],'fasta', alphabet = IUPAC.ambiguous_dna)
    saveCounts(ref)
    c = loadCounts(ref)