#!/usr/bin/python3

import sys, re

# Demand Python 3.
if sys.version_info[0] < 3:
    print ("Python 3 is required, but you are using Python %i.%i.%i") % (sys.version_info[0], sys.version_info[1], sys.version_info[2])
    sys.exit(1)

# expected arguments: ./collectMutations.py filewithmutations

def setupCodons():

    bases = ['T', 'C', 'A', 'G']
    codons = [a+b+c for a in bases for b in bases for c in bases]
    
    return codons
    
def setupAminoAcids()

    acids = "FYC*WLPHQIMTNKSRVADEG"
    acidList = list(acids)
    
    return aminoAcids

def codonTable(codons):
    """
    generate codon table, derived from NCBI standard code

    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
    Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
    Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
    Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
    """

    aminoAcids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codonTable = dict(zip(codons, amino_acids))

    return codonTable
    
def setupCounts(codons,aminoAcids):
    """
    All mutations are counted into a list of dictionaries, one for each position -
    list[0] = first position
    Each position carries a dictionary with counts for that position on DNA
    and/or protein level
    """
    # GFP is 720 nt long, barcodes extend further - for now make it longer
    countsDNA = [{triplet: 0 for triplet in codons} for k in range(798, step=3)]
    countsProtein = [{acid: 0 for triplet in aminoAcids} for k in range(250)]    
    return [countsDNA, countsProtein]

def translate(triplet):
# translate triplet DNA to protein according to standard code

    protein = codonTable.get(triplet)
    return protein


def classify(codonMutation):

    if "Rejected" in codonMutation:
        break
    elif len(codonMutation) == 4 and "---" not in codonMutation
        # this is a single triplet subsitution
        
        
def getPosition()

def count()




# put everything together

codonTable = setupCodons()

with open(sys.argv[1]) as f:
        
        while (True):
            # Read the next mutation description
            codonMutation = f.readline().rstrip().split()
            # Did we run out of file?
                if codonMutation == "":
                break;

            # do other stuff
