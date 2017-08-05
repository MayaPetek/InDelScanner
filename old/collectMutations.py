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
    
def setupAminoAcids():

    acids = "FYC*WLPHQIMTNKSRVADEG"
    aminoAcids = list(acids)
    
    return aminoAcids

def codonTable(codons):
    """
    generate codon table, derived from NCBI standard code

    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
    Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
    Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
    Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
    """
    codons.append("---")
    aminoAcids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG-'
    codonTable = dict(zip(codons, aminoAcids))

    return codonTable
    
def setupCounts(codons,aminoAcids):
    """
    All mutations are counted into a list of dictionaries, one for each position -
    list[0] = first position
    Each position carries a dictionary with counts for that position on DNA
    and/or protein level
    """
    # add deletions to codon and amino acid lists
    codons.append("---")
    aminoAcids.append("-")
    
    # GFP is 720 nt long, barcodes extend further - for now make it longer
    countsDNA = [{triplet: 0 for triplet in codons} for k in range(250)]
    countsProtein = [{acid: 0 for acid in aminoAcids} for k in range(250)]   
     
    return countsDNA, countsProtein
    
def translate(triplet):
# translate triplet DNA to protein according to standard code

    protein = codonTable.get(triplet)
    return protein

def containsMutation(codonLine):
    """
    Check if the mutation is Rejected or an empty line
    """
    return codonLine[0].isnumeric()
       
def classify(codonLine, codons):
    """
    Check for a valid mutation and add appropriate counter to point mutations.
    Currently ignores linked mutations.
    """
    if not containsMutation(codonLine):
        return "Rejected"
    
    
    if len(codonLine) == 4:
    # a single amino acid change
        # see if it's a valid mutation
        if (codonLine[1] or codonLine[3]) not in codons:
            return "Invalid mutation"
        else:
        # a substitution/deletion of the form: [index, triplet, ->, after]
            return "Point"
    else:
        return None


"""
Count mutations into a large table
"""
codons = setupCodons()
codonTable = codonTable(codons)
aminoAcids = setupAminoAcids()

countsDNA, countsProtein = setupCounts(codons,aminoAcids)

rejected = 0
invalid = 0

with open(sys.argv[1]) as f:
    for line in f:
        # ignore empty lines
        if not line.isspace():
            codonLine = line.rstrip().rstrip(";").split()
            if classify(codonLine, codons) == "Point":
                index = int( int(codonLine[0]) / 3)
                countsDNA[index][codonLine[3]] += 1
                countsProtein[index][translate(codonLine[3])] += 1
            elif classify(codonLine, codons) == "Invalid mutation":
                invalid += 1
            elif classify(codonLine, codons) == "Rejected":
                rejected +=1
            else:
                pass

# give results
print("Previously rejected mutations: ", rejected)
print("Invalid mutations, e.g. -TT: ", invalid)
print("List of triplet counts for valid mutations:\n", countsDNA)
print("List of amino acid counts for valid mutations:\n", countsProtein)
