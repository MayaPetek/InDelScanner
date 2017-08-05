#!/usr/bin/python3

import sys

# Demand Python 3.
if sys.version_info[0] < 3:
    print ("Python 3 is required, but you are using Python %i.%i.%i") % (sys.version_info[0], sys.version_info[1], sys.version_info[2])
    sys.exit(1)

import re
from stringUtils import getChar, getChars, setChar, stealLetter
from outputUtils import printErrors


# Count ones with non-multiple of 3 gaps
# Steal symbols across %3 gaps to align gaps and stuff.
# Also output the gap size.
# Reject if more than one gap. This happens rarely and only when the alignment is irredeemably broken.

MAX_ERRORS = 15


# Overwrite this many symbols at the start/end of the read string with "-". This apparently helps work around
# some flavours of aligner unreliability.

IGNORE_FIRST_N_SYMBOLS = 0;
IGNORE_LAST_N_SYMBOLS = 0;

# Require N end symbol of the read to match reference. If there are mutations at the very end of a read, it is
# impossible to tell what the real mutation was (substitution, misplaced InDel, combination)
# Matching 2 letters should be sufficient for alignment quality, 3 is a more cautious choice
# Matching 3 prevents problems with letter stealing

MATCH_N_END = 3;

PRINT_COLOURED_DIFF = len(sys.argv) >= 4 and sys.argv[3] == "DEBUG";

def readReference():
    """
    Read sequence from fasta file and add Ns at the end
    """
    with open(sys.argv[2]) as f:
    # Discard the first line of the file.
        referenceName = f.readline().rstrip()
        referenceSequence = f.readline().rstrip().upper()
        referenceSequence = [referenceSequence[i:i+3] for i in range(0, len(referenceSequence), 3)]
        referenceProtein = [codonTable.get(triplet) for triplet in referenceSequence]
        return referenceName, referenceSequence, referenceProtein

def readUntil(f, sentinel):
    """
    Read and append lines from f until a line containing `sentinel` is reached (or we run out of file).
    """
    out = ""
    line = f.readline().rstrip()
    while ((not line.startswith(sentinel)) and line != ""):
        out += line;
        line = f.readline().rstrip()

    return out;

def cullEnds(readMatch):
    # Overwrite the first few non-dashes with dashes.

    # Extract the entire aligned region.
    match = re.search('([ATGC]+(-+(?=[ATGC]))?)+', readMatch)

    # Trim it. We trim off the number we're asked to trim, as well as the number needed to trim off any partial
    # triplets.
    for i in range(0, IGNORE_FIRST_N_SYMBOLS):
        readMatch = setChar(readMatch, "-", match.start() + i);

    # Trim off the number from the end we're asked for. We can't trim to alignment now, since errors in the middle
    # might alter what we want to do.
    for i in range(0, IGNORE_LAST_N_SYMBOLS):
        readMatch = setChar(readMatch, "-", match.end() - 1 - i);

    return readMatch;

def getStartPoint(readMatch, refMatch):
    """
    Figure out where the first symbol of interest in the read is. This start point should be:
    - Not a -.
    - Not (yet) triplet aligned
    """
    startPoint = 0
    for i in range(0, len(refMatch)):
        if readMatch[i] != "-":
            startPoint = i
            break
    return startPoint

def getEndPoint(readMatch, refMatch):
    """
    Figure out where the last symbol of interest in the read is.
    """
    endPoint = len(refMatch)
    for i in range(0, len(refMatch)):
        if readMatch[-1-i] != "-":
            endPoint = len(refMatch) - i
            break
    return endPoint

def endMatch(readMatch, refMatch, startPoint, endPoint):
    """
    Aligner errors arise when mutations are at the ends of the read rather than in the middle.
    Trimming ends only shifts the problem. Instead require that read ends match the reference.
    Return False if either end doesn't match
    """

    # does the start mismatch
    start = getChars(readMatch, startPoint, MATCH_N_END).casefold() == getChars(refMatch, startPoint, MATCH_N_END).casefold()
    # does the end mismatch
    end = getChars(readMatch, endPoint - MATCH_N_END, MATCH_N_END).casefold() == getChars(refMatch, endPoint - MATCH_N_END, MATCH_N_END).casefold()
    return start and end

def alignStartPoint(startPoint):
    """
    Move start point to a triplet-aligned position
    """

    # Hop back to the last multiple of 3...
    while (startPoint % 3 != 0):
        startPoint -= 1;

    return startPoint

def copyFirstAcid(readMatch, refMatch, startPoint):
    # Handle the first amino acid (which we might only partially have in the read. We don't want to report it
    # as an error if some prefix is simply absent in the read.
    expectedAcid = getChars(refMatch, startPoint, 3);
    actualAcid = getChars(readMatch, startPoint, 3);

    # Steal letters from expectedAcid to fill in the prefix of actualAcid...
    for i in range(0, 3):
        if actualAcid[i] != "-":
            break;

        # Change character i of actualAcid to match expectedAcid.
        actualAcid = actualAcid[:i] + expectedAcid[i] + actualAcid[i + 1:]

    # Correct the read
    readMatch = readMatch[:startPoint] + actualAcid + readMatch[startPoint + 3:]

    return readMatch

def hasMultipleGaps(readMatch):
    """
    Return true iff readMatch contains multiple gaps.
    """

    # This regex looks for 3 distinct islands of letters. That implies at least two gaps.
    return re.search('[ATGC]+-+[ATGC]+-+[ATGC]+', readMatch) != None

def findGap(readMatch):
    """
    Get the start and end index (as a 2-tuple) of the gap in readMatch, or None if there isn't one.
    """

    # Find the gap, captured in group 1. This gives us the indexes in the string where the gap starts and
    # ends. (we define "gap" as "dashes in between letters". We previously used hasMultipleGaps to reject
    # any strings which have multiple gaps (which would confuse this regex)).
    match = re.search('[ATGC]+(-+)[ATGC]+', readMatch)
    if match == None:
        # No gap exists.
        return None

    return (match.start(1), match.end(1))


def gapAlign(readMatch, gap):
    """
    Perform the "letter-stealing" operation:
    - Find the gap (if there is one).
    - While the gap startpoint is not aligned to a triplet boundary, move letters from the end of the gap
      to the start of the gap.

    @param readMatch Aligned read to process.
    @param gap The gap, as returned by findGap.
    """

    old = readMatch

    # Shift letters from the end to the start...
    while gap[0] % 3 != 0:
        assert(readMatch[gap[0]] == "-")

        # This means we ran out of symbols to steal before we managed to align.
        # This indicates a gap at the end, and one we can't fix. Reject!
        if readMatch[gap[1]] == "-":
            return None

        c = readMatch[gap[1]]
        readMatch = setChar(readMatch, c, gap[0])
        readMatch = setChar(readMatch, "-", gap[1])

        # Shift the gap to the right...
        gap = (gap[0] + 1, gap[1] + 1)

    return readMatch


def setupCodonsAminoAcids():

    bases = ['T', 'C', 'A', 'G']
    codons = [a+b+c for a in bases for b in bases for c in bases]
    codons.append("---")
    aminoAcids = list("FYC*WLPHQIMTNKSRVADEG-")
    
    return codons, aminoAcids
    

def codonTable(codons):
    """
    generate codon table, derived from NCBI standard code

    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
    Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
    Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
    Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
    """

    acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG-'
    codonTable = dict(zip(codons, acids))

    return codonTable
 
    
def setupCounts(codons, aminoAcids, mutations):
    """
    All mutations are counted into a list of dictionaries, one for each position
    Each position carries a dictionary with counts for that position on DNA
    and/or protein level
    """
   
    # GFP is 720 nt long, barcodes extend further - for now make it longer
    countsDNA = [{triplet: 0 for triplet in codons} for k in range(len(referenceSequence))]
    countsProtein = [{acid: 0 for acid in aminoAcids} for k in range(len(referenceSequence))]
    countsTotal = [{mutation: 0 for mutation in mutations} for k in range(len(referenceSequence))]
     
    return countsDNA, countsProtein, countsTotal


def classifyPoint(mutation, codons):
    # each mutation is a dictionary with keys 'position', 'expected' and 'actual'
    if (mutation.get('expected') not in codons) or (mutation.get('actual') not in codons):
        return None
    elif mutation.get('actual') == "---":
        return "---"
    elif mutation.get('expected') == "---":
        return "INS"
    else:
        return "NNN"

    
def classify(errors, codons):
    """
    Check for a valid mutation and add appropriate counter to point mutations.
    Currently ignores linked mutations.
    """    
    for mutation in errors:
        if not classifyPoint(mutation, codons):
            return "Invalid mutation"
    
    # check that mutations are sequential
    actualPositions = [errors[n].get('position') for n in range(len(errors))]
    expectedPositions = [actualPositions[0] + 3*n for n in range(len(actualPositions))]
    if not actualPositions == expectedPositions:
        return "Non sequential"
    
    # classify each point mutation and add up the results
    # eg. "NNN ---" is substitution followed by deletion
    mutationTypes = []
    for mutation in errors:
        mutationTypes.append(classifyPoint(mutation, codons))

    return ''.join(mutationTypes)



def summarizeChanges(readMatch, refMatch, codonTable, countsDNA, countsProtein, codons, mutations):
    """
    Print out the ammino acid changes for a read.
    """
    errors = []

    # Reject reads with multiple gaps
    if hasMultipleGaps(readMatch):
        print("Rejected (multigap)")
        return

    readMatch = cullEnds(readMatch)
    startPoint = getStartPoint(readMatch, refMatch)
    endPoint = getEndPoint(readMatch, refMatch)

    # Reject reads with errors at the ends
    if not endMatch(readMatch, refMatch, startPoint, endPoint):
        print("Rejected (ends misalign)")
        return

    # Copy first acid
    startPoint = alignStartPoint(startPoint)
    readMatch = copyFirstAcid(readMatch, refMatch, startPoint)

    # Reject reads with gaps that aren't a multiple of 3 in length.
    gap = findGap(readMatch);
    if gap != None:
        gapSize = (gap[1] - gap[0])

        if gapSize % 3 != 0:
            print("Rejected (gap of size %i)" % gapSize)
            return

        # Perform letter-stealing across the gap (if any). The resulting modified read will be ready for naive
        # triplet-wise comparison.
        readMatch = gapAlign(readMatch, gap)
        if readMatch == None:
            print("Rejected (unalignable gap at end)")
            return

    # We've got to find and report several types of difference: insertion, deletion, and substitution.
    # We'll begin by running through the strings until we find the first 3-letter block that contains a "-"
    # in either string, or which has letters in both strings, but they differ.
    # After that, we will regard any "-" in readMatch as being part of a deletion error, any "-" in refMatch as
    # being part of an insertion error, until we reach the point where all remaining symbols in readMatch are
    # "-" (at which point we are finished).
    for i in range(startPoint, len(refMatch), 3):
        # Check if we've run out of symbols (the rest of the string is just dashes).
        suffix = readMatch[i:];
        if re.search('[ATGC]', suffix) == None:
            break

        # Read this triplet from both strings.
        expectedAcid = getChars(refMatch, i, 3);
        actualAcid = getChars(readMatch, i, 3);

        # Check if this is the last acid, and it's incomplete, ignore it.
        suffix = readMatch[i + 3:];
        if re.search('[ATGC]', suffix) == None and "-" in actualAcid:
            continue


        # Compare the triplets!
        if expectedAcid != actualAcid:
            errors.append({
                'position': i,
                'expected': expectedAcid,
                'actual': actualAcid
            });

    if len(errors) > MAX_ERRORS:
        print("Rejected (errcount)");
    elif len(errors) > 0:
        printErrors(errors, readMatch, refMatch, PRINT_COLOURED_DIFF)
        index = int( int(errors[0].get('position')) / 3)
        if classify(errors, codons) in mutations:
            countsTotal[index][classify(errors, codons)] += 1
        if classify(errors, codons) in ["NNN", "---"]:
            countsDNA[index][errors[0].get('actual')] += 1
            countsProtein[index][codonTable.get(errors[0].get('actual'))] += 1


def printCounts(inputFileName, referenceSequence, countsList, mutationNames, outputFileName):

    with open(inputFileName + '.'+ str(outputFileName), "w") as f:
        # write counts to file
        print("# POSITION WT ", '  '.join(mutationNames), file = f)
        for position in range(len(countsList)):
            values = [position+1, referenceSequence[position]]
            for entry in mutationNames:
                values.append(str(countsList[position].get(entry)))
            print(' '.join(str(number) for number in values), file = f)  


# Analyze the output of `pairwiseCheck.sh`, emitting a handy summary of ammino acid changes.
ctr = 0

codons, aminoAcids = setupCodonsAminoAcids()
codonTable = codonTable(codons)

referenceName, referenceSequence, referenceProtein = readReference()

mutations = ['NNN', '---', 'NNN---', 'NNNNNN','------', '---------','NNN------', 'NNN---------'] 

countsDNA, countsProtein, countsTotal = setupCounts(codons,aminoAcids, mutations)

with open(sys.argv[1]) as f:
    # Discard the first line of the file.
    f.readline();

    while (True):
        if PRINT_COLOURED_DIFF:
            print("\n\n####################################################################################################")
            print("Read: %i \n" % ctr);

        # Read the read and reference blocks for the next one...
        refMatch = readUntil(f, ">")
        readMatch = readUntil(f, referenceName)

        assert(len(refMatch) == len(readMatch));

        # Did we run out of file?
        if readMatch == "":
            break;

        summarizeChanges(readMatch, refMatch, codonTable, countsDNA, countsProtein, codons, mutations)

        ctr += 1;

inputFileName = str(sys.argv[1]).rstrip(".aln")
"""
printCounts(inputFileName, referenceSequence, countsDNA, codons, "countsDNA")
printCounts(inputFileName, referenceProtein, countsProtein, aminoAcids, "countsProtein")
printCounts(inputFileName, referenceSequence, countsTotal, mutations, "countsTotal")
"""
