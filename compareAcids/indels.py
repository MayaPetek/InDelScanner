#!/usr/bin/python3

import sys
import re
import time
import argparse

from pprint import pprint
from collections import defaultdict


from stringUtils import getChar, getChars
from outputUtils import printErrors

from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, AlignIO
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable


# Demand Python 3.
if sys.version_info[0] < 3:
    print ("Python 3 is required, but you are using Python %i.%i.%i") % (
           sys.version_info[0], sys.version_info[1], sys.version_info[2])
    sys.exit(1)


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
parser = argparse.ArgumentParser(description='A script to find and counts all 3 bp substitutions and 3/6/9 bp deletions in a gene from a multiple alignment')
parser.add_argument('-a','--alignment', help='Multiple sequence alignment',required=True)
parser.add_argument('-r','--reference', help='Reference fasta file with with the alignment was constructed', required=True)
parser.add_argument('-o','--output',help='Output file name', required=True)
parser.add_argument('-d','--debug', help='Turn on debugging', required=False, action="store_true")
args = parser.parse_args()


# Visual representation of how reads match reference
PRINT_COLOURED_DIFF = args.debug

# max number of mutations called in a read
MAX_ERRORS = 4

# How many nucleotides need to match at the end of a read for a valid alignment:
# - matching 2 should correct alignment errors, 3 avoids problems with InDel
#   repositioning
# - 3 also simplifies handling the first codon: it's either complete or it's OK
#   to move 1 or 2 bases over to the next triplet
MATCH_N_END = 3


def findEnds(read, ref):
    """
    Figure out where the first symbol of interest in the read is. This start point should be:
    - Not a -.
    - Not (yet) triplet aligned
    Similarly find the last non-dash symbol in the read.
    """
    ends = {"start": 0, "end": len(ref)}
    for i in range(0, len(ref)):
        if read[i] != "-":
            ends["start"] = i
            ends["aligned"] = i
            break

    for i in range(0, len(ref)):
        if read[-1-i] != "-":
            ends["end"]= len(ref) - i
            break

    while (ends.get("aligned") % 3 != 0):
        ends["aligned"] +=1      

    return ends

    
def endMatch(read, ref, ends):
    """
    Aligner errors arise when mutations are at the ends of the read rather than in the middle.
    Trimming ends only shifts the problem. Instead require that read ends match the reference.
    Return False if either end doesn't match
    """
    startRead = read[ends.get("start"):ends.get("start") + MATCH_N_END]
    startRef = ref[ends.get("start"):ends.get("start") + MATCH_N_END]
    endRead = read[ends.get("end") - MATCH_N_END:ends.get("end")]
    endRef = ref[ends.get("end") - MATCH_N_END:ends.get("end")]
    start = startRead == startRef
    end = endRead == endRef
    return start and end       


def hasMultipleGaps(read):
    """
    Return true iff read contains multiple gaps.
    """
    # This regex looks for 3 distinct islands of letters. That implies at least two gaps.
    return re.search('[ATGC]+-+[ATGC]+-+[ATGC]+', str(read)) != None


def findGap(read):
    """
    Get the start and end index (as a 2-tuple) of the gap in readMatch, or None if there isn't one.
    """
    # Find the gap, captured in group 1. This gives us the indexes in the string where the gap starts and
    # ends. (we define "gap" as "dashes in between letters". We previously used hasMultipleGaps to reject
    # any strings which have multiple gaps (which would confuse this regex)).
    match = re.search('[ATGC]+(-+)[ATGC]+', str(read))
    if match == None:
        # No gap exists.
        return None

    return (match.start(1), match.end(1))


def gapAlign(read, gap):
    """
    Perform the "letter-stealing" operation:
    - Find the gap (if there is one).
    - While the gap startpoint is not aligned to a triplet boundary, move letters from the end of the gap
      to the start of the gap.

    @param readMatch Aligned read to process.
    @param gap The gap, as returned by findGap.
    """
    
    if gap is None:
        return read
        
    movingGap = (gap[0], gap[1])
        
    # Shift letters from the end to the start...
    while movingGap[0] % 3 != 0:
        assert(read[movingGap[0]] == "-")

        # This means we ran out of symbols to steal before we managed to align.
        # This indicates a gap at the end, and one we can't fix. Reject!
        if read[movingGap[1]] == "-":
            return None
        
        read[movingGap[0]] = read[movingGap[1]]
        read[movingGap[1]] = "-"

        # Shift the gap to the right...
        movingGap = (movingGap[0] + 1, movingGap[1] + 1)
      
    return read

def verifyRead(read, ref, rejected):

    # Reject reads with multiple gaps
    if hasMultipleGaps(read):
        rejected['multigap'] +=1
        return

    ends = findEnds(read, ref)

    # Reject reads with errors at the ends
    if not endMatch(read, ref, ends):
        rejected['ends misalign'] += 1
        return

    # Reject reads with gaps that aren't a multiple of 3 in length.
    gap = findGap(read)

    if gap is not None:
        gapSize = (gap[1] - gap[0])
        if gapSize % 3 != 0:
            rejected['gap length ' + str(gapSize)] +=1
            return

    return True


def findErrors(read, ref, rejected):
    """
    We've got to find and report several types of difference: insertion, deletion, and substitution.
    We'll begin by running through the strings until we find the first 3-letter block that contains a "-"
    in either string, or which has letters in both strings, but they differ.
    After that, we will regard any "-" in read as being part of a deletion error, any "-" in ref as
    being part of an insertion error, until we reach the point where all remaining symbols in read are
    # "-" (at which point we are finished).
    """
    
   
    # Perform letter-stealing across the gap (if any). The resulting modified read will be ready for naive
    # triplet-wise comparison.
    
    ends = findEnds(read, ref)
    gap = findGap(read)

    if gap is None:
        errors = ["Sub"]
    else:
        errors = [str(gap[0]+1)]

    read = gapAlign(read, gap)
    
    if read is None:
        rejected['unalignable ends'] +=1
        return
    
    for i in range(ends.get("aligned"), ends.get("end"), 3):
        # Check if we've run out of symbols (the rest of the string is just dashes).
        suffix = str(read[i:])
        if re.search('[ATGC]', suffix) is None:
            break

        # Read this triplet from both strings.
        expectedCodon = getChars(str(ref), i, 3)
        actualCodon = getChars(str(read), i, 3)

        # Check if this is the last acid, and it's incomplete, ignore it.
        suffix = str(read[i + 3:])
        if re.search('[ATGC]', suffix) is None and "-" in actualCodon:
            continue

        # Compare the triplets! Use 1-based counting for biologists.
        if expectedCodon != actualCodon:
            errors.extend([str(i+1), expectedCodon, actualCodon])
              
    return errors


def prepareCounts(reference):
    """
    Generate all possible deletions, clasify them and set up counts
    Mutations in NGS reads are then compared against this list and if they fit,
    they are counted towars the final summary
    1. create a "mutation" read
    2. Convert mutation & reference to suitable RefSeq objects
    3. Call errors
    4. Append to count dictionary
    """
    deletion = [3, 6, 9]
    codons = list(CodonTable.unambiguous_dna_by_name["Standard"].forward_table.keys())
    codons.sort()
    ref = str(reference)
    validDels = []
    validSubs = []
    tmp = defaultdict(int)    # create keys as new mistakes are found
    
    for length in deletion:
        for i in range(len(reference) - length):
            possibleread = MutableSeq(ref[:i] + ("-"*length) + ref[i+length:],
                               reference.alphabet)
            if verifyRead(possibleread, reference, tmp) is None:
                continue   
            errors = findErrors(possibleread, reference, tmp)
            if errors == ["Sub"]:
                continue            
            validDels.append(tuple(errors))
            
    for i in range(0, len(reference) - 3,3):
        for triplet in codons:
            s1 = ("Sub", str(i), getChars(ref,i,3), triplet)
            s2 = ("Sub", str(i), ref[i:i+3], ref[i]+triplet[:2], str(i+3),
                  ref[i+3:i+6], triplet[2:]+ref[i+4:i+6])
            s3 = ("Sub", str(i), ref[i:i+3], ref[i:i+2]+triplet[:1], str(i+3),
                  ref[i+3:i+6], triplet[1:]+ref[i+5])           
            validSubs.extend((s1,s2,s3))
           
            
    # tmp collected some rejected counts, now start fresh
    rejected = defaultdict(int)
    
    for key in tmp:
        rejected[key] = 0
    
    return validDels, validSubs, rejected

def printDiff(errors, read, ref, ctr, PRINT_COLOURED_DIFF):
    print("\n\n#############################################################################")
    print("Read: %i \n" % ctr)
    printErrors(errors, read, ref, PRINT_COLOURED_DIFF)
    
    
def main():
    """
    1. Read reference file
    2. Scan over reference sequence to generate all possible mutations
    3. For each ref & read in multiple alignment:
        - verify the read is good quality
        - call the mutation
        - add to count table
    4. Print counts
    """
    time_0 = time.time()
    
    reference = SeqIO.read(args.reference,'fasta', alphabet = IUPAC.ambiguous_dna).seq.upper()

    validDels, validSubs, rejected = prepareCounts(reference)
    allcounts = defaultdict(int)
    
    time_1 = time.time()
    
    print("Completed setting up counts in", (time_1 - time_0) // 60, "min and", (time_1 - time_0) % 60, "s")
       
    read = None
    ref = None
    
    ctr=0 # counter for DEBUG printing
    
    # reading & looping over read/reference sequence in multiple sequence alignment
    # use AlignIO parser and keep sequence only, allowing it to change (important for gap shifts) 
    for pair in AlignIO.parse(args.alignment,"fasta", alphabet = IUPAC.ambiguous_dna, seq_count = 2):

        ctr += 1
        errors = []
                   
        ref = pair[0].seq.tomutable()
        read = pair[1].seq.tomutable()
       
        
        # is the read broken?
        if verifyRead(read, ref, rejected) == None:
            if PRINT_COLOURED_DIFF:
                printDiff(errors, read, ref, ctr, PRINT_COLOURED_DIFF)
            continue

        errors = findErrors(read, ref, rejected)

        if (len(errors)/3) > MAX_ERRORS:
            rejected['too many errors'] +=1
            if PRINT_COLOURED_DIFF:
                printDiff(errors, read, ref, ctr, PRINT_COLOURED_DIFF)
            
        # count up all errors in one massive dictionary
        elif len(errors) > 0:
            allcounts[tuple(errors)] += 1
         

    
    time_2 = time.time()
    print("Completed processing reads in", (time_2 - time_1) // 60, "min and", (time_2 - time_1) % 60, "s")
        
    # once all counts are done, filter out deletions & substitutions
        
    delcounts = {key:value for (key, value) in allcounts.items()
                if (key in validDels) and value > 0}

    subcounts = {key:value for (key, value) in allcounts.items()
                if (key in validSubs) and value > 0}

    time_3 = time.time()
    print("Extracted interesting mutations in", (time_3 - time_2) // 60, "min and", (time_3 - time_2) % 60, "s")


    # Print everything interesting to output
    
    with open(args.output + ".csv", "w") as output:
        print("Rejected", file = output)
        for (key, value) in rejected.items():
            print(key, value, sep=",", file = output)
        print("Substitutions", file = output)
        for (key, value) in subcounts.items():
            print(' '.join(key), value, sep=",", file = output)
        print("Deletions", file = output)
        for (key, value) in delcounts.items():
            print(' '.join(key), value, sep=",", file = output)

    with open(args.output + ".all.csv", "w") as output:
        for (key, value) in allcounts.items():
            print(' '.join(key), value, sep=",", file = output)
    
    with open(args.output + ".d3transposon.csv", "w") as output:
        for (key, value) in delcounts.items():
            if key.count('---') == 1:
                print(key[0], value, sep=",", file = output)    
      
    time_4 = time.time()
    print("Printed results in", (time_4 - time_3) // 60, "min and", (time_4 - time_3) % 60, "s")    

main()
