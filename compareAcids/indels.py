#!/usr/bin/python3

import sys, re
from stringUtils import getChar, getChars, setChar, stealLetter
from outputUtils import printErrors

from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, AlignIO
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable

# Demand Python 3.
if sys.version_info[0] < 3:
    print ("Python 3 is required, but you are using Python %i.%i.%i") % (sys.version_info[0], sys.version_info[1], sys.version_info[2])
    sys.exit(1)

"""
Arguments:
1.   Alignment of sequencing reads to a reference sequence in FASTA format. The
     reference comes before read sequence. (>Ref, seq, >Read, seq).
2.   Reference sequence in FASTA format. Same file that was used to create the 
     alignment, in particular the reference name here needs to match reference
     name in the alignment.
(3.) Optional: DEBUG print a coloured representation of mismatches
"""
PRINT_COLOURED_DIFF = len(sys.argv) >= 4 and sys.argv[3] == "DEBUG"

# max number of mutations called in a read
MAX_ERRORS = 4 
# how many nucleotides need to match at the end of a read for a valid alignment
# matching 2 should correct alignment errors, 3 avoids problems with InDel repositioning
# 3 also simplifies handling the first codon: it's either complete or it's OK
# to move 1 or 2 bases over to the next triplet
MATCH_N_END = 3 

def findEnds(read, ref):
    """
    Figure out where the first symbol of interest in the read is. This start point should be:
    - Not a -.
    - Not (yet) triplet aligned
    Similarly find the last non-dash symbol in the read.
    """
    ends = {"start": 0, "end": len(ref.seq)}
    for i in range(0, len(ref.seq)):
        if read.seq[i] != "-":
            ends["start"] = i
            ends["aligned"] = i
            break

    for i in range(0, len(ref.seq)):
        if read.seq[-1-i] != "-":
            ends["end"]= len(ref.seq) - i
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
    startRead = read.seq[ends.get("start"):ends.get("start") + MATCH_N_END]
    startRef = ref.seq[ends.get("start"):ends.get("start") + MATCH_N_END]
    endRead = read.seq[ends.get("end") - MATCH_N_END:ends.get("end")]
    endRef = ref.seq[ends.get("end") - MATCH_N_END:ends.get("end")]
    start = startRead == startRef
    end = endRead == endRef
    return start and end       


def hasMultipleGaps(read):
    """
    Return true iff read contains multiple gaps.
    """
    # This regex looks for 3 distinct islands of letters. That implies at least two gaps.
    return re.search('[ATGC]+-+[ATGC]+-+[ATGC]+', str(read.seq)) != None


def findGap(read):
    """
    Get the start and end index (as a 2-tuple) of the gap in readMatch, or None if there isn't one.
    """
    # Find the gap, captured in group 1. This gives us the indexes in the string where the gap starts and
    # ends. (we define "gap" as "dashes in between letters". We previously used hasMultipleGaps to reject
    # any strings which have multiple gaps (which would confuse this regex)).
    match = re.search('[ATGC]+(-+)[ATGC]+', str(read.seq))
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

    old = read
    readSequence = str(read.seq)

    # Shift letters from the end to the start...
    while gap[0] % 3 != 0:
        assert(readSequence[gap[0]] == "-")

        # This means we ran out of symbols to steal before we managed to align.
        # This indicates a gap at the end, and one we can't fix. Reject!
        if readSequence[gap[1]] == "-":
            return None

        c = readSequence[gap[1]]
        readSequence = setChar(readSequence, c, gap[0])
        readSequence = setChar(readSequence, "-", gap[1])

        # Shift the gap to the right...
        gap = (gap[0] + 1, gap[1] + 1)
    
    read = SeqRecord(Seq(readSequence, IUPAC.ambiguous_dna),id=old.id, name=old.name)
    
    return read

def verifyRead(read, ref):
    # Reject reads with multiple gaps
    if hasMultipleGaps(read):
        print("Rejected (multigap)") # will add a multigap counter later
        return

    ends = findEnds(read, ref)

    # Reject reads with errors at the ends
    if not endMatch(read, ref, ends):
        print("Rejected (ends misalign)") # will add counter later
        return

    # Reject reads with gaps that aren't a multiple of 3 in length.
    gap = findGap(read)
    if gap != None:
        gapSize = (gap[1] - gap[0])

        if gapSize % 3 != 0:
            print("Rejected (gap of size %i)" % gapSize)
            return

        # Perform letter-stealing across the gap (if any). The resulting modified read will be ready for naive
        # triplet-wise comparison.
        read = gapAlign(read, gap)
        if read.seq == None:
            print("Rejected (unalignable gap at end)")
            return

    return read, ends

def findErrors(read, ref, ends):
    """
    We've got to find and report several types of difference: insertion, deletion, and substitution.
    We'll begin by running through the strings until we find the first 3-letter block that contains a "-"
    in either string, or which has letters in both strings, but they differ.
    After that, we will regard any "-" in read as being part of a deletion error, any "-" in ref as
    being part of an insertion error, until we reach the point where all remaining symbols in read are
    # "-" (at which point we are finished).
    """
    errors = []
    
    for i in range(ends.get("aligned"), ends.get("end"), 3):
        # Check if we've run out of symbols (the rest of the string is just dashes).
        suffix = str(read.seq[i:])
        if re.search('[ATGC]', suffix) == None:
            break

        # Read this triplet from both strings.
        expectedCodon = getChars(str(ref.seq), i, 3)
        actualCodon = getChars(str(read.seq), i, 3)

        # Check if this is the last acid, and it's incomplete, ignore it.
        suffix = str(read.seq[i + 3:])
        if re.search('[ATGC]', suffix) == None and "-" in actualCodon:
            continue

        # Compare the triplets!
        if expectedCodon != actualCodon:
            errors.append({
                'position': i,
                'expected': expectedCodon,
                'actual':   actualCodon
            })


    if len(errors) > MAX_ERRORS:
        print("Rejected (errcount)");
    elif len(errors) > 0:
        printErrors(errors, str(read.seq), str(ref.seq), PRINT_COLOURED_DIFF)
    
    return errors


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
    reference = SeqIO.read(sys.argv[2],'fasta', alphabet = IUPAC.ambiguous_dna)

    # reading & looping over read/reference sequence in multiple sequence alignment
    aln = SeqIO.parse(sys.argv[1],"fasta", alphabet = IUPAC.ambiguous_dna)
    
    ctr=1 # counter for DEBUG printing
    
    while True:
        if PRINT_COLOURED_DIFF:
            print("\n\n#############################################################################")
            print("Read: %i \n" % ctr)
        ref = next(aln)
        read = next(aln)
      
        # is the read broken?
        if verifyRead(read, ref) == None:
            continue
            
        read, ends = verifyRead(read, ref)
        errors = findErrors(read, ref, ends)
        
        ctr += 1

main()
