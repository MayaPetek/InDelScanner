#!/usr/bin/python3

import re
import pickle

from collections import defaultdict

from Bio.Seq import MutableSeq
from Bio.Data import CodonTable

from stringUtils import getChar, getChars

class allowedDict(dict):
    def insert(self, errors):
        # errors is a tuple given by findErrors
        try:
            self[errors] += 1
        except KeyError:
            pass
        return


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
        if read[-1 - i] != "-":
            ends["end"] = len(ref) - i
            break

    while (ends.get("aligned") % 3 != 0):
        ends["aligned"] += 1

    return ends


def endMatch(read, ref, ends, MATCH_N_END):
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

    @param read Aligned read to process.
    @param gap The gap, as returned by findGap.
    """

    if gap is None:
        return read

    movingGap = (gap[0], gap[1])

    # Shift letters from the end to the start...
    while movingGap[0] % 3 != 0:
        assert (read[movingGap[0]] == "-")

        # This means we ran out of symbols to steal before we managed to align.
        # This indicates a gap at the end, and one we can't fix. Reject!
        if read[movingGap[1]] == "-":
            return None

        read[movingGap[0]] = read[movingGap[1]]
        read[movingGap[1]] = "-"

        # Shift the gap to the right...
        movingGap = (movingGap[0] + 1, movingGap[1] + 1)

    return read


def verifyRead(read, ref, rejected, MATCH_N_END):
    # Reject reads with multiple gaps
    if hasMultipleGaps(read):
        rejected['multigap'] += 1
        return

    ends = findEnds(read, ref)

    # Reject reads with errors at the ends
    if not endMatch(read, ref, ends, MATCH_N_END):
        rejected['ends misalign'] += 1
        return

    # Reject reads with gaps that aren't a multiple of 3 in length.
    gap = findGap(read)

    if gap is not None:
        gapSize = (gap[1] - gap[0])
        if gapSize % 3 != 0:
            rejected['gap length ' + str(gapSize)] += 1
            return

    return True


def classifyPoint(expectedCodon, actualCodon, codons):
    # Compare two codons and return 's' if both are letters, 'd' if actual is a deletion and 'i' if insertion
    valid_codons = set(codons + ['---'])
    if {expectedCodon,actualCodon}.issubset(valid_codons):  # both codons are valid
        if actualCodon == "---":
            return 'd'
        elif expectedCodon == "---":
            return 'i'
        else:
            return 's'
    else:
            return 'b'


def findErrors(read, ref, rejected, codons, MAX_ERROR_INDEX):
    """
    @ read, ref: MutableSeq objects
    @ rejected: defaultdict(int) for counting bad reads
    @ codons: all valid codons incl. stop codons, but not '---'
    Returns None for broken reads and a tuple containing errors otherwise

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

    errors = [""]
    # the first position gives error classification, eg. 'dd' or 'sd'

    read = gapAlign(read, gap)

    if read is None:
        rejected['unalignable ends'] += 1
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
            break

        if i > MAX_ERROR_INDEX:
            if errors[0].find('i') == -1:
                break

        # Compare the triplets! Continue with 0-based counts.
        if expectedCodon != actualCodon:
            t = classifyPoint(expectedCodon, actualCodon, codons)
            errors[0] += t
            errors.extend([str(i), expectedCodon, actualCodon])


    return tuple(errors)

def prepareCounts(reference):
    """
    @param reference: a SeqRecord object

    Generate all possible deletions, clasify them and set up counts
    Mutations in NGS reads are then compared against this list and if they fit,
    they are counted towards the final summary
    1. create a "mutation" read
    3. Call errors
    4. Append to count dictionary
    """
    deletion = [3, 6, 9]
    codons = list(CodonTable.unambiguous_dna_by_name["Standard"].forward_table.keys())
    codons += CodonTable.unambiguous_dna_by_name["Standard"].stop_codons

    interesting = ('s', 'ss', 'd', 'sd', 'dd', 'sdd', 'ddd', 'sddd')
    codons.sort()

    counts = {t:allowedDict() for t in interesting}
    rejected = defaultdict(int)  # create keys as new mistakes are found

    ref = reference.seq.upper()
    # Substitutions use getChars to slice strings, need reference to be a string
    r = str(ref)

    for length in deletion:
        for i in range(len(ref) - length):
            possibleread = MutableSeq(r[:i] + ("-" * length) + r[i + length:],
                                      ref.alphabet)
            if not verifyRead(possibleread, ref, rejected):
                continue
            errors = findErrors(possibleread, ref, rejected, codons)
            try:
                counts[errors[0]][errors] = 0
            except KeyError:
                if errors[0] == '':
                    pass
                else:
                    print(errors)
                    raise

    for t in codons:
        for i in range(len(ref) - 3):
            possibleread = MutableSeq(r[:i] + t + r[i + 3:],
                                      ref.alphabet)
            if not verifyRead(possibleread, ref, rejected):
                continue
            errors = findErrors(possibleread, ref, rejected, codons)
            try:
                counts[errors[0]][errors] = 0
            except KeyError:
                if errors[0] == '':
                    pass
                else:
                    print(errors)
                    raise

    return counts


def saveCounts(reference):
    counts = prepareCounts(reference)
    with open(reference.name + '.counts.p', 'wb') as f:
        pickle.dump(counts, f)

def loadCounts(reference):
    with open(reference.name + '.counts.p', 'rb') as f:
        counts = pickle.load(f)
    return counts