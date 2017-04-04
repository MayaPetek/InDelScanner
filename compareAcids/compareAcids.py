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
# Reject if more than one gap.

# If more than this many errors are found, we "reject" a read as being too broken.
# TODO: Most of these are when something like this happens:
# ------------TTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGAC
# CATTATCAACAAAATACTCCAATTGGC---------------------------------
# ------------------------------------------------------------
# -----------------------------------------------------------C
# You can probably tweak the parameters of mafft to stop it from doing this, improving
# your data!

MAX_ERRORS = 10;

# Overwrite this many symbols at the start/end of the read string with "-". This apparently helps work around
# some flavours of aligner unreliability.
IGNORE_FIRST_N_SYMBOLS = 0;
IGNORE_LAST_N_SYMBOLS = 0;

PRINT_COLOURED_DIFF = len(sys.argv) >= 3 and sys.argv[2] == "DEBUG";


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

# TODO: Interpret change types.

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
    - Triplet-aligned
    - Not a -.
    """
    startPoint = 0
    for i in range(0, len(refMatch)):
        if readMatch[i] != "-":
            startPoint = i
            break

    # Zip forwards to the next multiple of 3...
    while (startPoint % 3 != 0):
        startPoint += 1

    return startPoint

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


def summarizeChanges(readMatch, refMatch):
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

    # Reject reads with gaps that aren't a multiple of 3 in length.
    gap = findGap(readMatch);
    if gap != None:
        gapSize = (gap[1] - gap[0])

        if gapSize % 3 != 0:
            print("Rejected (gap of size %i)" % gapSize);
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
    else:
        printErrors(errors, readMatch, refMatch, PRINT_COLOURED_DIFF);


# Analyze the output of `pairwiseCheck.sh`, emitting a handy summary of ammino acid changes.
ctr = 0
with open(sys.argv[1]) as f:
    # Discard the first line of the file.
    f.readline();

    while (True):
        if PRINT_COLOURED_DIFF:
            print("\n\n####################################################################################################")
            print("Read: %i \n" % ctr);

        # Read the read and reference blocks for the next one...
        refMatch = readUntil(f, ">");
        readMatch = readUntil(f, ">Reference");

        assert(len(refMatch) == len(readMatch));

        # Did we run out of file?
        if readMatch == "":
            break;

        summarizeChanges(readMatch, refMatch);

        ctr += 1;
