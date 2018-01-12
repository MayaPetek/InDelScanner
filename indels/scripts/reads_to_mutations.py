#!/usr/bin/python3

import sys
import pickle
import re
import argparse
from collections import defaultdict

from Bio import SeqIO, AlignIO
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq, translate

from indels.ind import trim_read, findEnds, endMatch, findGap, gapAlign
from indels.output import print_coloured_diff


def indel_len(sequence, start):
    l = 0
    while sequence[start + l] == '-':
        l += 1
    return l

def find_DNA_diff(read, ref):
    """
    @ read, ref: MutableSeq objects
    :return errors - tuple (position, expected triplet, actual triplet, ) / none if broken read

    Letter by letter report mutations in NGS read, all counts 1- based in result (code in 0-count).
    - substitution: 78C = nt 78 in reference is changed to C
    - deletions: 78d6 = 6 nt deleted after 78: 1-78, d6, 85-end
    - insertion: 78iATC = after nt 78 inserted seq. ATC
    """

    if read is None:
        return

    # No gap realignment at this point

    # quality control that there are no mutations at ends of reads
    ends = findEnds(read, ref)
    if not endMatch(read, ref, ends):
        return

    # scan read & reference letter by letter, counting position in reference
    # reads have been trimmed so that reference starts @ 0
    dna_errors = []
    ref_index = ends.get('start')
    i = ends.get('start')

    while i < ends.get('end'):
        # check for differences
        if read[i] == ref[i]:
            ref_index += 1
            i += 1

        elif read[i] == '-':
            # start of a deletion
            l = indel_len(read, i)
            # now we know the length of a deletion, check for frameshifts
            if l % 3 == 0:
                dna_errors += [str(ref_index) + 'd' + str(l)]  # deletion length l starting at ref_index in 0-count
                i += l
                ref_index += l
            else:
                dna_errors += [str(ref_index) + 'f']
                break

        elif ref[i] == '-':
            # start of an insertion
            l = indel_len(ref, i)
            # check for frameshifts
            if l % 3 == 0:
                dna_errors += [str(ref_index) + 'i' + str(read[i:i+l])]
                print(dna_errors)
                i += l
            else:
                dna_errors += [str(ref_index) + 'f']
                break

        else:
            # substitution
            dna_errors += [str(ref_index + 1) + str(read[i])]
            i += 1
            ref_index += 1

    return tuple(dna_errors)


def find_protein_diff(read, ref):

    # quality control
    if read is None:
        return
    ends = findEnds(read, ref)
    if not endMatch(read, ref, ends):
        return

    # scan reference triplet by triplet
    # move letters when encountering an indel
    prot_errors = []
    i = ends.get('aligned')
    ref_index = int(ends.get('aligned') / 3) # reference amino acid index

    while i <= ends.get('end'):
        if read is None:
            break
        ref_codon = ref[i:i+3]
        read_codon = read[i:i+3]

        if '-' in read_codon:  # found a deletion
            # Check if this is the last acid, and it's incomplete, ignore it.
            if re.search('[ATGC]', str(read[i + 3:])) is None:
                break

            if read_codon == '---':  # single codon deletion
                prot_errors += [str(ref_index) + 'd']
                i += 3
                ref_index += 1

            else:  # check it's not a frame shift
                gap = findGap(read[i - 1:])
                l = gap[1] - gap[0]
                if l % 3 != 0:
                    prot_errors.append(str(ref_index) + 'f')
                    break
                # realign gap and repeat loop at same position to compare the codons
                gap = (gap[0] + i - 1, gap[1] + i - 1)
                read = gapAlign(read, gap)
                continue

        elif '-' in ref_codon:  # found an insertion
            gap = findGap(ref[i-1:])
            l = gap[1] - gap[0]
            if l % 3 != 0:
                prot_errors.append(str(ref_index) + 'f')
                break
            if gap[0] == 1:  # insertion after codon
                insertion = read[gap[0] + i - 1:gap[1] + i -1]
                prot_errors.append(str(ref_index) + 'i' + str(translate(insertion)) )
                i += l
                ref_index += 1
            else:  # realign gap and repeat loop at same position to compare the codons
                gap = (gap[0] + i - 1, gap[1] + i - 1)
                ref = gapAlign(ref, gap)
                continue

        elif translate(read_codon) != translate(ref_codon):  # must be a substitution
            prot_errors.append(str(ref_index) + str(translate(read_codon)) )
            i += 3
            ref_index += 1

        else:
            i += 3
            ref_index += 1

    return tuple(prot_errors)


def countOneAlignment(args):
    """
    Don't bother with expected/allowed mutations, just find everything and filter later
    Final format: {DNA error: [(protein error), fraction,
    1. Read reference file
    2. Scan over reference sequence to generate all possible mutations
    3. For each ref & read in multiple alignment:
        - verify the read is good quality
        - call the mutation
        - add to count table
    4. Print counts
    """
    # Demand Python 3.
    if sys.version_info[0] < 3:
        print("Python 3 is required, but you are using Python %i.%i.%i") % (
            sys.version_info[0], sys.version_info[1], sys.version_info[2])
        sys.exit(1)

    reference = SeqIO.read(args.reference, 'fasta', alphabet=IUPAC.ambiguous_dna)

    ctr = 0 # counter for DEBUG printing
    one_lane_counts = defaultdict(int)
    rejected = defaultdict(int)

    
    # reading & looping over read/reference sequence in multiple sequence alignment
    # use AlignIO parser and keep sequence only, allowing it to change (important for gap shifts) 
    for pair in AlignIO.parse(args.alignment, "fasta", alphabet=IUPAC.ambiguous_dna, seq_count=2):

        ctr += 1

        # both read and ref are MutableSeq
        ref = pair[0].seq.tomutable()
        read = pair[1].seq.tomutable()
        read = MutableSeq(str(read).replace('N', '.'), read.alphabet)
        readname = pair[1].id

        # trim sequencing read to reference
        ref, read = trim_read(ref, read)

        dna_errors = find_DNA_diff(read, ref) # errors = a tuple
        prot_errors = find_protein_diff(read, ref)
        print(dna_errors, prot_errors)
        one_lane_counts[dna_errors] += 1

    return one_lane_counts


if __name__ == "__main__":
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
    parser = argparse.ArgumentParser(
        description='Finds all 3 bp substitutions and 3/6/9 bp deletions in a gene')
    parser.add_argument('-a', '--alignment', help='Multiple sequence alignment', required=True)
    parser.add_argument('-r', '--reference', help='Reference fasta file with with the alignment was constructed',
                        required=True)
    parser.add_argument('-o', '--output', help='Output file name', required=False)
    parser.add_argument('-d', '--debug', help='Turn on debugging', required=False, action="store_true") # Visual
    args = parser.parse_args()

    one_lane_counts = countOneAlignment(args)

    # How many nucleotides need to match at the end of a read for a valid alignment:
    # - matching 2 should correct alignment errors, 3 avoids problems with InDel
    #   repositioning
    # - 3 also simplifies handling the first codon: it's either complete or it's OK
    #   to move 1 or 2 bases over to the next triplet

    # valid_counts, bad_counts, rejected = countOneAlignment(MAX_ERRORS=4, MATCH_N_END=3)
    #
    # with open(args.output + '.valid_counts.p', 'wb') as f:
    #     pickle.dump(valid_counts, f)
    # with open(args.output + '.bad_counts.p', 'wb') as f:
    #     pickle.dump(bad_counts, f)

    n = 0
    for error, count in one_lane_counts.items():
        if count > 0:
            n +=1
            print(error, count)
    print('Found', n, 'mutations')