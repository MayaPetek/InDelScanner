"""
All functions for finding substitutions, insertions & deletions
"""

import re
import glob
import os
import csv
import pickle
import numpy as np

from Bio.Emboss.Applications import NeedleallCommandline

from collections import defaultdict

import Bio.AlignIO
import Bio.SeqIO
import Bio.Data.CodonTable
import Bio.Alphabet.IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq, Seq


"""
STRING UTILS
"""


def get_char(target, index):
    """
    Get a char from a string, or "-" if `index` is out of range.
    """
    if index < 0 or index >= len(target):
        return "-"

    return target[index]


def get_chars(target, index, n):
    """
    Get a sequence of `n` chars from `target` at `index`, safely, via `getChar`
    """

    out = ""
    for i in range(0, n):
        out += get_char(target, index + i)

    return out


def set_char(target_str, char, i):
    """
    Set the character at index `i` in `targetStr` to `char`, returning the new string.
    """

    return target_str[:i] + char + target_str[i + 1:]


"""
INPUT
Read in data, prepare dictionaries, get experimental scores
"""


def get_codons(withdeletions=True):
    """
    :type withdeletions: bool
    """
    codons = list(Bio.Data.CodonTable.unambiguous_dna_by_name["Standard"].forward_table.keys())
    codons += Bio.Data.CodonTable.unambiguous_dna_by_name["Standard"].stop_codons
    if withdeletions:
        codons = codons + ['---']
    codons.sort()

    return codons


def prepare_counts(reference):
    """
    @param reference: a SeqRecord object

    Generate all possible deletions, clasify them and set up counts
    Mutations in NGS reads are then compared against this list and if they fit,
    they are counted towards the final summary
    1. Create a "mutation" read
    3. Find the mutation
    4. Append to count dictionary
    """
    deletion = [3, 6, 9]
    codons = get_codons(withdeletions=False)

    interesting = ('s', 'ss', 'd', 'sd', 'dd', 'sdd', 'ddd', 'sddd')
    counts = {}

    rejected = defaultdict(int)  # create keys as new mistakes are found

    ref = reference.seq.upper()
    # Substitutions use getChars to slice strings, need reference to be a string
    r = str(ref)

    # MAKE DELETIONS
    for length in deletion:
        for i in range(len(ref) - length):
            possibleread = MutableSeq(r[:i] + ("-" * length) + r[i + length:],
                                      ref.alphabet)
            if not verifyRead(possibleread, ref, rejected, MATCH_N_END=3):
                continue
            errors = find_dna_mutations(possibleread, ref, rejected)

            # add mutation to simplified counts dictionary
            t = classify_mutation(errors, style='dna')
            if t in interesting:
                counts[errors] = 0
            else:
                print(t, errors)
    print("Finished making deletions")

    # MAKE SUBSTITUTIONS
    for c in codons:
        for i in range(len(ref) - 3):
            possibleread = MutableSeq(r[:i] + c + r[i + 3:], ref.alphabet)
            if not verifyRead(possibleread, ref, rejected, MATCH_N_END=3):
                continue
            errors = find_dna_mutations(possibleread, ref, rejected)

            t = classify_mutation(errors, style='dna')
            if t in interesting:
                counts[errors] = 0
            else:
                print(t, errors)
    print("Finished making substutions")

    return counts


def get_counts(reference, suffix):
    """
    Prepare a dictionary collecting all interesting counts
    :param reference: Bio.Seq object
    :param suffix: string, file ending name for the counts file
    :return: dictionary {errors: 0} for all interesting muttions
    """
    try:
        with open(reference.name + suffix, 'rb') as f:
            valid_counts = pickle.load(f)
        print("Imported counts")
    except IOError:
        print("Making counts")
        valid_counts = prepare_counts(reference)
        with open(reference.name + suffix, 'wb') as f:
            pickle.dump(valid_counts, f)
        print("Finished making counts")

    return valid_counts


def set_up_total(reference):
    """
    Find the pre-prepared counts for gene reference, load it and prepare a dictionary
    :param reference: fasta file. Expect a corresponding reference.counts.p, output of prepareCounts & saveCounts.
    :return: dictionary describes everything know about a mutation
             total['type]['counts'/'depth'//'protein'/'exp_activity'/'pred_activity']
    """
    valid_counts = get_counts(reference, '.counts.p')
    total = {}
    i = 0
    for errors in valid_counts.keys():
        total[errors] = {'counts': {}, 'depth': {}, 'index': i}
        i += 1
    print("Set up total")

    return total


def get_sequencing_data(files, total, experimental):
    """
    The input file in args.file lists where counts & depth data is located:
    line = name,counts file loc,depth file loc, depth 2 file loc, activity (H/N/MM/L)
    Read in all mutations, predict activity and place them in total.
    :param total: at first empty dictionary with keys for all expected mutations
    :return: added all known information about expected mutations
    """
    # input file need: path to counts file, fraction name
    with open(files, 'r') as f:
        for line in  f.readlines():
            # name = H/N/L/MM, denotes activity fraction
            name, count_loc, depth1_loc, depth2_loc = line.rstrip().split(',')
            depth = calculate_depth(depth1_loc, depth2_loc)
            # load individual counts
            with open(count_loc, 'rb') as p:
                counts = pickle.load(p)

            # read all data into t
            i = 0
            for errors in counts.keys():
                if len(errors) == 0:
                    continue
                total[errors]['dna_position'] = int(errors[0])
                total[errors]['depth'][name] = average_depth(errors, depth)
                total[errors]['counts'][name] = counts[errors]
                if counts[errors] != 0:
                    i +=1
                total[errors]['protein_mutation'] = mutation_to_protein_notation(errors)
                total[errors]['dna_type'] = classify_mutation(errors, style='dna')
                total[errors]['protein_type'] = classify_mutation(errors, style='protein')
                if errors in experimental:
                    total[errors]['exp_activity'] = experimental[errors]['activity']
                else:
                    total[errors]['exp_activity'] = None

            print("Found", i, "mutations in", name, "fraction")
    return total


def get_total(files_name, experimental, reference=None):
    """
    Retrieve a 'total' dictionary with all saved data
    :param files_name:
    :return:
    """
    try:
        with open(files_name + '.total_no_prediction.p', 'rb') as f:
            total = pickle.load(f)
        print("Imported collected counts")
    except IOError:
        print("Importing sequencing data")
        total = set_up_total(reference)
        total = get_sequencing_data(files_name, total, experimental)
        with open(files_name + '.total_no_prediction.p', 'wb') as f:
            pickle.dump(total, f)
        print("Finished importing sequencing data")

    return total


"""
INDIVIDUAL READS
Check quality, match ends, find positions of indels
"""


def findEnds(read, ref, start_offset):
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

    while ends.get("aligned") % 3 != start_offset % 3:
        ends["aligned"] += 1

    return ends


def endMatch(read, ref, ends, MATCH_N_END=3):
    """
    Aligner errors arise when mutations are at the ends of the read rather than in the middle.
    Trimming ends only shifts the problem. Instead require that read ends match the reference.
    Return False if either end doesn't match
    How many nucleotides need to match at the end of a read for a valid alignment:
    - matching 2 should correct alignment errors, 3 avoids problems with InDel repositioning
    - 3 simplifies the first codon: it's either complete or it's OK to move 1 or 2 bases over to the next triplet
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


def gapAlign(read, gap, start_offset):
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
    newread = read

    # Shift letters from the end to the start...
    while movingGap[0] % 3 != start_offset % 3:
        assert (newread[movingGap[0]] == "-")

        # This means we ran out of symbols to steal before we managed to align.
        # This indicates a gap at the end, and one we can't fix. Reject!
        if read[movingGap[1]] == "-":
            return None

        newread[movingGap[0]] = newread[movingGap[1]]
        newread[movingGap[1]] = "-"

        # Shift the gap to the right...
        movingGap = (movingGap[0] + 1, movingGap[1] + 1)

    return newread


def verifyRead(read, ref, rejected, MATCH_N_END=3):
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


"""
FINDING ERRORS
"""


def find_dna_mutations(read, ref, rejected, MAX_ERROR_INDEX=720):
    """
    @ read, ref: MutableSeq objects
    @ rejected: defaultdict(int) for counting bad reads
    @ codons: all valid codons incl. stop codons and '---' if looking for deletions
    :return errors - tuple (position, expected triplet, actual triplet, ) / none if broken read

    We've got to find and report several types of difference: insertion, deletion, and substitution.
    We'll begin by running through the strings until we find the first 3-letter block that contains a "-"
    in either string, or which has letters in both strings, but they differ.
    After that, we will regard any "-" in read as being part of a deletion error, any "-" in ref as
    being part of an insertion error, until we reach the point where all remaining symbols in read are
    # "-" (at which point we are finished).
    """

    # Perform letter-stealing across the gap (if any). The resulting modified read will be ready for naive
    # triplet-wise comparison.

    # quality control that there are no mutations at ends of reads
    ends = findEnds(read, ref)
    if not endMatch(read, ref, ends):
        return
    gap = findGap(read)

    errors = []
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
        expectedCodon = get_chars(str(ref), i, 3)
        actualCodon = get_chars(str(read), i, 3)

        # Check if this is the last acid, and it's incomplete, ignore it.
        suffix = str(read[i + 3:])
        if re.search('[ATGC]', suffix) is None and "-" in actualCodon:
            break

        if i > MAX_ERROR_INDEX:
            e = tuple(errors)
            t = classify_mutation(e, style='dna')
            if t.find('i') == -1:
                break

        # Compare the triplets! Continue with 0-based counts.
        if not re.match(actualCodon, expectedCodon):
            errors.extend([str(i), expectedCodon, actualCodon])

    return tuple(errors)


def classify_mutation(errors, style='dna'):
    """
    For each triplet (position, expected, actual) in errors decide whether it's s/d/i/broken
    :param errors: tuple
    :param codons_with_deletions: list of all valid codons including '---'
    :return: 'sd' on dna or protein level
    """
    assert style in ('dna', 'protein')
    e = []
    if style == 'dna':
        for i in range(0, len(errors), 3):
            e.append(classify_point_dna(errors[i+1], errors[i+2]))
    elif style == 'protein':
        for i in range(0, len(errors), 3):
            e.append(classify_point_protein(errors[i+1], errors[i+2]))

    classification = ''.join(e)
    if classification == '':
        return 'wt'
    else:
        return classification


def mutation_to_protein_notation(errors):
    """
    This attempts to translate both codons and compares them to find the effect of mutation on the protein
    :param errors: tuple
    :param codons_with_deletions: all valid codons incl. '---'
    :return: S22A/K23Δ . Δ=amino acid deletion.
    """

    # use the function for translating point mutations
    protein = []
    dna_effect = classify_mutation(errors, style='dna')
    if dna_effect == 'wt':
        return ''
    elif 'b' in dna_effect:
        return 'frameshift'
    elif 'i' in dna_effect:
        return 'insertion'

    for i in range(0, len(errors), 3):
        dna_pos = int(errors[i])
        prot_pos = str(int(dna_pos / 3 + 1))
        prot_ref = translate_point(errors[i + 1])
        prot_mut = translate_point(errors[i + 2])

        if prot_ref == prot_mut:
            continue
        protein.append(prot_ref + prot_pos + prot_mut)

    return '/'.join(protein)


# # From preferences.py
#
# def readDepth(depthfile):
#     with open(depthfile) as f:
#         d = csv.reader(f, delimiter=',')
#         h = next(d)
#         depth = defaultdict(dict)
#         for line in d:
#             depth[int(line[0])] = {h[i]: int(line[i]) for i in [1, 2, 3, 4]}
#     return depth
#
#
# """ Counts files comes in order of: Rejected, Substitutions, Deletions"""
#
#
# def readRejected(name, filelocation, rejected):
#     """ Read in file and add counts to corresponding index"""
#     with open(filelocation, newline='\n') as f:
#         rr = csv.reader(f, delimiter=',')  # gives list of strings
#         for line in rr:
#             if line[0] == "Rejected":
#                 continue
#             elif line[0] == "Substitutions":
#                 break
#             else:
#                 rejected[line[0]][str(name)] = int(line[1])
#     return rejected
#
#
# def readDeletions(name, filelocation, dels):
#     """ Read in file and add counts to corresponding index"""
#     with open(filelocation, newline='\n') as f:
#         rr = csv.reader(f, delimiter=',')  # gives list of strings
#         for line in rr:
#             if "Deletions" in line[0]:
#                 for line in rr:
#                     if line[0] == "Deletions":
#                         continue
#                     else:
#                         dels[line[0]][str(name)] = int(line[1])
#
#     return dels
#
#
# def typesubdel(mut):  # mut is a string
#
#     m = mut.split()
#     pos = int(m.pop(0))
#     name = ' '.join(m)
#
#     subdel = ""
#     for i in range(0, len(m), 3):
#         if m[i + 2] == '---':
#             subdel += "d"
#         else:
#             subdel += "s"
#
#     return pos, name, subdel
#
#
# def adjustFreq(depth, dels, hml):
#     # start with three interesting fractions & empty adjusted frequencies dict
#     freq = defaultdict(dict)
#
#     for mutation, library in dels.items():
#         pos, name, subdel = typesubdel(mutation)
#         if int(library.get("N", 0)) < 4:
#             continue
#         elif depth[pos]["N"] < 10000:
#             continue
#         else:
#             # we have found a well characterized mutation
#             freq[pos]["Name"] = name
#             freq[pos]["Residue"] = int(1 + (int(mutation.split(' ')[1]) - 1) / 3)
#             freq[pos]["SubDel"] = subdel
#             for entry in hml:
#                 freq[pos][entry] = (library.get(entry, 0) * depth[pos][entry]) / (
#                     library.get("baseline") * library.get("N") * depth[pos]["N"])
#     return freq
#
#
# def intensity(freq, hml):
#     for pos, v in freq.items():
#         t = {k: v[k] for k in hml}
#         m = max(t, key=lambda x: t[x])
#         # max_key = max(v, key=lambda k: v[k])
#         if t[m] > 1:
#             freq[pos]["class"] = m
#         else:
#             freq[pos]["class"] = 0
#     return freq
#


"""
POINT FUNCTIONS
Everything to do with a signle amino acid - conversion, classification, comparison
"""

def classify_point_dna(expected, actual):
    """
    Compare two codons and classify whether the mutation is a substition/insertion/deletion on DNA level
    Assumes the two codons are different
    :param expectedCodon: string
    :param actualCodon: string
    :return: one letter string d/i/s/b
    """
    # Compare two codons and return 's' if both are letters, 'd' if actual is a deletion and 'i' if insertion
    cdn = set(get_codons(withdeletions=True))
    if {expected, actual}.issubset(cdn):  # both codons are valid
        if actual == "---":
            return 'd'
        elif expected == "---":
            return 'i'
        elif expected == actual:
            return ''
        else:
            return 's'
    else:
            return 'b'


def translate_point(triplet):
    """
    Take a codon and translate into single letter amino acid, including deletions
    :param triplet: string
    :return: single letter string
    """
    codon = Seq(triplet, alphabet=Bio.Alphabet.IUPAC.ambiguous_dna)
    if triplet == '---':
        aa = 'Δ'
    else:
        aa = str(codon.translate())

    return aa


def classify_point_protein(expected, actual):
    """
    Generate first position of errors, ie. 'sdd'/'d', classification.
    Assume errors is a well-behaved mutation, that is s/d only
    :param expected: dna triplet
    :param actual: dna triplet string
    :return: string, 's' or 'd'
    """
    prot_ref = translate_point(expected)
    prot_mut = translate_point(actual)

    if prot_ref == 'Δ':
        return 'i'
    elif prot_mut == 'Δ':
        return 'd'
    elif prot_ref != prot_mut:
        return 's'
    else:
        return ''

    # dna_ref = Seq(expected, alphabet=Bio.Alphabet.IUPAC.ambiguous_dna)
    # dna_mut = Seq(actual, alphabet=Bio.Alphabet.IUPAC.ambiguous_dna)
    # prot_ref = str(dna_ref.translate())
    #
    # if str(dna_mut) == '---':
    #     return 'd'
    # else:
    #     prot_mut = str(dna_mut.translate())
    #
    # if prot_ref != prot_mut:
    #     return 's'
    # else:
    #     return ''

"""
DEPTH as generated by samtools
"""

def calculate_depth(depth_1, depth_2):
    """
    Collect samtools depth output into a list. 2nd column = 1-based position, 3rd column = coverage.
    Samtools gives two separate files for assembled and unassembled reads
    :param depth_1: output of samtools depth, tab delimited
    :param depth_2: same for other set of reads
    :return: a list of ints with coverage per position
    """

    depth = []
    # depth of assembled and unassembled reads is in two separate files
    with open(depth_1, 'r') as f:
        for line in f.readlines():
            l = line.split()
            # depth in third column, 'samtools depth' output
            depth.append(int(l[2]))
    # open second file and add the count to same position in depth list
    with open(depth_2, 'r') as f:
        for line in f.readlines():
            l = line.split()
            i = int(l[1]) - 1
            depth[i] += int(l[2])
    return depth

def average_depth(errors, depth):
    """
    Use errors description to find which nucleotides it fits in, average coverage over that sequence.
    Includes positions in the mutations codons that don't change themselves.
    :param errors: a tuple describing a set of mutations
    :param depth: a list giving coverage per nucleotide position, eg. depth[16]=3378
    :return: a float giving average coverage for this mutation
    """
    # to allow non-consecutive mutations
    e_depth = []
    for i in range(0, len(errors), 3):
        start = int(errors[i])
        e_depth += depth[start:start+3]
    avg = sum(e_depth) / len(errors)

    return avg

"""
SANGER SEQUENCING
"""

# Useful function to make directories (and pass if dir already exists)
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError:
        if os.path.isdir(path):
            pass
        else: raise


def trim_fastq_biopython(in_file, out_file, q_cutoff=50, consec=6, id=None, rc=False):
    """
    Trim a FASTQ file and write out the trimmed sequences as a FASTQ file.

    Only processes the sequence with identifer string rec.  If id
    is None, takes first sequence.
    """
    # Load in sequences using Bio.SeqIO.  We'll keep the result as a dict.
    with open(in_file, 'rU') as f:
        seqs = Bio.SeqIO.to_dict(Bio.SeqIO.parse(f, 'fastq'))

    # Pull out the id we want
    if id is None:
        key, seq = seqs.popitem()
    else:
        try:
            seq = seqs[id]
        except KeyError:
            raise KeyError('id not found in input file')

    # Get Boolean array for good quality
    q_good = np.array(seq.letter_annotations['phred_quality']) >= q_cutoff

    # Find first set of consec good bases
    i = 0
    while i < len(q_good) - consec and not q_good[i:i + consec].all():
        i += 1

    # Find last set of consec good bases
    j = len(q_good)
    while j >= consec and not q_good[j - consec:j].all():
        j -= 1

    # Write out trimmed sequence
    with open(out_file, 'w') as f:
        if rc:
            # if the sequence is from a reverse primer, reverse-complement first
            rc = seq[i:j].reverse_complement(id=True, name=True, description=True)
            Bio.SeqIO.write(rc, f, 'fastq')
        else:
            Bio.SeqIO.write(seq[i:j], f, 'fastq')


def ab1_to_fastq(fname, rc=False, q_cutoff=50):

    # Get the prefix of the file
    prefix, suffix = os.path.splitext(fname)

    # Convert and trim FASTQ file
    Bio.SeqIO.convert(fname, 'abi', prefix + '.fastq', 'fastq')
    trim_fastq_biopython(prefix + '.fastq', prefix + '_trimmed.fastq', rc=rc, q_cutoff=q_cutoff)

    return prefix + '_trimmed.fastq'


def convert_ab1(ab1_dir, rc=False):
    """
    Convert all *.ab1 files in input directory to trimmed fastq sequences.
    Write result to new fastq files in 'fastq' directory.
    :param ab1_dir: Path to directory containing ab1 files.
    :return:
    """
    fastq_dir = 'fastq'
    mkdir_p(fastq_dir)

    # Get all the .ab1 files in the directory (glob module is convenient!)
    listing = glob.glob(os.path.join(ab1_dir, '*.ab1'))

    # Go through each file in the chromatogram directory and convert it to FASTQ
    for fname in listing:
        # Get the prefix of the file
        prefix = os.path.splitext(os.path.split(fname)[-1])[0]

        # Make the name of the output FASTQ file
        fastq_fname = os.path.join(fastq_dir, prefix + '.fastq')

        # Use Biopython to convert file format
        Bio.SeqIO.convert(fname, 'abi', fastq_fname, 'fastq')

    # Get a list of fastq in the chromatogram directory
    listing = glob.glob(os.path.join(fastq_dir, '*.fastq'))

    # Loop through and trim all FASTQ files
    for fname in listing:
        if 'trimmed' not in fname:
            prefix, suffix = os.path.splitext(fname)
            trim_fastq_biopython(fname, prefix + '_trimmed.fastq', rc=True)


def needle_align(fqname, argsref):
    """
    Use the Emboss Needle package to align fastq read to reference, return trimmed reads from the alignment
    :param fqname: name of fastq file
    :param argsref: name of reference file
    :return:
    """
    prefix, suffix = os.path.splitext(fqname)
    outname = prefix + '.aln'
    needle_cline = NeedleallCommandline(r'/opt/emboss/bin/needleall', bsequence=fqname, asequence=argsref,
                                        verbose=False, gapopen=15, gapextend=0.5,
                                        outfile=outname, aformat='fasta')
    needle_cline()

    # the alignment contains only two sequences, so use Bio.AlignIO.read
    pair = Bio.AlignIO.read(outname, "fasta", alphabet=Bio.Alphabet.IUPAC.ambiguous_dna)
    ref = pair[0].seq.tomutable()
    read = pair[1].seq
    read = MutableSeq(str(read).replace('N', '.'), read.alphabet)
    id = pair[1].id
    ref, read = trim_read(ref, read)

    return ref, read, id


def trim_read(ref, read):
    """
    In Sanger sequencing the read will likely extend beyond the gene reference. Trim both to reference length.
    :param ref: MutableSeq for the Sanger read
    :param read: MutableSeq with the reference sequence
    :return: trimmed ref & read
    """
    start, end = 0, len(ref)

    # trim start
    for i in range(len(ref)):
        if ref[i] != '-':
            start = i
            break
    for i in range(len(ref), 1, -1):
        if ref[i-1] != '-':
            end = i
            break

    return ref[start:end], read[start:end]
