#!/usr/bin/env python3
"""
1. Convert all *.ab1 files into trimmed fastq
2. Run needle to generate alignment
3. Trim so the read does not extend beyond reference
4. Find and print mutations

"""
import argparse
import glob
import os
import csv

import Bio
from Bio import SeqIO, AlignIO
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC
from Bio.Emboss.Applications import NeedleallCommandline

import numpy as np


def trim_fastq_biopython(in_seq, out_file, q_cutoff=50, consec=5, id=None, rc=False):
    """
    Trim a FASTQ file and write out the trimmed sequences as a FASTQ file.

    Only processes the sequence with identifer string rec.  If id
    is None, takes first sequence.
    """
    # Load in sequences using Bio.SeqIO.  We'll keep the result as a dict.
    seqs = Bio.SeqIO.to_dict(in_seq)

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
    """
    Convert one ab1 file to trimmed fastq
    :param fname:
    :param rc:
    :param q_cutoff:
    :return:
    """

    # Get the prefix of the file
    prefix, suffix = os.path.splitext(fname)

    # Convert and trim FASTQ file
    Bio.SeqIO.convert(fname, 'abi', prefix + '.fastq', 'fastq')
    trim_fastq_biopython(prefix + '.fastq', prefix + '_trimmed.fastq', rc=rc, q_cutoff=q_cutoff)

    return prefix + '_trimmed.fastq'


def convert_ab1(ab1_dir, rc=False):
    """
    Convert all *.ab1 files in input directory to trimmed fastq sequences.
    Write result to new *.fq files in 'fastq' directory.
    :param ab1_dir: Path to directory containing ab1 files.
    :return:
    """

    # Get all the .ab1 files in the directory (glob module is convenient!)
    listing = glob.glob(os.path.join(ab1_dir, '*.ab1'))

    # Go through each file in the chromatogram directory and convert it to FASTQ
    for fname in listing:
        # Get the prefix of the file
        prefix = os.path.splitext(os.path.split(fname)[-1])[0]

        # Make the name of the output FASTQ file
        fastq_fname = prefix + '.fastq'

        # Use Biopython to convert file format, then trim the results
        record = SeqIO.read(fname, 'abi')
        Bio.SeqIO.convert(fname, 'abi', fastq_fname, 'fastq')

        trim_fastq_biopython(fname, prefix + '.fq', rc)

        # Remove the temporary full fastq file
        os.remove(fastq_fname)


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


def needle_align(fqname, reffile):
    """
    Use the Emboss Needle package to align fastq read to reference, return trimmed reads from the alignment
    :param fqname: name of fastq file
    :param argsref: name of reference file
    :return:
    """
    prefix, suffix = os.path.splitext(fqname)
    outname = prefix + '.aln'
    needle_cline = NeedleallCommandline(r'/opt/emboss/bin/needleall', bsequence=fqname, asequence=reffile,
                                        verbose=False, gapopen=15, gapextend=0.5,
                                        outfile=outname, aformat='fasta')
    needle_cline()

    # the alignment contains only two sequences, so use Bio.AlignIO.read
    pair = Bio.AlignIO.read(outname, "fasta", alphabet=Bio.Alphabet.IUPAC.ambiguous_dna)
    ref = pair[0].seq.tomutable()
    read = pair[1].seq.tomutable()

    return trim_read(ref, read)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reference', help='FASTA file with reference sequence', required=True)
    parser.add_argument('-d', '--directory', help='Directory containing ab1 files', required=True)
    parser.add_argument('--input', help='CSV file documenting samples', required=True)
    parser.add_argument('-o', '--output', help='Name of output CSV file', required=True)
    args = parser.parse_args()

    convert_ab1(args.directory, rc=True)

    with open(args.output, 'w') as out:
        fieldnames = ['Sample', 'DNA_hgvs', 'Protein_short', 'Protein_tuple']
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        with open(args.input, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row['Use'] != 'no':
                    # split sample name, read in fastq, run needle, find mutations, write output
                    sample, primer, suffix = row['File'].split('.')
                    fqname = '.'.join([sample, primer, 'fq'])
                    ref, read = needle_align(fqname, args.reference)
                    print(ref, read, sep='\n')



