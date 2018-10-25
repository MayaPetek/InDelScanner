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
import numpy as np

import Bio
from Bio import SeqIO, AlignIO
from Bio.Alphabet import IUPAC
from Bio.Emboss.Applications import NeedleallCommandline

from ind import trim_read
from composition import find_dna_hgvs, find_protein_diff


def trim_fastq_biopython(in_seq, out_file, rc, q_cutoff=50, consec=5):
    """
    Trim one FASTQ file and write out the trimmed sequences as a FASTQ file.
    """

    # Get Boolean array for good quality
    q_good = np.array(in_seq.letter_annotations['phred_quality']) >= q_cutoff

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
            rev_seq = in_seq[i:j].reverse_complement(id=True, name=True, description=True)
            Bio.SeqIO.write(rev_seq, f, 'fastq')
        else:
            Bio.SeqIO.write(in_seq[i:j], f, 'fastq')


def convert_ab1(ab1_dir, rc):
    """
    Wrapper function.
    Convert all *.ab1 files in input directory to trimmed fastq sequences.
    Write result to new *.fq files in 'fastq' directory.
    :param ab1_dir: Path to directory containing ab1 files.
    :param rc: True/False - reverse complement the sequence
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
        Bio.SeqIO.convert(fname, 'abi', fastq_fname, 'fastq')
        record = SeqIO.read(fastq_fname, 'fastq')
        trim_fastq_biopython(record, prefix + '.fq', rc)

        # Remove the temporary full fastq file
        os.remove(fastq_fname)


def needle_align(fqname, reffile):
    """
    Use the Emboss Needle package to align fastq read to reference, return trimmed reads from the alignment
    :param fqname: name of fastq file
    :param reffile: name of reference file
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
    parser.add_argument('--input', help='CSV file documenting samples, encoded with Unicode-8', required=True)
    parser.add_argument('-o', '--output', help='Name of output CSV file', required=True)
    parser.add_argument('-s', '--start_offset', help='Number of nt before starting ATG (integer)', required=False,
                        default=3, type=int)
    parser.add_argument('-e', '--end_trail', help='Number of nt after end of gene (integer)', required=False, default=3,
                        type=int)
    parser.add_argument('--debug', help='Turn on debugging', required=False, action="store_true")  # Visual
    args = parser.parse_args()

    # convert_ab1(args.directory, rc=True)

    with open(args.output, 'w') as out:
        fieldnames = ['Sample', 'DNA_hgvs', 'Protein_short', 'Protein_tuple']
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()
        with open(args.input, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row['Use'] != 'no':
                    # split sample name, read in fastq, run needle, find mutations, write output
                    sample, primer, suffix = row['File'].split('.')
                    fqname = '.'.join([sample, primer, 'fq'])
                    try:
                        ref, read = needle_align(fqname, args.reference)
                    except Bio.Application.ApplicationError:
                        print('Sample {0} cannot align.'.format(sample))
                        continue
                    except ValueError:
                        print('Sample {0} is empty.'.format(sample))
                        continue

                    sref, sread = trim_read(ref, read)

                    dna_hgvs = find_dna_hgvs(sread, sref, row['Background'], args.debug, args.start_offset, args.end_trail)
                    prot_errors, prot_short = find_protein_diff(sread, sref, args.debug, args.start_offset,
                                                                args.end_trail)
                    writer.writerow({'Sample': sample, 'DNA_hgvs': dna_hgvs,
                                     'Protein_short': prot_short, 'Protein_tuple': prot_errors})
                    print(sample, prot_short, prot_errors, dna_hgvs)
