#!/usr/bin/python3

"""
Usage:
generate_baseline_d369i3.py S6-6.fa
1. Generate all possible 3/6/9 bp deletions and +3, +6 bp NNN insertions within full sequence
(from reference S6-6.fa, which is just the gene plus minimal flanking sequence)
2. Print in fasta format

3. Process with needle to see which mutations are distinct / redundant
4. Run PTE_composition.py -f [location of aln] -s 6 -e 3 -o S6.base_d369i36
"""

import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable


# Demand Python 3.
if sys.version_info[0] < 3:
    print("Python 3 is required, but you are using Python %i.%i.%i") % (
           sys.version_info[0], sys.version_info[1], sys.version_info[2])
    sys.exit(1)

def main():
    """
    Read reference, generate all possible sub/del reads and write them to file
    """
    reference = SeqIO.read(sys.argv[1],'fasta', alphabet = IUPAC.ambiguous_dna)
    sequence = str(reference.seq.upper())
    substitution = [3]
    deletion = [3,6,9]
    insertion = [6]
    codons = list(CodonTable.unambiguous_dna_by_name["Standard"].forward_table.keys())
    codons.sort()

    print(codons)
    
    start_offset = 9 # number of nt before A in ATG
    end_offset = 6 # number of nt to ignore at the end
    
    # GFP parameters: start 9, end 6
    
    # len = 720
    # final index = 719
    # offset = 6
    # final deletion = index 711, 712, 713
    
    # deletion removing first ATG: keep sequence [:9] and [12:]
    # deletion removing last codon: keep sequence [:len_seq - len_deletion - end_offset] and [i+len_deletion:]

    with open(reference.name + ".base_d369" + ".fa", "w") as output:
        for length in deletion:
            for i in range(start_offset, len(sequence) - length - end_offset):
                newseq = sequence[:i] + sequence[i + length:]
                record = SeqRecord(Seq(newseq, reference.seq.alphabet),
                                   id="d" + str(length) + '_' + str(i), description="")
                SeqIO.write(record, output, "fasta")

    print('Finished generating deletions')

    # Now take the same code as for deletions, but insert in substitutions

    with open(reference.name + ".base_s3" + ".fa", "w") as output:
        for i in range(start_offset, len(sequence) - 3 - end_offset):
                for c in codons: 
                    newseq = sequence[:i] + c + sequence[i + 3:]
                    record = SeqRecord(Seq(newseq, reference.seq.alphabet),
                                       id='s3_' + str(i) + '_' + str(c), description="")
                    SeqIO.write(record, output, "fasta")

    print('Finished generating 3 bp substitutions')
                
#     with open(reference.name + ".base_s" + str(length) + ".fa", "w") as output:
#         for i in range(start_offset - 1, len(sequence) - length - end_offset):
#             for triplet_1 in codons:
#                 newseq = sequence[:i] + sequence[i + length:]
#                 record = SeqRecord(Seq(newseq, reference.seq.alphabet),
#                                    id="d" + str(length) + '_' + str(i), description="")
#                 SeqIO.write(record, output, "fasta")

#     with open(reference.name + ".base_i3.fa", "w") as output:
#         for i in range(start_offset - 1, len(sequence) - 3 - end_offset):
#             for triplet_1 in codons:
#                     newseq = sequence[:i] + triplet_1 + sequence[i:]
#                     record = SeqRecord(Seq(newseq, reference.seq.alphabet),
#                                        id="insertion6", description="")
#                     SeqIO.write(record, output, "fasta")

#     with open(reference.name + ".base_i6.fa", "w") as output:
#         for i in range(start_offset - 1, len(sequence) - 6 - end_offset):
#             for triplet_1 in codons:
#                 for triplet_2 in codons:
#                     newseq = sequence[:i] + triplet_1 + triplet_2 + sequence[i:]
#                     record = SeqRecord(Seq(newseq, reference.seq.alphabet),
#                                        id="insertion6", description="")
#                     SeqIO.write(record, output, "fasta")


main()
