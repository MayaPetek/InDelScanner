from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO

constant_seq = Seq("TCAGGAACTGCAGCTAGC", IUPAC.ambiguous_dna)

with open("output_fasta", 'w') as f_out:
    for primer_seq in SeqIO.parse("fwd_primers.fa", 'fasta'):
        new_seq = constant_seq + primer_seq
        print(new_seq)
        SeqIO.write(new_seq, f_out, 'fasta')
