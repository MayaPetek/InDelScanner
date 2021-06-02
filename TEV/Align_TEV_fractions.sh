#!/bin/bash

set -e

source_dir='$1'
ref_loc_fw='/home/mp/InDelScanner/TEV/wtTEV.fa'
ref_loc_rv=''

mkdir $1/alignments

for filename in target_dir/*fw*.fa; do
    [ -e "$filename" ] || continue
    # get the basename
    frac_name=${filanme%.fa}
    # set up the alignment
    opt/emboss/bin/needleall -gapopen 15 -gapextend 0.5 -asequence $ref_loc_fw -supper1 \
-bsequence filename -supper2 -aformat3 fasta -outfile $fracname.aln -errfile needle.err

done