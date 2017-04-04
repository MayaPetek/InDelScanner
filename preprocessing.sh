#!/bin/bash

sam="$1"
reference=`sed -n '2p' $2`

# Takes properly aligned matches and extracts reads, throws away anything containing N
cat "$sam" | grep -E "^\S+[	](99|147|83|163)" | cut -d '	' -f 10 | grep -v N > $sam.reads

# Removes reads that fully match refence and keeps the interesting ones; reference is GFP8

while read i;  do
    grep --quiet --color=always $i <<<"$reference"
if [[ $? -eq 0 ]]; then
    true;
else
   echo -e ">read\n$i" >> $sam.unmatchedReads.fa;
fi

done < $sam.reads

# next need to do alignment with needle-all

needleall -gapopen 10 -gapextend 0.5 -asequence $2 -bsequence $sam.unmatchedReads.fa -aformat3 fasta -outfile $sam.aln.tmp

# it names all reference reads with the name of reference file, but compareAcids requires them to be called >Reference

sed s/GFP8/Reference/ $sam.aln.tmp > $sam.aln

# cleanup

rm $sam.aln.tmp $sam.reads
