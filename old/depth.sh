#!/bin/bash

forwardReads="$1"
reverseReads="$2"
referenceFileBase="$3"
referenceSequence=`sed -n '2p' $3.fa`
referenceName=`sed -n '1p' $3.fa`
baseName="$4"

# Assemble reads that overlap in the middle, removing reads containing N
# most N-containing reads have multiple poor quality calls
/opt/pear-0.9.10-bin-64/pear-0.9.10-bin-64 -f $forwardReads -r $reverseReads -o $baseName --keep-original --min-overlap 5 --min-assembly-length 0 --quality-threshold 10 --max-uncalled-base 0.01
# gives output in the form:
# - $basename.assembled.fastq
# - $basename.unassembled.forward.fastq and $basename.unassembled.reverse.fastq
rm $baseName.discarded.fastq

echo "Read assembly complete"

# create two SAM files and extract properly mapped reads that don't contain N
echo "Creating SAM files..."
bowtie2 -x $referenceFileBase -1 $baseName.unassembled.forward.fastq -2 $baseName.unassembled.reverse.fastq -S $baseName.unassembled.sam
bowtie2 -x $referenceFileBase -U $baseName.assembled.fastq -S $baseName.assembled.sam
rm $baseName.unassembled.forward.fastq $baseName.unassembled.reverse.fastq $baseName.assembled.fastq
echo "Bowtie2 alignment complete"

# sort the sam files to get depth per position
samtools sort -o $baseName.unassembled.sorted.sam $baseName.unassembled.sam
samtools depth -a -m 1000000 $baseName.unassembled.sorted.sam > $baseName.unassembled.depth.txt
samtools sort -o $baseName.assembled.sorted.sam $baseName.assembled.sam
samtools depth -a -m 1000000 $baseName.assembled.sorted.sam > $baseName.assembled.depth.txt
echo "Depth calculated"

rm $baseName.assembled.sam $baseName.unassembled.sam
