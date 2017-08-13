#!/bin/bash

# example input: $acids/preprocessing.sh f.fastq r.fastq ../references/eGFP test
forwardReads="$1"
reverseReads="$2"
referenceFileBase="$3"
referenceSequence=`sed -n '2p' $3.fa`
referenceName=`sed -n '1p' $3.fa`
baseName="$4"

# Assemble reads that overlap in the middle, removing reads containing N
# most N-containing reads have multiple poor quality calls
/opt/pear-0.9.10-bin-64/pear-0.9.10-bin-64 -f $forwardReads -r $reverseReads -o $baseName \
--keep-original --min-overlap 5 --min-assembly-length 0 --quality-threshold 10 --max-uncalled-base 0.01
# gives output in the form:
# - $basename.assembled.fastq
# - $basename.unassembled.forward.fastq and $basename.unassembled.reverse.fastq
rm $baseName.discarded.fastq

echo "Read assembly complete"

# create two SAM files
echo "Creating SAM files..."
bowtie2 -x $referenceFileBase -1 $baseName.unassembled.forward.fastq -2 $baseName.unassembled.reverse.fastq \
-S $baseName.unassembled.sam --gbar 3 --rdg 15,1 --no-unal
bowtie2 -x $referenceFileBase -U $baseName.assembled.fastq -S $baseName.assembled.sam \
--gbar 3 --rdg 15,1 --no-unal
rm $baseName.unassembled.forward.fastq $baseName.unassembled.reverse.fastq $baseName.assembled.fastq

# sort the sam files to get depth per position
samtools sort -o $baseName.unassembled.sorted.sam $baseName.unassembled.sam
samtools depth -a -m 1000000 $baseName.unassembled.sorted.sam > $baseName.unassembled.depth.txt
samtools sort -o $baseName.assembled.sorted.sam $baseName.assembled.sam
samtools depth -a -m 1000000 $baseName.assembled.sorted.sam > $baseName.assembled.depth.txt
echo "Bowtie2 alignment complete"

# Takes properly aligned matches and extracts reads
echo "Extracting reads from SAM files..."
cat $baseName.assembled.sam | grep -E "^\S+[	](0|16)"  > $baseName.reads.fastq
cat $baseName.unassembled.sam | grep -E "^\S+[	](99|147|83|163)"  >> $baseName.reads.fastq
#rm $baseName.assembled.sam $baseName.unassembled.sam
echo "Extracting proper reads complete"

# next need to do alignment with needle-all

needleall -gapopen 15 -gapextend 0.5 -asequence $referenceFileBase.fa -supper1 \
-bsequence $baseName.reads.fastq -supper2 \
-aformat3 fasta -outfile $baseName.aln -errfile $baseName.err
echo "Multiple to one alignment complete"

echo "Calling mutations"
$acids/compareAcids/indels.py -a $baseName.aln -r $referenceFileBase.fa -o $baseName
