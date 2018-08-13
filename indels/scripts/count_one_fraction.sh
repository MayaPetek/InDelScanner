#!/bin/bash
#========
# USAGE: count_one_fraction.sh forward.fastq.gz reverse.fastq.gz reference.fa output_name activity
#========

forwardReads="$1"
reverseReads="$2"
referenceName="$3"
referenceFileBase=${3%.fa}
referenceSequence=$(sed -n '2p' $3)
baseName="$4"
activity="$5"


set -e

export PATH=$PATH:/home/maya/Install/Acids/indels

# Check if an alignment already exists

echo -e "\nStarting with raw reads..."
# Assemble reads that overlap in the middle, removing reads containing N
# most N-containing reads have multiple poor quality calls
/opt/pear-0.9.10-bin-64/pear-0.9.10-bin-64 -f $forwardReads -r $reverseReads -o $baseName.$activity \
--keep-original --min-overlap 5 --min-assembly-length 0 --quality-threshold 15 --max-uncalled-base 0.01
# gives output in the form:
# - $basename.$activity.assembled.fastq
# - $basename.$activity.unassembled.forward.fastq and $basename.unassembled.reverse.fastq
echo -e "Read assembly complete"

# create two SAM files and extract properly mapped reads that don't contain N
echo -e "\nMapping valid reads to reference..."
/opt/bowtie2-2.3.0/bowtie2-build $3 $referenceFileBase
echo -e "Mapping unassembled reads\n"
/opt/bowtie2-2.3.0/bowtie2 -x $referenceFileBase -1 $baseName.$activity.unassembled.forward.fastq -2 $baseName.$activity.unassembled.reverse.fastq -S $baseName.$activity.unassembled.sam
echo -e "\nMapping assembled reads"
/opt/bowtie2-2.3.0/bowtie2 -x $referenceFileBase -U $baseName.$activity.assembled.fastq -S $baseName.$activity.assembled.sam


# sort the sam files to get depth per position
echo -e "Calculating coverage per position..."
samtools sort -o $baseName.$activity.unassembled.sorted.sam $baseName.$activity.unassembled.sam
samtools depth -a -m 1000000 $baseName.$activity.unassembled.sorted.sam > $baseName.$activity.unassembled.depth.txt
samtools sort -o $baseName.$activity.assembled.sorted.sam $baseName.$activity.assembled.sam
samtools depth -a -m 1000000 $baseName.$activity.assembled.sorted.sam > $baseName.$activity.assembled.depth.txt
echo -e "Bowtie2 alignment complete"

# Takes properly aligned matches and extracts reads
# I want columns 1 (name) and 10 (read)

echo -e "\nExtracting reads from SAM files..."
if [ -s $baseName.$activity.reads ]
then
    rm $baseName.$activity.reads
fi
cat $baseName.$activity.assembled.sam | grep -E "^\S+[	](0|16)" | cut -f 1,10  > $baseName.$activity.reads
cat $baseName.$activity.unassembled.sam | grep -E "^\S+[	](99|147|83|163)" | cut -f 1,10  >> $baseName.$activity.reads

echo -e "Extracting proper reads complete"

# Replace all N with . (regex wildcard)
sed -i -e 's/N/./g' $baseName.$activity.reads

# Take second field and compare it to reference
# Removes reads that fully match reference and keeps the interesting ones; reference is from second line of fasta file
# Output is in fasta so needle accepts it as input

echo -e "\nFiltering out reads that fully match reference..."
# Check if the output file already exists and remove it
if [ -s $baseName.$activity.interestingReads.fa ]
then
    rm $baseName.$activity.interestingReads.fa
fi

while read name seq
do
    if ! grep --ignore-case --quiet $seq <<<"$referenceSequence"; then
        echo -e ">"$name'\n'$seq >> $baseName.$activity.interestingReads.fa;
#    grep --ignore-case --quiet $seq <<<"$referenceSequence" || exit_code=$?
#    if (( exit_code == 1 )); then
#        echo -e ">"$name'\n'$seq >> $baseName.interestingReads.fa;
    fi
done < $baseName.$activity.reads
echo -e "Interesting reads have been filtered out"

# next need to do alignment with needle-all
if [ -s $baseName.$activity.aln ]
then
    rm $baseName.$activity.aln
fi

/opt/emboss/bin/needleall -gapopen 15 -gapextend 0.5 -asequence $referenceFileBase.fa -supper1 \
-bsequence $baseName.$activity.interestingReads.fa -supper2 \
-aformat3 fasta -outfile $baseName.$activity.aln -errfile $baseName.$activity.err
echo -e "Multiple to one alignment complete\n"


echo -e "\nRemoving intermediate files"
rm $baseName.$activity.discarded.fastq
rm $baseName.$activity.unassembled.forward.fastq $baseName.$activity.unassembled.reverse.fastq $baseName.$activity.assembled.fastq
rm $baseName.$activity.assembled.sam $baseName.$activity.unassembled.sam
rm $baseName.$activity.reads $baseName.$activity.err
echo -e "Done with set $baseName"
