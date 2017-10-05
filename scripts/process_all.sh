#!/usr/bin/env bash
# ARGUMENTS: TOTAL_FILE_NAME
total_file="$1"
set -e
# set up file for collector
if [ -s $total_file ]; then
    truncate -s 0 $total_file
fi
# process individual fractions
# $acids/count_one_fraction.sh ../fastq_reads/eGFPN-R1.fastq.gz ../fastq_reads/eGFPN-R2.fastq.gz ../references/eGFP.fa eN N $total_file
$acids/count_one_fraction.sh ../fastq_reads/GFP8H-R1.fastq.gz ../fastq_reads/GFP8H-R2.fastq.gz ../references/GFP8.fa 8H H $total_file
$acids/count_one_fraction.sh ../fastq_reads/GFP8MM-R1.fastq.gz ../fastq_reads/GFP8MM-R2.fastq.gz ../references/GFP8.fa 8MM MM $total_file
$acids/count_one_fraction.sh ../fastq_reads/GFP8L-R1.fastq.gz ../fastq_reads/GFP8L-R2.fastq.gz ../references/GFP8.fa 8L L $total_file
# collect everything together
# $acids/indels/multiple_activity_fractions.py -r ../references/eGFP.fa -f $total_file -o eGFP.total.p
