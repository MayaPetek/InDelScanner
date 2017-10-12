#!/usr/bin/env bash
# ARGUMENTS: TOTAL_FILE_NAME
total_file="$1"
acids=/home/maya/Install/Acids
set -e
 set up file for collector
if [ -s $total_file ]; then
    truncate -s 0 $total_file
fi
# process individual fractions
$acids/scripts/count_one_fraction.sh ../fastq_reads/eGFPN-R1.fastq.gz ../fastq_reads/eGFPN-R2.fastq.gz ../references/eGFP.fa eN N $total_file
$acids/scripts/count_one_fraction.sh ../fastq_reads/eGFPH-R1.fastq.gz ../fastq_reads/eGFPH-R2.fastq.gz ../references/eGFP.fa eH H $total_file
$acids/scripts/count_one_fraction.sh ../fastq_reads/eGFPMM-R1.fastq.gz ../fastq_reads/eGFPMM-R2.fastq.gz ../references/eGFP.fa eMM MM $total_file
$acids/scripts/count_one_fraction.sh ../fastq_reads/eGFPL-R1.fastq.gz ../fastq_reads/eGFPL-R2.fastq.gz ../references/eGFP.fa eL L $total_file
# collect everything together
$acids/scripts/multiple_activity_fractions.py -r ../references/eGFP.fa -f $total_file -e ../Stephane_eGFP.csv
