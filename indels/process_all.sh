#!/usr/bin/env bash
acids=/home/maya/Install/Acids/indels
set -e

# process individual fractions
$acids/scripts/count_one_fraction.sh ../fastq_reads/S6-d1-R1.fastq.gz ../fastq_reads/S6-d1-R2.fastq.gz ../references/full_fragment.fa S6 d3
$acids/scripts/count_one_fraction.sh ../fastq_reads/S6-d2-R1.fastq.gz ../fastq_reads/S6-d2-R2.fastq.gz ../references/full_fragment.fa S6 d6
$acids/scripts/count_one_fraction.sh ../fastq_reads/S6-d3-R1.fastq.gz ../fastq_reads/S6-d3-R2.fastq.gz ../references/full_fragment.fa S6 d9

$acids/scripts/count_one_fraction.sh ../fastq_reads/S6-i1-R1.fastq.gz ../fastq_reads/S6-i1-R2.fastq.gz ../references/full_fragment.fa S6 i3
$acids/scripts/count_one_fraction.sh ../fastq_reads/S6-i2-R1.fastq.gz ../fastq_reads/S6-i2-R2.fastq.gz ../references/full_fragment.fa S6 i6
$acids/scripts/count_one_fraction.sh ../fastq_reads/S6-i3-R1.fastq.gz ../fastq_reads/S6-i3-R2.fastq.gz ../references/full_fragment.fa S6 i9

#$acids/scripts/one_fraction_local.sh ../fastq_reads/S6-s1-R1.fastq.gz ../fastq_reads/S6-s1-R2.fastq.gz ../references/S6short.fa S6-s1 s1 $total_file
#
#$acids/scripts/one_fraction_local.sh ../../fastq_reads/eGFPL-R1.fastq.gz ../../fastq_reads/eGFPL-R2.fastq.gz ../../references/eGFP.fa eGFP L
#$acids/scripts/one_fraction_local.sh ../../fastq_reads/eGFPH-R1.fastq.gz ../../fastq_reads/eGFPH-R2.fastq.gz ../../references/eGFP.fa eGFP H
#$acids/scripts/one_fraction_local.sh ../../fastq_reads/eGFPMM-R1.fastq.gz ../../fastq_reads/eGFPMM-R2.fastq.gz ../../references/eGFP.fa eGFP MM
#$acids/scripts/one_fraction_local.sh ../../fastq_reads/eGFPN-R1.fastq.gz ../../fastq_reads/eGFPN-R2.fastq.gz ../../references/eGFP.fa eGFP N
##
#$acids/scripts/one_fraction_local.sh ../../fastq_reads/GFP8L-R1.fastq.gz ../../fastq_reads/GFP8L-R2.fastq.gz ../../references/GFP8.fa GFP8 L
##$acids/scripts/one_fraction_local.sh ../../fastq_reads/GFP8H-R1.fastq.gz ../../fastq_reads/GFP8H-R2.fastq.gz ../../references/GFP8.fa GFP8 H
#$acids/scripts/one_fraction_local.sh ../../fastq_reads/GFP8MM-R1.fastq.gz ../../fastq_reads/GFP8MM-R2.fastq.gz ../../references/GFP8.fa GFP8 MM
#$acids/scripts/one_fraction_local.sh ../../fastq_reads/GFP8N-R1.fastq.gz ../../fastq_reads/GFP8N-R2.fastq.gz ../../references/GFP8.fa GFP8 N

# convert mutations
#$acids/scripts/reads_to_mutations.py -f . -r /home/maya/Documents/fitness_landscapes/GFP/references/eGFP.fa
#$acids/scripts/multiple_activity_fractions.py -r ./S6.fa -f $total_file -c

#$acids/scripts/one_fraction_local.sh ../fastq_reads/S6-d1-R1.fastq.gz ../fastq_reads/S6-d1-R2.fastq.gz ../references/out_fragment.fa out_fragment d3 $total_file
#$acids/scripts/one_fraction_local.sh ../fastq_reads/S6-d1-R1.fastq.gz ../fastq_reads/S6-d1-R2.fastq.gz ../references/in_fragment.fa in_fragment d3 $total_file