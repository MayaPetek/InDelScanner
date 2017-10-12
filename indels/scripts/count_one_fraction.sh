#!/bin/bash
#========
# USAGE: count_one_fraction.sh forward.fastq.gz reverse.fastq.gz reference.fa output_name activity totalFile
#========

forwardReads="$1"
reverseReads="$2"
referenceName="$3"
referenceFileBase=${3%.fa}
referenceSequence=$(sed -n '2p' $3)
baseName="$4"
activity="$5"
totalFile="$6"

set -e

export PATH=$PATH:/home/maya/Install/Acids/indels

# Check if an alignment already exists
if [ -s $baseName.aln ]; then
    echo -e "\nAn alignment with this output name already exists. Calling mutations now!"
    /home/maya/Install/Acids/scripts/count_one_lane.py -a $baseName.aln -r $referenceFileBase.fa -o $baseName
    echo -e "$activity,$baseName.valid_counts.p,$baseName.assembled.depth.txt,$baseName.unassembled.depth.txt" >> $totalFile
    exit

else
    echo -e "\nStarting with raw reads..."
    # Assemble reads that overlap in the middle, removing reads containing N
    # most N-containing reads have multiple poor quality calls
    /opt/pear-0.9.10-bin-64/pear-0.9.10-bin-64 -f $forwardReads -r $reverseReads -o $baseName \
    --keep-original --min-overlap 5 --min-assembly-length 0 --quality-threshold 15 --max-uncalled-base 0.01
    # gives output in the form:
    # - $basename.assembled.fastq
    # - $basename.unassembled.forward.fastq and $basename.unassembled.reverse.fastq
    echo -e "Read assembly complete"

    # create two SAM files and extract properly mapped reads that don't contain N
    echo -e "\nMapping valid reads to reference..."
    echo -e "Mapping unassembled reads\n"
    bowtie2 -x $referenceFileBase -1 $baseName.unassembled.forward.fastq -2 $baseName.unassembled.reverse.fastq -S $baseName.unassembled.sam
    echo -e "\nMapping assembled reads"
    bowtie2 -x $referenceFileBase -U $baseName.assembled.fastq -S $baseName.assembled.sam


    # sort the sam files to get depth per position
    echo -e "Calculating coverage per position..."
    samtools sort -o $baseName.unassembled.sorted.sam $baseName.unassembled.sam
    samtools depth -a -m 1000000 $baseName.unassembled.sorted.sam > $baseName.unassembled.depth.txt
    samtools sort -o $baseName.assembled.sorted.sam $baseName.assembled.sam
    samtools depth -a -m 1000000 $baseName.assembled.sorted.sam > $baseName.assembled.depth.txt
    echo -e "Bowtie2 alignment complete"

    # Takes properly aligned matches and extracts reads
    # I want columns 1 (name) and 10 (read)

    echo -e "\nExtracting reads from SAM files..."
    if [ -s $baseName.reads ]
    then
        rm $baseName.reads
    fi
    cat $baseName.assembled.sam | grep -E "^\S+[	](0|16)" | cut -f 1,10  > $baseName.reads
    cat $baseName.unassembled.sam | grep -E "^\S+[	](99|147|83|163)" | cut -f 1,10  >> $baseName.reads

    echo -e "Extracting proper reads complete"

    # Replace all N with . (regex wildcard)
    sed -i -e 's/N/./g' $baseName.reads

    # Take second field and compare it to reference
    # Removes reads that fully match reference and keeps the interesting ones; reference is from second line of fasta file
    # Output is in fasta so needle accepts it as input

    echo -e "\nFiltering out reads that fully match reference..."
    # Check if the output file already exists and remove it
    if [ -s $baseName.interestingReads.fa ]
    then
        rm $baseName.interestingReads.fa
    fi

    while read name seq
    do
        if ! grep --ignore-case --quiet $seq <<<"$referenceSequence"; then
            echo -e ">"$name'\n'$seq >> $baseName.interestingReads.fa;
    #    grep --ignore-case --quiet $seq <<<"$referenceSequence" || exit_code=$?
    #    if (( exit_code == 1 )); then
    #        echo -e ">"$name'\n'$seq >> $baseName.interestingReads.fa;
        fi
    done < $baseName.reads
    echo -e "Interesting reads have been filtered out"

    # next need to do alignment with needle-all
    if [ -s $baseName.aln ]
    then
        rm $baseName.aln
    fi

    needleall -gapopen 15 -gapextend 0.5 -asequence $referenceFileBase.fa -supper1 \
    -bsequence $baseName.interestingReads.fa -supper2 \
    -aformat3 fasta -outfile $baseName.aln -errfile $baseName.err
    echo -e "Multiple to one alignment complete\n"


    echo -e "\nCalling mutations"
    time /home/maya/Install/Acids/scripts/count_one_lane.py -a $baseName.aln -r $referenceFileBase.fa -o $baseName
    echo -e "$activity,$baseName.valid_counts.p,$baseName.assembled.depth.txt,$baseName.unassembled.depth.txt" >> $totalFile
    echo -e "Found the mutations"

    echo -e "\nRemoving intermediate files"
    rm $baseName.discarded.fastq
    rm $baseName.unassembled.forward.fastq $baseName.unassembled.reverse.fastq $baseName.assembled.fastq
    rm $baseName.assembled.sam $baseName.unassembled.sam
    rm $baseName.reads $baseName.err
    echo -e "Done with set $baseName"
fi