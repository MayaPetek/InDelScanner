## Deletion fitness landscapes

#### Wrapper for multiple files
TODO: Some way of collecting all files from one experiment, automate individual processing, bring back for comparison

#### Preprocessing
- If using paired end reads: merge reads with PEAR. Take the opportunity to filter out very broken data.
- Align all reads against reference.
- While we have SAM files, take the opportunity to calculate depth / position.
- Extract reads that are correctly mapped, keep the name.
- Throw away reads the fully match reference. This is faster than NW alignment for all reads later.
- Feed interesting reads to EMBOSS Needleman-Wunsch aligner. I trust it better than samtools to correctly position InDels.

#### Counting substitutions, deletions and combinations for each sequencing file
- Load a dictionary of all interesting mutations we're considering
- Read a read+reference into a SeqRecord
TODO: Fastq in needle should help this have a sensible name, maybe use it for debugging? Print name of read there
- Various checks that the read is not broken
    Barcodes: As long as the read does not contain insertions, the barcode is ignored and does not contribute to detected mutations
- If the mutation is defined as interesting, figure out what kind it is and add to valid_counts dictionary
- Add the mutation to a dictionary counting everything
- Save both dictionaries for later viewing.
