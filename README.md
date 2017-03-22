# Proteins 'n' stuff

A repo of miscellaneous scripts to help solve a very specific string problem for some reason.

### Data preparation steps

Throw away all the stuff we don't like
This yields only "good enough" matches, apparently. All seem to have been already rejigged to be
forward-only.
99, 83 is fwd
147, 163 is backward
```bash
cat sam | grep -E "^\S+[   ](99|147|83|163)" > filteredSam
```

Extract reads (Just takes the 10th column of the file)
```bash
cat filteredSam | cut -d '  ' -f 10 > reads
```

Replace "N" with ".", since "." is regex for "anything". This causes those symbols to match any symbol
in the result.
```bash
cat reads | sed -e 's/N/./g' > nreads
``

If you want to exclude things with any Ns instead, swap in this command:
```bash
# cat reads | grep -v N > reads
```

Extract list of imperfect matches. You'll want to paste this into a shell script file and invoke it, because
of the heredocument used to speed it up (which isn't allowed inline).
```bash
while read i;  do
    grep --color=always $i <<EOF
ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCCGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTCTGACGTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCACATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCCTTCAAAGATGACGGGACCTACAAGACGCGTGCCGAGGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGTCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAAAAGCTT
EOF
if [[ $? -eq 0 ]]; then
    true;
else
   echo $i >> unmatchedReads;
fi

done < nreads

# TODO: Er, something about the dots?! :P
```


The above takes a while. If you've already done this _not_ excluding "N"-containing reads, you can prune any entries in
unmatchedReads that were casued by "N"-containing reads now by doing this:
```bash
cat unmatchedReads | grep -v "\." > prunedUnmatchedReads
```
(and then maybe overwrite the original file).


Now we want to perform pairwise sequence alignment between every line of `unmatchedReads` and the reference string.
```bash
for i in $(cat unmatchedReadsWithoutN); do ./pairwiseCheck.sh $i >> alignments; done;
```

Then we feed the alignments into the Python script that puts it all together and prints the summary
You should look at compareAcids.py: it has other notes in it.

```
cs compareAcids
./compareAcids.py ../alignments
```
