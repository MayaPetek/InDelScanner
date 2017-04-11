#!/bin/bash
# single codon substitutions have 4 fields; field 4 is letters
# two codon substittions have 8 fields; fields 4 & 8 are letters

input=$1
output=$2

# count occurence of mutations in raw mutation file

cat $input | sort -n | uniq -c - $output.count

# split according to type of mutation: 3 bp substitution, 3/6/9 bp deletion, and on codon boundary or not

awk '
BEGIN { FS=" "; OFS=" ";}
    { if (NF==5 && $3 !~ /-/ && $5 !~ /-/ )  print $2, $5, $1 > "'$output'.mutations" }
    { if (NF==9 && $6-$2==3 && $3 !~ /-/ && $5 !~ /-/ && $7 !~ /-/ && $9 !~ /-/ ) print $2, $5, $6, $9, $1> "'$output'.mutations" }

    { if (NF==5 && $3 !~ /-/ && $5 ~ /---/ )  print $2, "del3", $1 > "'$output'.mutations" }
    { if (NF==9 && $6-$2==3 && $3 !~ /-/ && $5 !~ /-/ && $7 !~ /-/ && $9 ~ /---/ ) print $2, $5, "del3", $1 > "'$output'.mutations" }

    { if (NF==9 && $6-$2==3 && $3 !~ /-/ && $5 ~ /---/ && $7 !~ /-/ && $9 ~ /---/ )  print $2, "del6", $1 > "'$output'.mutations" }
    { if (NF==13 && $6-$2==3 && $10-$6==3 && $3 !~ /-/ && $5 !~ /-/ && $7 !~ /-/ && $9 ~ /---/ && $11 !~ /-/ && $13 ~ /---/ ) print $2, $5, "del6", $1 > "'$output'.mutations" }

    { if (NF==13 && $6-$2==3 && $10-$6==3 && $3 !~ /-/ && $5 ~ /---/ && $7 !~ /-/ && $9 ~ /---/ && $11 !~ /-/ && $13 ~ /---/   )  print $2, "del9", $1 > "'$output'.mutations" }
    { if (NF==17 && $6-$2==3 && $10-$6==3 && $14-$10==3 && $3 !~ /-/ && $5 !~ /-/ && $7 !~ /-/ && $9 ~ /---/ && $11 !~ /-/ && $13 ~ /---/ && $15 !~ /-/ && $17 ~ /---/ ) print $2, $5, "del9", $1 > "'$output'.mutations" }
' $output.count

rm $output.count
