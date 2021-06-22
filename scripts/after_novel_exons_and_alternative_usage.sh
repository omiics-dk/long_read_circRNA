#!/bin/bash

cat *.circRNA_exon_usage.txt | sort -k 5,5 | sort -k 1,1 | uniq > circRNA_exon_usage.txt
# Filter to only keep rows with 5 columns. These are the real hits:
cat circRNA_exon_usage.txt | awk 'NF==5{print}{}' > circRNA_exon_usage_filter.txt

# Get the alternatively used exons
cat circRNA_exon_usage.txt | awk '$2>9' | awk '$4<0.9' | awk '$4>0.1' > circRNA_alternative_exon_usage.txt
wc -l circRNA_exon_usage.txt circRNA_alternative_exon_usage.txt

mkdir DELETE
mv core_circ_grep_10reads.sh-*.out DELETE
mv *.circRNA_exon_usage.txt DELETE


