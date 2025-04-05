#!/bin/bash

sample=$1

cd /scratch/project_2006608/Methylation/WW_data/$sample"_contigs"

contigs=$(find . -type f -name "*gff" | sed 's/_basemods\.gff//g')

find . -type f -name "*_basemods.gff" | sed 's/_basemods\.gff//g' | sed 's/\.\///g' > contig_names.txt
touch mod_counts.txt

for c in $contigs; do
        count=$(tail -n+5 $c"_basemods.gff" | wc -l);
        echo $count >> mod_counts.txt;
done

paste -d "\t" contig_names.txt mod_counts.txt > $sample"_mod_counts.txt"

sed -i '1i contig\tmod_count' $sample"_mod_counts.txt"