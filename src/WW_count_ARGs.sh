#!/bin/bash

sample=$1

cd /scratch/project_2006608/Methylation/WW_data/$sample"_contigs"

# Filter by alignment length
awk '$4 >= 100' $sample"_resfinder_out.txt" > tmp && mv tmp $sample"_resfinder_out.txt"

contigs=$(find . -type f -name "*gff" | sed 's/_basemods\.gff//g' | sed 's/\.\///g')

touch ARG_counts.txt
touch names.txt

for c in $contigs; do
        echo $c >> names.txt;
        count=$(grep -o "$c" $sample"_resfinder_out.txt" | wc -l);
        echo $count >> ARG_counts.txt;
done

paste -d "\t" names.txt ARG_counts.txt > $sample"_ARG_counts.txt"

# Add headers
sed -i '1i contig\tARG_count' $sample"_ARG_counts.txt"