#!/bin/bash

cd /scratch/project_2006608/Methylation/HAMBI_data/contigs

contigs=$(ls *gff | sed 's/_basemods\.gff//g')

ls *_basemods.gff | sed 's/_basemods\.gff//g' > contig_names.txt
touch mod_counts.txt

for c in $contigs; do
        count=$(tail -n+5 $c"_basemods.gff" | wc -l);
        echo $count >> mod_counts.txt;
done

paste -d "\t" contig_names.txt mod_counts.txt > HAMBI_mod_counts.txt