#!/bin/bash

# Load tools
module load seqkit/2.5.1

accessions=$(less ../WW_data/erm_F_UMAP/updated_contig_names.txt)

for a in $accessions;
do
        line=$(grep "$a" ../WW_data/erm_F_UMAP/erm_F_location_flanking.txt)
        if [[ -n "$line" ]]; then
                startFlank=$(echo $line | cut -d' ' -f 4)
                endFlank=$(echo $line | cut -d' ' -f 5)
                seqkit subseq -r $startFlank:$endFlank ../WW_data/erm_F_UMAP/*$a"_bakta_out"/*$a".fna" > ../WW_data/erm_F_UMAP/extracted_data/$a".fasta"
        fi
done