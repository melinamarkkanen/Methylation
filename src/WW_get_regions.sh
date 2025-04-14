#!/bin/bash

touch ../WW_data/erm_F_UMAP/startFlank.txt
touch ../WW_data/erm_F_UMAP/endFlank.txt

while read -a line
do
  	contig=${line[0]}
        start=${line[1]}
        end=${line[2]}
        startFlank=$((${line[1]} - 7500))
        endFlank=$((${line[2]} + 7500))
        echo $startFlank >> ../WW_data/erm_F_UMAP/startFlank.txt
        echo $endFlank >> ../WW_data/erm_F_UMAP/endFlank.txt

        # Combine into one file
        paste -d '\t' ../WW_data/erm_F_UMAP/erm_F_location.txt ../WW_data/erm_F_UMAP/startFlank.txt ../WW_data/erm_F_UMAP/endFlank.txt > ../WW_data/erm_F_UMAP/erm_F_location_flanking.txt

done < $1