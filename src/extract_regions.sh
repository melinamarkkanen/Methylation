#!/bin/bash

# Load tools
module load seqkit/2.5.1

accessions=$(less sample_names.txt)

for a in $accessions;
do
        line=$(grep "$a" locations_flanking.txt)
        if [[ -n "$line" ]]; then
                startFlank=$(echo $line | cut -d' ' -f 4)
                endFlank=$(echo $line | cut -d' ' -f 5)
                seqkit subseq -r $startFlank:$endFlank *$a"_bakta_out"/*$a".fna" > extracted_data/$a".fasta"
        fi
done