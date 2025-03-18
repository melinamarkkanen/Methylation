#!/bin/bash

# Load tools
module load seqkit/2.5.1

# Go to folder
cd /scratch/project_2006608/Methylation

# Define variables and point data
cluster=$1
samples=("bcAd1023T--bcAd1023T" "bcAd1037T--bcAd1037T" "bcAd1039T--bcAd1039T" "bcAd1046T--bcAd1046T" "bcAd1063T--bcAd1063T")
lista=$(cat notebooks/UMAP_HAMBI_top150_features/top150_$cluster"_contigs.txt")

# Run
for sample in "${samples[@]}";do
        mkdir -p HAMBI_data/MAGs/$cluster/$sample"_contigs"/
        for i in $lista;do
                grep -A 1 -f <(echo "$i") HAMBI_data/metagenomic_assembly/$sample"_contigs.fasta" > HAMBI_data/MAGs/$cluster/$sample"_contigs"/$i".fasta"
        done

        # Remove empty ones
        find HAMBI_data/MAGs/$cluster/$sample"_contigs"/ -size 0 -delete
done

# Check lengths
seqkit fx2tab --length --name --header-line HAMBI_data/MAGs/$cluster/*_contigs/*.fasta >> HAMBI_data/MAGs/$cluster/$cluster"_lengths.txt"