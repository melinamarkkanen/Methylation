#!/bin/bash

# Load tools
module load seqkit/2.5.1

# Go to folder
cd /scratch/project_2006608/Methylation

# Define variables and point data
sample=$1
cluster=$2
lista=$(cat notebooks/UMAP_WW/$sample"_"$cluster"_contigs_above100.txt")

# Run
mkdir -p WW_data/MAGs/$sample/$cluster/above100_contigs/
for i in $lista;do
        grep -A 1 -f <(echo "$i") WW_data/$sample/$sample"_contigs.fasta" > WW_data/MAGs/$sample/$cluster/above100_contigs/$i".fasta"
done

# Check lengths
seqkit fx2tab --length --name --header-line WW_data/MAGs/$sample/$cluster/*_contigs/*.fasta >> WW_data/MAGs/$sample/$cluster/$cluster"_lengths.txt"