#!/bin/bash

# Set barcode variable
sample=$1

# Load the SeqKit
module load seqkit/2.5.1

# Set contig variable

contigs=$(ls WW_data/$sample"_contigs"/*.fasta | sed 's/WW_data\/${sample}"_contigs"//g' | sed 's/\.fasta//')

# Complete contigs
for c in $contigs; do
        echo $c > WW_data/$sample/search
        cat WW_data/$sample/search | xargs samtools view -F 4 -b WW_data/$sample/$sample"_contigs_mapped.bam" > WW_data/$sample"_contigs"/$c".bam"
done
