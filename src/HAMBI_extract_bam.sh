#!/bin/bash

# Set barcode variable
barcode=$1

# Set contig variable
contigs=$(ls HAMBI_data/contigs/*fasta | sed 's/contigs\///g' | sed 's/\.fasta//')

# Load the SeqKit
module load seqkit/2.5.1

# Run
for c in $contigs; do
    echo $c > HAMBI_data/metagenomic_assembly/search
    cat HAMBI_data/metagenomic_assembly/search | xargs samtools view -F 4 -b HAMBI_data/metagenomic_assembly/$barcode"_mapped.bam" > HAMBI_data/contigs/${barcode}"_"${c}".bam"
done