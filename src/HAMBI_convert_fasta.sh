#!/bin/bash

# Set barcode variable
barcode=$1

awk '/^S/{print ">"$2;print $3}' HAMBI_data/metagenomic_assembly/${barcode}".asm.bp.p_ctg.gfa" > HAMBI_data/metagenomic_assembly/${barcode}".fasta"

# Rename to have the sample name in fron of contig ID
sed "s/^>/&${barcode}_/g" HAMBI_data/metagenomic_assembly/${barcode}".fasta" > HAMBI_data/metagenomic_assembly/${barcode}"_contigs.fasta"