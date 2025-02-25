# Method validation
## Metagenome assemmbly of wastewater communities
```
# Set variable
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p doc/sample_names.txt)

# Run
hifiasm_meta -t 4 -o assemblies/${sample} data/${sample}.hifi_reads.fastq.gz

# Convert to fasta
awk '/^S/{print ">"$2;print $3}' ${sample}.p_ctg.gfa > ${sample}.p_ctg.fa
```

## Methylation analysis of wastewater communities


## UMAP


## Random Forest Classifier


## Sequence logos vs. MultiMotifMaker


## MultiMotifMaker of clustered contigs