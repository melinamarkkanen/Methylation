# Method validation
## Metagenome assembly of wastewater communities by hifiasm_meta (https://github.com/chhylp123/hifiasm):
```
# Run
hifiasm_meta -t 4 -o assemblies/${sample} data/${sample}.hifi_reads.fastq.gz

# Convert to fasta
awk '/^S/{print ">"$2;print $3}' ${sample}.p_ctg.gfa > ${sample}.p_ctg.fa
```








## Jotain
```
cd Methylation/WW_data

# Test preanalysis
module load snakemake/7.17.1

snakemake --profile workflow/profile --use-envmodules \
	--snakefile workflow/Snakefile_WW_preanalysis --use-singularity --singularity-args "--bind /scratch/project_2006608/Methylation_Viikki_HiFi/data/" -np

# seraavaks WW methylation, test it like this but run with the queue system
module load snakemake/7.17.1

snakemake --profile workflow/profile --use-envmodules \
	--snakefile workflow/Snakefile_WW_methylation_analysis --use-singularity -np
```

## *Counter for the number of contigs and .bam*

## *cleaner.sh*







## Methylation analysis of wastewater communities
&nbsp;
&nbsp;
&nbsp;
## UMAP
&nbsp;
&nbsp;
&nbsp;
## Random Forest Classifier
&nbsp;
&nbsp;
&nbsp;
## Sequence logos vs. MultiMotifMaker
&nbsp;
&nbsp;
&nbsp;
## MultiMotifMaker of clustered contigs