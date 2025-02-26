# Method validation
## Metagenome assembly of wastewater communities by hifiasm_meta (https://github.com/chhylp123/hifiasm):
```
# Run
hifiasm_meta -t 4 -o assemblies/EFF1 EFF1.hifi_reads.fastq.gz

# Convert to fasta
awk '/^S/{print ">"$2;print $3}' EFF1.p_ctg.gfa > EFF1/EFF1_contigs.fasta
```
&nbsp;
&nbsp;
## Methylation analysis
### Preparations
```
# Extract contigs into separate folders
cd Methylation/WW_data/EFF1_contigs
module load biokit
seqretsplit EFF1_contigs.fasta

# For the Snakemake the assembly files must be in specific folders
cd Methylation/WW_data/EFF1
```
### Running ```workflow/Snakefile_WW_preanalysis```: *(should we add snakemake dry run draws?)*
- aligning HiFi reads **with** kinetics tags (.bam) for the methylation analysis
- preparing the required files for the downstream steps of the methylation analysis
```
snakemake --profile workflow/profile --use-envmodules \
	--snakefile workflow/Snakefile_WW_preanalysis --use-singularity --singularity-args "--bind /scratch/project_2006608/Methylation_Viikki_HiFi/data/" -np
```
### Running ```workflow/Snakefile_WW_methylation_analysis``` using HyperQueue: *(should we add snakemake dry run draws?)*
```
# Dry run
module load snakemake
snakemake --profile workflow/profile --use-envmodules \
	--snakefile workflow/Snakefile_WW_methylation_analysis --use-singularity -np

# Actual run in HyperQueue
sbatch sbatch-hq-sing.sh

# The HyperQueue creates a ton of job report files, remove them in screen by:
./cleaner.sh

# Count the obtained .gff files
cd /scratch/project_2006608/Methylation/WW_data/EFF1_contigs
find . -type f -name "*gff" | wc -l
find . -type f -name "*csv" | wc -l
# Compare to the number of contigs
less EFF1_contigs.fasta | grep ">" | wc -l
```
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