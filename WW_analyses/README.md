# Analysis steps for the wastewater community data
&nbsp;
## Preparatory analysis
### Metagenomic assembly
```
hifiasm_meta -t 4 -o assemblies/EFF1 EFF1.hifi_reads.fastq.gz

# Convert to fasta
awk '/^S/{print ">"$2;print $3}' EFF1.p_ctg.gfa > EFF1/EFF1_contigs.fasta

# Extract contigs into separate folders
cd Methylation/WW_data/EFF1_contigs
module load biokit
seqretsplit EFF1_contigs.fasta

# For the Snakemake the assembly files must be in specific folders
cd Methylation/WW_data/EFF1
```
### The .bam with kinetic tags were created previously with this [workflow](https://github.com/melinamarkkanen/sul4_project/blob/main/workflow/Snakefile_HiFi_reads)

### Run [Snakefile_WW_preanalysis_analysis](./../workflow/Snakefile_WW_preanalysis_analysis) for:
- aligning HiFi reads **with** kinetics tags (.bam) for the methylation analysis
- preparing the required files for the downstream steps of the methylation analysis
```
snakemake --profile workflow/profile --use-envmodules \
	--snakefile workflow/Snakefile_WW_preanalysis --use-singularity --singularity-args "--bind /scratch/project_2006608/Methylation_Viikki_HiFi/data/" -np
```
## Preliminary analysis for the methylation detection
### Run [Snakefile_WW_methylation_analysis](./../workflow/Snakefile_WW_methylation_analysis) using HyperQueue [script](./../src/sbatch-hq-sing.sh):
- obtain the .gff files for downatream analyses
```
# Dry run
module load snakemake
snakemake --profile workflow/profile --use-envmodules \
	--snakefile workflow/Snakefile_WW_methylation_analysis --use-singularity -np

# Actual run in HyperQueue
sbatch sbatch-hq-sing.sh
```
&nbsp;
## Generate Position Weight Matrices (PWM)
### [scoring_matrices_WW.py](./../src/scoring_matrices_WW.py)
- to filter data, the methylation types that have less than **50** detected sites are filled with 0 matrices which increased the models performance
- the scoring matrices are then flattened to feature matrices. The flattened feature matrices are then used to train the random forest model to predict the taxonomic classification of the contigs. 
```
# Generate the matrices (interactive session)
module load python-data
cd /scratch/project_2006608/Methylation
python3 src/scoring_matrices_WW.py WW_data/EFF1_contigs WW_data/EFF1_matrices_top50
...
```
### Create merged_data.tsv files (separate for each sample)
```
# Combine modification types
cd WW_data/EFF1_matrices/flattened

# make copies of the .tsv files....
cp m6A.tsv cp_m6A.tsv
cp m4C.tsv cp_m4C.tsv
cp modified_base.tsv cp_modified_base.tsv

# Add shared row names
less m4C.tsv | wc -l
yes "EFF1" | head -n 60373 > sample

# Remove extra contig names
cut -f 2-  m6A.tsv > tmp && mv tmp m6A.tsv
cut -f 2-  modified_base.tsv > tmp && mv tmp modified_base.tsv
# Paste
paste m4C.tsv modified_base.tsv m6A.tsv sample > EFF1_concat_matrices.tsv

# Check number of columns
awk -F'\t' '{print NF; exit}' EFF1_concat_matrices.tsv
# 494

# Check row names
#less EFF1_concat_matrices.tsv | cut -f 1 | head

# Add the final column headers
## *This is just given now based on previous analyses, could be a function to make it*

# The Common_id etc starts at 494-
cat ../../header.tsv | cut -f 1-493 > tmp
cat tmp EFF1_concat_matrices.tsv > temp_EFF1_concat_matrices.tsv && mv temp_EFF1_concat_matrices.tsv EFF1_concat_matrices.tsv

# add header 'sample' as the last column name

# Check again
less EFF1_concat_matrices.tsv  | cut -f 1 | head

# Check number of cols
awk -F'\t' '{print NF; exit}' EFF1_concat_matrices.tsv

# Rename according to included modifications
mv EFF1_concat_matrices.tsv EFF1_concat_matrices_top50.tsv
```

## Taxonomical composition
```
# Run Sylph 
apptainer_wrapper exec sylph profile db/gtdb_database.syldb HiFi_fastq/*.hifi_reads.fastq.gz -t $SLURM_CPUS_PER_TASK > Sylph_out/Sylph_WW.tsv

# Run sylph-tax
apptainer_wrapper exec sylph-tax taxprof Sylph_out/Sylph_WW.tsv -t GTDB_r214 -o Sylph_out/prefix_

# Merge samples
apptainer_wrapper exec sylph-tax merge *.sylphmpa --column relative_abundance -o Sylph_WW_merged.txt

# Get genus
awk '$1 ~ "clade_name" || $1 ~ "g__" {print $0}' Sylph_HiFi_merged.txt | grep -v "|t__" | grep -v "|s__"  > Sylph_HiFi_merged_genus_full.txt

# Edit
sed -i 's/HiFi_fastq\///g' Sylph_WW_merged_genus_full.txt
sed -i 's/\.hifi_reads\.fastq\.gz//g' Sylph_WW_merged_genus_full.txt
sed -i 's/--bcAd10[0-9][0-9]T//g' Sylph_WW_merged_genus_full.txt
```