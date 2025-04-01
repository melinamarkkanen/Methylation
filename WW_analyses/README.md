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
- MISTÄ NE KINETIC READIT TULEE
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
## Create scoring matrices and flattened feature matrices
### Position Weight Matrices (PWM) 
- to filter data, the methylation types that have less than **100** or **200** !! detected sites are filled with 0 matrices which increased the models performance
- the scoring matrices are then flattened to feature matrices. The flattened feature matrices are then used to train the random forest model to predict the taxonomic classification of the contigs. 
```
# Generate the matrices (interactive session)
module load python-data
cd /scratch/project_2006608/Methylation
python3 src/scoring_matrices_WW.py WW_data/EFF1_contigs WW_data/EFF1_matrices_top20
python3 src/scoring_matrices_WW.py WW_data/EFF1_contigs WW_data/EFF1_matrices_top100
python3 src/scoring_matrices_WW.py WW_data/EFF1_contigs WW_data/EFF1_matrices_top200
...
```

## Create merged_data_EFF1.tsv
(## Create merged_data_WW.tsv)
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

*************
# Combine all
cd WW_data
cat EFF1_matrices/flattened/EFF1_concat_matrices.tsv EFF2_matrices/flattened/EFF2_concat_matrices.tsv EFF3_matrices/flattened/EFF3_concat_matrices.tsv INF1_matrices/flattened/INF1_concat_matrices.tsv INF2_matrices/flattened/INF2_concat_matrices.tsv INF3_matrices/flattened/INF3_concat_matrices.tsv SLU1_matrices/flattened/SLU1_concat_matrices.tsv SLU2_matrices/flattened/SLU2_concat_matrices.tsv SLU3_matrices/flattened/SLU3_concat_matrices.tsv > merged_data_WW.tsv
*************

# Check number of columns
awk -F'\t' '{print NF; exit}' EFF1_concat_matrices.tsv
# 494

# Check row names
less EFF1_concat_matrices.tsv | cut -f 1 | head

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
mv EFF1_concat_matrices.tsv EFF1_concat_matrices_top100.tsv
```

## Analyse by sample in Jupiter
## Extract clusters
### eg.: ./WW_get_contigs_for_MAGs.sh <sample> <cluster> <above>
./WW_get_contigs_for_MAGs.sh EFF1 C1 100