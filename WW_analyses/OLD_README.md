########## OLD ###############
# Metagenome assembly of wastewater communities by hifiasm_meta (https://github.com/chhylp123/hifiasm):
```
# Run
hifiasm_meta -t 4 -o assemblies/EFF1 EFF1.hifi_reads.fastq.gz

# Convert to fasta
awk '/^S/{print ">"$2;print $3}' EFF1.p_ctg.gfa > EFF1/EFF1_contigs.fasta
```
&nbsp;
&nbsp;
# Preliminary analysis for Methylation detection
```
# Extract contigs into separate folders
cd Methylation/WW_data/EFF1_contigs
module load biokit
seqretsplit EFF1_contigs.fasta

# For the Snakemake the assembly files must be in specific folders
cd Methylation/WW_data/EFF1
```
## Running ```workflow/Snakefile_WW_preanalysis```: *(should we add snakemake dry run draws?)*
- aligning HiFi reads **with** kinetics tags (.bam) for the methylation analysis
- preparing the required files for the downstream steps of the methylation analysis
```
snakemake --profile workflow/profile --use-envmodules \
	--snakefile workflow/Snakefile_WW_preanalysis --use-singularity --singularity-args "--bind /scratch/project_2006608/Methylation_Viikki_HiFi/data/" -np
```
## Running ```workflow/Snakefile_WW_methylation_analysis``` using HyperQueue: *(should we add snakemake dry run draws?)*
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
# Creation of Position Weight Matrices (PWM) (example sample EFF1)
- to filter data, the methylation types that have less than **100**!! (50?) detected sites are filled with 0 matrices which increased the models performance
- the scoring matrices are then flattened to feature matrices. The flattened feature matrices are then used to train the random forest model to predict the taxonomic classification of the contigs. 
```
# Generate the matrices (interactive session)
module load python-data
cd /scratch/project_2006608/Methylation
python3 src/scoring_matrices_WW.py WW_data/EFF1_contigs WW_data/EFF1_matrices_top100
...
```

## Create merged_data_EFF1.tsv
(## Create merged_data_WW.tsv)
```
# Combine modification types
cd WW_data/EFF1_matrices/flattened
cd WW_data/INF1_matrices/flattened
cd WW_data/SLU1_matrices/flattened
cd WW_data/EFF3_matrices/flattened
cd WW_data/EFF2_matrices/flattened
cd WW_data/INF2_matrices/flattened
cd WW_data/INF3_matrices/flattened
cd WW_data/SLU2_matrices/flattened
cd WW_data/SLU3_matrices/flattened

# make copies of the .tsv files....
cp m6A.tsv cp_m6A.tsv
cp m4C.tsv cp_m4C.tsv
cp modified_base.tsv cp_modified_base.tsv

# Add shared row names
less m4C.tsv | wc -l
#yes "EFF1" | head -n 60373 > sample
#yes "INF1" | head -n 64737 > sample
#yes "SLU1" | head -n 74387 > sample
#yes "EFF3" | head -n 54829 > sample
#yes "EFF2" | head -n 59694 > sample
#yes "INF2" | head -n 62223 > sample
#yes "INF3" | head -n 62794 > sample
#yes "SLU2" | head -n 76128 > sample
yes "SLU3" | head -n 73646 > sample

# Remove extra contig names
cut -f 2-  m6A.tsv > tmp && mv tmp m6A.tsv
cut -f 2-  modified_base.tsv > tmp && mv tmp modified_base.tsv
# Paste
#paste m4C.tsv modified_base.tsv m6A.tsv sample > EFF1_concat_matrices.tsv
#paste m4C.tsv modified_base.tsv m6A.tsv sample > INF1_concat_matrices.tsv
#paste m4C.tsv modified_base.tsv m6A.tsv sample > SLU1_concat_matrices.tsv
#paste m4C.tsv modified_base.tsv m6A.tsv sample > EFF3_concat_matrices.tsv
#paste m4C.tsv modified_base.tsv m6A.tsv sample > EFF2_concat_matrices.tsv
#paste m4C.tsv modified_base.tsv m6A.tsv sample > INF2_concat_matrices.tsv
#paste m4C.tsv modified_base.tsv m6A.tsv sample > INF3_concat_matrices.tsv
#paste m4C.tsv modified_base.tsv m6A.tsv sample > SLU2_concat_matrices.tsv
paste m4C.tsv modified_base.tsv m6A.tsv sample > SLU3_concat_matrices.tsv

# Check number of columns
#awk -F'\t' '{print NF; exit}' EFF1_concat_matrices.tsv
#awk -F'\t' '{print NF; exit}' INF1_concat_matrices.tsv
#awk -F'\t' '{print NF; exit}' SLU1_concat_matrices.tsv
#awk -F'\t' '{print NF; exit}' EFF3_concat_matrices.tsv
#awk -F'\t' '{print NF; exit}' EFF2_concat_matrices.tsv
#awk -F'\t' '{print NF; exit}' INF2_concat_matrices.tsv
#awk -F'\t' '{print NF; exit}' INF3_concat_matrices.tsv
#awk -F'\t' '{print NF; exit}' SLU2_concat_matrices.tsv
awk -F'\t' '{print NF; exit}' SLU3_concat_matrices.tsv
# 494

# Check row names
#less EFF1_concat_matrices.tsv | cut -f 1 | head
#less INF1_concat_matrices.tsv | cut -f 1 | head
#less SLU1_concat_matrices.tsv | cut -f 1 | head
#less EFF3_concat_matrices.tsv | cut -f 1 | head
#less EFF2_concat_matrices.tsv | cut -f 1 | head
#less INF2_concat_matrices.tsv | cut -f 1 | head
#less INF3_concat_matrices.tsv | cut -f 1 | head
#less SLU2_concat_matrices.tsv | cut -f 1 | head
less SLU3_concat_matrices.tsv | cut -f 1 | head

# Add the final column headers
## *This is just given now based on previous analyses, could be a function to make it*

# The Common_id etc starts at 494-
cat ../../header.tsv | cut -f 1-493 > tmp
#cat tmp EFF1_concat_matrices.tsv > temp_EFF1_concat_matrices.tsv && mv temp_EFF1_concat_matrices.tsv EFF1_concat_matrices.tsv
#cat tmp INF1_concat_matrices.tsv > temp_INF1_concat_matrices.tsv && mv temp_INF1_concat_matrices.tsv INF1_concat_matrices.tsv
#cat tmp SLU1_concat_matrices.tsv > temp_SLU1_concat_matrices.tsv && mv temp_SLU1_concat_matrices.tsv SLU1_concat_matrices.tsv
#cat tmp EFF3_concat_matrices.tsv > temp_EFF3_concat_matrices.tsv && mv temp_EFF3_concat_matrices.tsv EFF3_concat_matrices.tsv
#cat tmp EFF2_concat_matrices.tsv > temp_EFF2_concat_matrices.tsv && mv temp_EFF2_concat_matrices.tsv EFF2_concat_matrices.tsv
#cat tmp INF2_concat_matrices.tsv > temp_INF2_concat_matrices.tsv && mv temp_INF2_concat_matrices.tsv INF2_concat_matrices.tsv
#cat tmp INF3_concat_matrices.tsv > temp_INF3_concat_matrices.tsv && mv temp_INF3_concat_matrices.tsv INF3_concat_matrices.tsv
#cat tmp SLU2_concat_matrices.tsv > temp_SLU2_concat_matrices.tsv && mv temp_SLU2_concat_matrices.tsv SLU2_concat_matrices.tsv
cat tmp SLU3_concat_matrices.tsv > temp_SLU3_concat_matrices.tsv && mv temp_SLU3_concat_matrices.tsv SLU3_concat_matrices.tsv

# add header 'sample' as the last column name

# Check again
#less EFF1_concat_matrices.tsv  | cut -f 1 | head
#less INF1_concat_matrices.tsv  | cut -f 1 | head
#less SLU1_concat_matrices.tsv  | cut -f 1 | head
#less EFF3_concat_matrices.tsv  | cut -f 1 | head
#less EFF2_concat_matrices.tsv  | cut -f 1 | head
#less INF2_concat_matrices.tsv  | cut -f 1 | head
#less INF3_concat_matrices.tsv  | cut -f 1 | head
#less SLU2_concat_matrices.tsv  | cut -f 1 | head
less SLU3_concat_matrices.tsv  | cut -f 1 | head

# Check number of cols
#awk -F'\t' '{print NF; exit}' EFF1_concat_matrices.tsv
#awk -F'\t' '{print NF; exit}' INF1_concat_matrices.tsv
#awk -F'\t' '{print NF; exit}' SLU1_concat_matrices.tsv
#awk -F'\t' '{print NF; exit}' EFF3_concat_matrices.tsv
#awk -F'\t' '{print NF; exit}' EFF2_concat_matrices.tsv
#awk -F'\t' '{print NF; exit}' INF2_concat_matrices.tsv
#awk -F'\t' '{print NF; exit}' INF3_concat_matrices.tsv
#awk -F'\t' '{print NF; exit}' SLU2_concat_matrices.tsv
awk -F'\t' '{print NF; exit}' SLU3_concat_matrices.tsv

# Rename according to included modifications
#mv EFF1_concat_matrices.tsv EFF1_concat_matrices_top100.tsv
#mv INF1_concat_matrices.tsv INF1_concat_matrices_top100.tsv
#mv SLU1_concat_matrices.tsv SLU1_concat_matrices_top100.tsv
#mv EFF3_concat_matrices.tsv EFF3_concat_matrices_top100.tsv
#mv EFF2_concat_matrices.tsv EFF2_concat_matrices_top100.tsv
#mv INF2_concat_matrices.tsv INF2_concat_matrices_top100.tsv
#mv INF3_concat_matrices.tsv INF3_concat_matrices_top100.tsv
#mv SLU2_concat_matrices.tsv SLU2_concat_matrices_top100.tsv
mv SLU3_concat_matrices.tsv SLU3_concat_matrices_top100.tsv
```
