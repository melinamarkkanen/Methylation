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
- to filter data, the methylation types that have less than **100**!! detected sites are filled with 0 matrices which increased the models performance
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

# make copies of the .tsv files....
cp m6A.tsv cp_m6A.tsv
cp m4C.tsv cp_m4C.tsv
cp modified_base.tsv cp_modified_base.tsv

# Add shared row names
less m4C.tsv | wc -l
#yes "EFF1" | head -n 60373 > sample
#yes "INF1" | head -n 64737 > sample
yes "SLU1" | head -n 74387 > sample

# Remove extra contig names
cut -f 2-  m6A.tsv > tmp && mv tmp m6A.tsv
cut -f 2-  modified_base.tsv > tmp && mv tmp modified_base.tsv
# Paste
#paste m4C.tsv modified_base.tsv m6A.tsv sample > EFF1_concat_matrices.tsv
#paste m4C.tsv modified_base.tsv m6A.tsv sample > INF1_concat_matrices.tsv
paste m4C.tsv modified_base.tsv m6A.tsv sample > SLU1_concat_matrices.tsv

# Check number of columns
#awk -F'\t' '{print NF; exit}' EFF1_concat_matrices.tsv
#awk -F'\t' '{print NF; exit}' INF1_concat_matrices.tsv
awk -F'\t' '{print NF; exit}' SLU1_concat_matrices.tsv
# 494

# Check row names
#less EFF1_concat_matrices.tsv | cut -f 1 | head
#less INF1_concat_matrices.tsv | cut -f 1 | head
less SLU1_concat_matrices.tsv | cut -f 1 | head

# Add the final column headers
## *This is just given now based on previous analyses, could be a function to make it*

# The Common_id etc starts at 494-
cat ../../header.tsv | cut -f 1-493 > tmp
#cat tmp EFF1_concat_matrices.tsv > temp_EFF1_concat_matrices.tsv && mv temp_EFF1_concat_matrices.tsv EFF1_concat_matrices.tsv
#cat tmp INF1_concat_matrices.tsv > temp_INF1_concat_matrices.tsv && mv temp_INF1_concat_matrices.tsv INF1_concat_matrices.tsv
cat tmp SLU1_concat_matrices.tsv > temp_SLU1_concat_matrices.tsv && mv temp_SLU1_concat_matrices.tsv SLU1_concat_matrices.tsv

# add header 'sample' as the last column name

# Check again
#less EFF1_concat_matrices.tsv  | cut -f 1 | head
#less INF1_concat_matrices.tsv  | cut -f 1 | head
less SLU1_concat_matrices.tsv  | cut -f 1 | head

# Check number of cols
#awk -F'\t' '{print NF; exit}' EFF1_concat_matrices.tsv
#awk -F'\t' '{print NF; exit}' INF1_concat_matrices.tsv
awk -F'\t' '{print NF; exit}' SLU1_concat_matrices.tsv
# Rename according to included modifications
#mv EFF1_concat_matrices.tsv EFF1_concat_matrices_top100.tsv
#mv INF1_concat_matrices.tsv INF1_concat_matrices_top100.tsv
mv SLU1_concat_matrices.tsv SLU1_concat_matrices_top100.tsv
```
## Let's create mod count data also for the WW data set
```
cd src/
./WW_count_gff_lines.sh EFF1
sed -i '1i contig\tmod_count' EFF1_mod_counts.txt
```
## Create ARG annottaions to be visualized in UMAP
### WW_blastn_resfinder.sh
```
# Load the tools
module load biokit

# Go to dir
cd /scratch/project_2006608/Methylation/HAMBI_data/metagenomic_assembly

# Set variable
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p ../sample_names.txt)

# Run
blastn -query EFF1_contigs.fasta \
        -subject ../../db/resfinder_db/all.fsa \
        -out EFF1_resfinder_out.txt -outfmt 6 \
		-perc_identity 90 -max_target_seqs 5
```
## Filter & count ARGs / contig
```
awk '$4 >= 100' EFF1_resfinder_out.txt" > tmp && mv tmp EFF1_resfinder_out.txt

cp EFF1_resfinder_out.txt raw_EFF1_resfinder_out.txt

##### First manually (set # and then remove those lines) check hits that are duplicates

sed -i '/^#/d' EFF1_resfinder_out.txt
sed -i '/^#/d' EFF2_resfinder_out.txt
sed -i '/^#/d' EFF3_resfinder_out.txt

sed -i '/^#/d' INF2_resfinder_out.txt

# Then run
cd src/
./WW_count_ARGs.sh EFF1
```

## Attach ARG names
```
# Extract columns
cut -f 1-2 EFF1_resfinder_out.txt > ARG_names.txt
# Remove the accession
sed -i 's/\(.*_.*\)_.*$/\1/' ARG_names.txt
# Check those that have multiple
cut -f 1 ARG_names.txt | sort | uniq -d

##### Manually fix those contigs that have multiple

# Add contig names of those where there is no ARG
cut -f 1 ARG_names.txt > contig_ARG_names.txt
grep -v -f contig_ARG_names.txt contig_names.txt > remaining_contig_names.txt
cat ARG_names.txt remaining_contig_names.txt > EFF1_ARG_names.txt
sed -i '1i contig\tARG_name' EFF1_ARG_names.txt
```

## Attach contig lengths
```
module load seqkit/2.5.1
seqkit fx2tab --length --name --header-line EFF1_contigs.fasta > EFF1_contigs_lengths.txt
```


## Analyse by sample in Jupiter
## Extract clusters
### eg.: ./WW_get_contigs_for_MAGs.sh <sample> <cluster> <above>
```
./WW_get_contigs_for_MAGs.sh EFF1 C1
```

## Check ARGs in HiFi MAG Pipeline results
```
nano search.txt # (contig names excel)
grep -f search.txt EFF1_resfinder_out.txt
```






