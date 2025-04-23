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



## Clinker for erm(F)_3 which are present in the UMAP
# Get fasta (make sure that contigs are found in this folder)
```
nano EFF2_lista.txt

for i in $(less EFF2_lista.txt);do grep -A 1 -f <(echo "$i") EFF2_contigs.fasta > "EFF2_"$i".fasta";done
for i in $(less INF1_lista.txt);do grep -A 1 -f <(echo "$i") INF1_contigs.fasta > "INF1_"$i".fasta";done
for i in $(less INF3_lista.txt);do grep -A 1 -f <(echo "$i") INF3_contigs.fasta > "INF3_"$i".fasta";done
for i in $(less SLU1_lista.txt);do grep -A 1 -f <(echo "$i") SLU1_contigs.fasta > "SLU1_"$i".fasta";done
for i in $(less SLU2_lista.txt);do grep -A 1 -f <(echo "$i") SLU2_contigs.fasta > "SLU2_"$i".fasta";done
for i in $(less SLU3_lista.txt);do grep -A 1 -f <(echo "$i") SLU3_contigs.fasta > "SLU3_"$i".fasta";done


# Run Bakta
WW_bakta.sh

# Get locations of erm(F)_3 based on bakta and then ResFinder 
cd /scratch/project_2006608/Methylation/WW_data/erm_F_UMAP
cat *_bakta_out/*l.tsv | grep "erm(F)" > erm_F_location.txt

# Tidy
awk '{print $1,$3,$4}' erm_F_location.txt > tmp && mv tmp erm_F_location.txt
sed -i 's/ /\t/g' erm_F_location.txt

# Check if we have all erm(F) location
cut -f1 erm_F_location.txt | sort > ID_erm_F_location.txt

sed 's/^[A-Z][A-Z][A-Z][0-9]_//g' contig_names.txt | sort > sorted_contig_names.txt
comm -3 sorted_contig_names.txt ID_erm_F_location.txt > non_matching.txt

        s1.ctg015106l (duplicated)
s2862.ctg004065l (-> BLAST) # 78468   79266
s31155.ctg040907l (-> BLAST) # 20099   20898
        s37987.ctg058715l (duplicated)
        s3.ctg008624l (duplicated)
        s4.ctg007717l (duplicated)
        s4.ctg054975l (duplicated)
        s4.ctg055934l (duplicated)
s4.ctg060290l (-> BLAST) # 24349   25073
s854.ctg000934l (-> BLAST) # 100396  101192
s8619.ctg012505l (-> BLAST) # 24847   25648

# Blastn in interactive
blastn -query SLU1_s8619.ctg012505l.fasta \
        -subject ../../db/resfinder_db/all.fsa \
        -out SLU1_s8619.ctg012505l_resfinder_out.txt -outfmt 6

# Fix the duplicates so that there is the start of the first aand end of the second

cd src/
./WW_get_regions.sh ../WW_data/erm_F_UMAP/erm_F_location.txt

# replace negative with one
cd /scratch/project_2006608/Methylation/WW_data/erm_F_UMAP
less erm_F_location_flanking.txt | grep "-"

sed -i 's/-[0-9][0-9][0-9]/1/g' erm_F_location_flanking.txt
sed -i 's/-[0-9][0-9]/1/g' erm_F_location_flanking.txt
less erm_F_location_flanking.txt | grep "-"

rm startFlank.txt
rm endFlank.txt

# Create new accession list
cut -f 1 erm_F_location_flanking.txt  > updated_contig_names.txt
less updated_contig_names.txt | sort | uniq > tmp && mv tmp updated_contig_names.txt

mkdir extracted_data

cd src
./WW_extract_regions.sh

# Run Bakta for extracted_data/
src/WW_bakta.sh


# Proovframe
cd extracted_data
mkdir proovframe

# gather
cat *bakta_out/*.fna > proovframe/all_erm_F.fasta
cat *bakta_out/*l.faa > proovframe/reference_proteins.faa

module load cdhit/4.8.1
cd-hit -i proovframe/reference_proteins.faa -o proovframe/reference_proteins_proovframe.faa -c 0.90

module load diamond/2.0.15
/projappl/project_2006608/proovframe/bin/proovframe map -a proovframe/reference_proteins_proovframe.faa -o proovframe/raw-seqs.tsv proovframe/all_erm_F.fasta

/projappl/project_2006608/proovframe/bin/proovframe fix proovframe/all_erm_F.fasta proovframe/raw-seqs.tsv -o proovframe/erm_F_proovframe.fasta

# chgeck lenghths
seqkit fx2tab --length --header-line --name proovframe/erm_F_proovframe.fasta > proovframe/erm_F_proovframe_lengths.txt

cd proovframe
seqretsplit erm_F_proovframe.fasta


# Run Bakta for extracted_data/
src/WW_bakta.sh

export SING_IMAGE=/projappl/project_2006608/containers/clinker-py:0.0.27.sif
apptainer_wrapper exec clinker filt_clinker_in/*.gff3 \
	-p filt_erm_F_UMAP_clinker.html \
	-o filt_erm_F_UMAP_clinker \
	-j $SLURM_CPUS_PER_TASK

# Manually delete irrelevatnt or similar?
s4.ctg054887l
...
```

### Sylph
```
# head -n 1 Sylph_HiFi_merged.txt > Sylph_HiFi_merged_genus.txt (?)
awk '$1 ~ "clade_name" || $1 ~ "g__" {print $0}' Sylph_HiFi_merged.txt | grep -v "|t__" | grep -v "|s__"  > Sylph_HiFi_merged_genus_full.txt
```

### geNomad for ermF contigs
```
apptainer exec --bind $PWD:$PWD,$DB_PATH:$DB_PATH \
        /projappl/project_2006608/containers/genomad:1.8.0.sif genomad end-to-end \
	--cleanup --splits 8 $contig geNomad_out/$contig".fasta" $DB_PATH
```


### fARGene analysis
```
# Go to dir
cd /scratch/project_2006608/Methylation/WW_data

# Export program
export PATH="/projappl/project_2006608/containers/fargene/bin:$PATH"

# Set temp dir
export TMPDIR="/scratch/project_2006608/Methylation/tmp_dir"

# Run
fargene -i fARGene_in/*.fasta \
        --hmm-model /scratch/project_2006608/Methylation/db/fargene/fargene_analysis/models/aminoglycoside_model_i.hmm \
        --score 100 \
        -o fARGene_aminoglycoside_model_i_out -p $SLURM_CPUS_PER_TASK \
```

## Prepare results for UMAP
