# Method validation
&nbsp;
## *Folder structure?*
&nbsp;
## Metagenome assemmbly of HAMBI community
### Creating database ```HAMBI_genomes.fasta``` for community members with WGS data
```
# Download available HAMBI assemblies from NCBI (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1047486/)
cd HAMBI_data/WGS_data
cat *fna > HAMBI_genomes.fasta
```
### Running ```workflow/Snakefile_HAMBI_preanalysis```: *(should we add snakemake dry run draws?)*
- generation of HiFi reads **without** kinetics tags (fastq.gz) for the assembly
- metagenomic assembly of the community
- BLASTn search between the metagenomic assemblies and WGS data
- generation of HiFi reads **with** kinetics tags (.bam) for the methylation analysis
```
module load snakemake
snakemake --profile workflow/profile --use-envmodules --use-singularity \
        --snakefile workflow/Snakefile_HAMBI_preanalysis --use-singularity \
        --singularity-args "--bind /scratch/project_2006608/Veera_PacBio/Raw_subread_data/m64145_231126_001443.subreads.bam" -np
```
&nbsp;
&nbsp;
## Creating ```HAMBI_labels.txt```
```
# Combine results
cd HAMBI_data/metagenomic_assembly
cat *txt > blast_out.txt

# Summarize by one hit/contig
awk '!seen[$1]++' blast_out.txt > besthit_blast_out.txt

# Filter by the alignment length of >= 1000 bp
awk '$4 >= 1000' besthit_blast_out.txt > filtered_blast_out.txt

# Bind the species names
cd HAMBI_data/WGS_data
less HAMBI_genomes.fasta | grep ">" | sed 's/>//' > WGS_ID.txt
sed -i 's/\.1 /\.1\t/g' WGS_ID.txt
awk 'NR==FNR {a[$1]=$0; next} $2 in a {print $0, a[$2]}' WGS_ID.txt ../metagenomic_assembly/filtered_blast_out.txt > ../metagenomic_assembly/HAMBI_taxa_names.txt
```
### Modify ```HAMBI_taxa_names.txt``` further in excel, save as ```HAMBI_labels.txt```:
- modify the labels according to GTDB nomenclature
- fill in column headers (see below)
- transfer back to Puhti for further modifications
```
contig_id       d       p	c	o	f	g	s	element	str	All	Domain to species	Domain to genus	Domain to family	Domain to order	Domain to class	Domain to phylum
---------------------------     -----------     --------------------    ------------------      -------------------	------------------	-------------	----------------------------	--------------	----------	--------------------------------	----------------------------	------------------	---------------------------	-------------------------	---------	-------------------------
bcAd1023T--bcAd1023T_ptg000001l	Bacteria	Pseudomonadota	Gammaproteobacteria	Burkholderiales	Burkholderiaceae_B	Comamonas	Comamonas testosteroni_C	chromosome	HAMBI_0403	Bacteria_Pseudomonadota_Gammaproteobacteria_Burkholderiales_Burkholderiaceae_B_Comamonas_Comamonas_testosteroni_C_chromosome_HAMBI_0403	Bacteria_Pseudomonadota_Gammaproteobacteria_Burkholderiales_Burkholderiaceae_B_Comamonas_Comamonas_testosteroni_C	Bacteria_Pseudomonadota_Gammaproteobacteria_Burkholderiales_Burkholderiaceae_B_Comamonas	Bacteria_Pseudomonadota_Gammaproteobacteria_Burkholderiales_Burkholderiaceae_B	Bacteria_Pseudomonadota_Gammaproteobacteria_Burkholderiales	Bacteria_Pseudomonadota_Gammaproteobacteria	Bacteria_Pseudomonadota
bcAd1023T--bcAd1023T_ptg000003l	Bacteria	Pseudomonadota	Gammaproteobacteria	Enterobacterales	Enterobacteriaceae	Kluyvera	Kluyvera intermedia	plasmid unnamed	HAMBI_1299	Bacteria_Pseudomonadota_Gammaproteobacteria_Enterobacterales_Enterobacteriaceae_Kluyvera_Kluyvera_intermedia_plasmid_unnamed_HAMBI_1299	Bacteria_Pseudomonadota_Gammaproteobacteria_Enterobacterales_Enterobacteriaceae_Kluyvera_Kluyvera_intermedia	Bacteria_Pseudomonadota_Gammaproteobacteria_Enterobacterales_Enterobacteriaceae_Kluyvera	Bacteria_Pseudomonadota_Gammaproteobacteria_Enterobacterales_Enterobacteriaceae	Bacteria_Pseudomonadota_Gammaproteobacteria_Enterobacterales	Bacteria_Pseudomonadota_Gammaproteobacteria	Bacteria_Pseudomonadota
…																																					
```
- count contigs missing values
```
cd HAMBI_data/metagenomic_assembly
cat *contigs.fasta | grep ">" | sed 's/>//' > contig_IDs.txt
cut -f 1 HAMBI_labels.txt > label_IDs.txt
grep -v -f label_IDs.txt contig_IDs.txt > missing_IDs.txt
less missing_IDs.txt | wc -l
# 669
```

- create the contigs with no label and fill in 'NA'
- run ```./create_file.sh```
```
#!/bin/bash

# Function to create the file
create_file() {
    local n_rows=$1      # Number of rows
    local n_columns=$2   # Number of columns
    local string=$3      # The string to repeat
    local filename=$4    # Output file name

    # Create the file with the specified name
    for ((i = 0; i < n_rows; i++))
    do
        # Create a line with the string repeated n_columns times, separated by tabs
        line=$(for ((j = 0; j < n_columns; j++)); do echo -n "$string"; if (( j < n_columns-1 )); then echo -n -e "\t"; fi; done)
        # Write the line to the file
        echo -e "$line" >> "$filename"
    done
}

# Create a file with 699 rows and 16 columns, repeating the string "NA"
create_file 669 16 "NA" "missing_values.txt"
```
- append them to the original file
```
paste -d '\t' missing_IDs.txt missing_values.txt > missing_data.txt
cat missing_data.txt >> HAMBI_labels.txt
```




&nbsp;
&nbsp;
## Methylation analysis of HAMBI community
### Extract contigs for methylation analysis
```
## Extract metagenome assembled contigs from each samples into individual files in ```contigs/```
cd HAMBI_data/metagenomic_assembly
cp *_contigs.fasta ../contigs
cd ../contigs/
module load biokit
seqretsplit bcAd1023T--bcAd1023T_contigs.fasta
seqretsplit bcAd1037T--bcAd1037T_contigs.fasta
seqretsplit bcAd1039T--bcAd1039T_contigs.fasta
seqretsplit bcAd1046T--bcAd1046T_contigs.fasta
seqretsplit bcAd1063T--bcAd1063T_contigs.fasta
rm *_contigs.fasta

# Rename the extracted contigs to match the original names
rename bcad1023t--bcad1023t bcAd1023T--bcAd1023T *.fasta
rename bcad1037t--bcad1037t bcAd1037T--bcAd1037T *.fasta
rename bcad1039t--bcad1039t bcAd1039T--bcAd1039T *.fasta
rename bcad1046t--bcad1046t bcAd1046T--bcAd1046T *.fasta
rename bcad1063t--bcad1063t bcAd1063T--bcAd1063T *.fasta
```
### Running ```workflow/Snakefile_HAMBI_methylation_analysis```: *(should we add snakemake dry run draws?)*
- align the HiFi reads with kinetics to the assemblies
- run ipSummary to obtain the .gff files for downatream analyses
```
# Run first only the first rule all
module load snakemake
snakemake --profile workflow/profile --use-envmodules --use-singularity \
        --snakefile workflow/Snakefile_HAMBI_methylation_analysis --use-singularity -np

# With two latter **rule all**s add flag ```--keep-going``` to avoid crashing after every contig with no kinetic reads mapped (=failed ipdSummary)
 module load snakemake
snakemake --profile workflow/profile --use-envmodules --use-singularity \
        --snakefile workflow/Snakefile_HAMBI_methylation_analysis --use-singularity --keep-going -np
```




&nbsp;
&nbsp;
## Create scoring matrices and flattened feature matrices
### Position Weight Matrices (PWM) 
- to filter data, the methylation types that have less than 20 detected sites are filled with 0 matrices which increased the models performance
- the scoring matrices are then flattened to feature matrices. The flattened feature matrices are then used to train the random forest model to predict the taxonomic classification of the contigs. 
```
# Prepare folders with .fasta & .gff files for each sample
cd /scratch/project_2006608/Methylation/HAMBI_data/contigs
mkdir bcAd1023T--bcAd1023T
cp bcAd1023T--bcAd1023T_*fasta bcAd1023T--bcAd1023T
cp bcAd1023T--bcAd1023T_*gff bcAd1023T--bcAd1023T

# Generate the matrices (interactive session)
module load python-data
cd /scratch/project_2006608/Methylation
python3 src/scoring_matrices.py HAMBI_data/contigs/bcAd1023T--bcAd1023T HAMBI_data/bcAd1023T_matrices

# Clean subfolders
cd /scratch/project_2006608/Methylation/HAMBI_data/contigs
rm -r bcAd1023T--bcAd1023T/
...
```
### Combine PWM with ```HAMBI_labels.txt```
```
# Combine modification types
cd HAMBI_data/bcAd1023T_matrices/flattened  # repeat for all samples

# Add shared row names
less m4C.tsv | cut -f 1 > common_id
# Remove extra contig names
cut -f 2-  m6A.tsv > tmp && mv tmp m6A.tsv
cut -f 2-  modified_base.tsv > tmp && mv tmp modified_base.tsv
# Paste
paste common_id m4C.tsv modified_base.tsv m6A.tsv > concat_matrices.tsv

# Combine all
cd HAMBI_data
cat bcAd1023T_matrices/flattened/concat_matrices.tsv bcAd1037T_matrices/flattened/concat_matrices.tsv bcAd1039T_matrices/flattened/concat_matrices.tsv bcAd1046T_matrices/flattened/concat_matrices.tsv bcAd1063T_matrices/flattened/concat_matrices.tsv > merged_data.tsv

# Check number of columns
awk -F'\t' '{print NF; exit}' merged_data.tsv
# 494

# Reorder PWM data according to taxa names
head -n 1 metagenomic_assembly/HAMBI_labels.txt > reord_merged_data.tsv
awk -F'\t' '{print NF; exit}' reord_merged_data.tsv
# 17
awk 'FNR==NR {x2[$1] = $0; next} $1 in x2 {print x2[$1]}' metagenomic_assembly/HAMBI_labels.txt merged_data.tsv >> reord_merged_data.tsv

# Check row names
less reord_merged_data.tsv | cut -f 1 | head
less merged_data.tsv  | cut -f 1 | head

# Remove the contig_id col from merged_data.tsv
cut -f 2- merged_data.tsv > tmp && mv tmp merged_data.tsv

# Add the final column headers
## *This is just given now based on previous analyses, could be a function to make it*

# The Common_id etc starts at 494-
cat header.tsv | cut -f 1-493 > tmp
cat tmp merged_data.tsv > temp_merged_data.tsv && mv temp_merged_data.tsv merged_data.tsv

# Check again
less reord_merged_data.tsv | cut -f 1 | head
less merged_data.tsv  | cut -f 1 | head

# Combine
paste merged_data.tsv reord_merged_data.tsv > tmp && mv tmp merged_data.tsv

# Check number of cols
awk -F'\t' '{print NF; exit}' merged_data.tsv

# There are still empty columns in the end, remove them
cut -f 1-510 merged_data.tsv > tmp && mv tmp merged_data.tsv
```
&nbsp;
&nbsp;
## Run Random Forest analyses ```notebooks/Random_forest_HAMBI```
### Sequence -20...20: 150 most important features
### Sequence -5...5: All features
```
cd /scratch/project_2006608/Methylation/HAMBI_data
# List of patterns to match against column names
## 5...-5
cut -f 1,17-27,58-68,99-109,140-150,181-191,222-232,263-273,304-314,345-355,386-396,427-437,468-478,494-510 merged_data.tsv > merged_data_short_5.tsv
## 10...-10
cut -f 1,12-34,53-73,94-114,135-155,176-196,217-237,258-278,299-319,340-360,381-401,422-442,463-483,494-510 merged_data.tsv > merged_data_short_10.tsv
```
&nbsp;
&nbsp;
&nbsp;

## Create Sequence logos


## UMAP
### Build MAGs according to contigs clustered by UMAP
#### Get contigs
```
cd ../src
./HAMBI_get_contigs_for_MAGs.sh C1A
```

#### Combine small contigs sample-wise
```
cd /scratch/project_2006608/Methylation/HAMBI_data/MAGs
# Check samples
ls */*contigs/
# Check contig lengths
tail -n 50 C4*/C4*txt
# Combine
## C4
cat C4/bcAd1023T--bcAd1023T_contigs/*fasta > C4/bcAd1023T--bcAd1023T_contigs/bcAd1023T--bcAd1023T_C4.fasta
cat C4/bcAd1037T--bcAd1037T_contigs/*fasta > C4/bcAd1037T--bcAd1037T_contigs/bcAd1037T--bcAd1037T_C4.fasta
cat C4/bcAd1039T--bcAd1039T_contigs/*fasta > C4/bcAd1039T--bcAd1039T_contigs/bcAd1039T--bcAd1039T_C4.fasta
cat C4/bcAd1046T--bcAd1046T_contigs/*fasta > C4/bcAd1046T--bcAd1046T_contigs/bcAd1046T--bcAd1046T_C4.fasta
cat C4/bcAd1063T--bcAd1063T_contigs/*fasta > C4/bcAd1063T--bcAd1063T_contigs/bcAd1063T--bcAd1063T_C4.fasta

# Remove these:
# bcAd1023T--bcAd1023T_ptg000007c	4577733
# bcAd1037T--bcAd1037T_ptg000008c	4577734
# bcAd1039T--bcAd1039T_ptg000024c	4577740
# bcAd1046T--bcAd1046T_ptg000003c	4579072
# bcAd1063T--bcAd1063T_ptg000008c	4577747

seqkit grep -v -n -f list.txt bcAd1023T--bcAd1023T_contigs/bcAd1023T--bcAd1023T_C4.fasta > tmp && mv tmp bcAd1023T--bcAd1023T_contigs/bcAd1023T--bcAd1023T_C4.fasta

seqkit grep -v -n -f list.txt bcAd1037T--bcAd1037T_contigs/bcAd1037T--bcAd1037T_C4.fasta > tmp && mv tmp bcAd1037T--bcAd1037T_contigs/bcAd1037T--bcAd1037T_C4.fasta

seqkit grep -v -n -f list.txt bcAd1039T--bcAd1039T_contigs/bcAd1039T--bcAd1039T_C4.fasta > tmp && mv tmp bcAd1039T--bcAd1039T_contigs/bcAd1039T--bcAd1039T_C4.fasta

seqkit grep -v -n -f list.txt bcAd1046T--bcAd1046T_contigs/bcAd1046T--bcAd1046T_C4.fasta > tmp && mv tmp bcAd1046T--bcAd1046T_contigs/bcAd1046T--bcAd1046T_C4.fasta

seqkit grep -v -n -f list.txt bcAd1063T--bcAd1063T_contigs/bcAd1063T--bcAd1063T_C4.fasta > tmp && mv tmp bcAd1063T--bcAd1063T_contigs/bcAd1063T--bcAd1063T_C4.fasta
```



#### Run CheckM2
````
# Create sample list
cd /scratch/project_2006608/Methylation/HAMBI_data/MAGs/C4
ls -d *contigs/ | sed 's/_contigs\///g' > sample_names.txt


# Set the variable
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_names.txt)

# Set temp dir
export TMPDIR="/scratch/project_2006608/Methylation/tmp_dir"

# Run
apptainer exec --bind $PWD:$PWD,$TMPDIR:/scratch/project_2006608/Methylation/tmp_dir,$CHECKM2DB:/scratch/project_2006608/Methylation_Viikki_HiFi/db/CheckM2_database/uniref100.KO.1.dmnd /projappl/project_2006608/containers/checkm2:1.0.1.sif checkm2 predict --input $sample"_contigs"/*.fasta \
        --output-directory $sample"_CheckM2_out" --extension fasta --threads 6 --force \
        --database_path /scratch/project_2006608/Methylation_Viikki_HiFi/db/CheckM2_database/uniref100.KO.1.dmnd

# Summarize
cat HAMBI_data/MAGs/C7/*_CheckM2_out/quality_report.tsv | cut -f 1-3 | grep -v "Name"
```

#### Check if plasmid or ARGs
```
cat HAMBI_labels.txt | grep -f C7_id.txt | grep -v "chromoso" | cut -f 1,9
cat *res* | grep -f C7_id.txt | cut -f 1-4
```



# Run GTDB-Tk
```
# Go to folder
cd /scratch/project_2006608/Methylation/HAMBI_data/MAGs/C1A

# Set variable
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_names.txt)

# Load the environment and variables
export PATH="/projappl/project_2006608/GTDB-Tk/bin:$PATH"
export GTDBTK_DATA_PATH=/scratch/project_2006608/GTDB-Tk/release220/

# Run
gtdbtk classify_wf --genome_dir $sample"_contigs" -x fasta \
        --out_dir $sample"_GTDB_out" --skip_ani_screen --cpus $SLURM_CPUS_PER_TASK

cat *GTDB_out/g*tsv | cut -f 1,5
```

### Check if the contig length matters in clustering in UMAP
```
# Load
module load seqkit/2.5.1

cd /scratch/project_2006608/Methylation/HAMBI_data/metagenomic_assembly

# Set environmental variable
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_names.txt)

# Run
seqkit fx2tab --length --name --header-line $sample".fasta" > $sample"_lengths.txt"

# Add sample name in the front 
sed -i -e 's/^/bcAd1023T--bcAd1023T_/' bcAd1023T--bcAd1023T_lengths.txt
sed -i -e 's/^/bcAd1037T--bcAd1037T_/' bcAd1037T--bcAd1037T_lengths.txt
sed -i -e 's/^/bcAd1039T--bcAd1039T_/' bcAd1039T--bcAd1039T_lengths.txt
sed -i -e 's/^/bcAd1046T--bcAd1046T_/' bcAd1046T--bcAd1046T_lengths.txt
sed -i -e 's/^/bcAd1063T--bcAd1063T_/' bcAd1063T--bcAd1063T_lengths.txt

# Extract those that are in notebooks/UMAP_HAMBI_top150_features/top150_contigs.txt
cd /scratch/project_2006608/Methylation/notebooks/UMAP_HAMBI_top150_features
cat ../../HAMBI_data/metagenomic_assembly/*_lengths.txt | grep -f top150_contigs.txt > top150_contigs_lengths.txt

# Reorder to match
awk 'NR==FNR {order[$1]=NR; next} {print order[$1], $0}' top150_contigs.txt top150_contigs_lengths.txt | sort -n | cut -d' ' -f2- > tmp
tr ' ' '\t' < tmp > top150_contigs_lengths.tsv
sed -i '1s/^/contig\tcontig_length\n/' top150_contigs_lengths.tsv
```

### Check if the the count of modifications matters in clustering in UMAP
```
cd src/
./HAMBI_count_gff_lines.sh

# Extract those that are in notebooks/UMAP_HAMBI_top150_features/top150_contigs.txt
cd /scratch/project_2006608/Methylation/notebooks/UMAP_HAMBI_top150_features
cat ../../HAMBI_data/contigs/HAMBI_mod_counts.txt | grep -f top150_contigs.txt > top150_contigs_mod_counts.txt

# Reorder to match
awk 'NR==FNR {order[$1]=NR; next} {print order[$1], $0}' top150_contigs.txt top150_contigs_mod_counts.txt | sort -n | cut -d' ' -f2- > tmp
tr ' ' '\t' < tmp > top150_contigs_mod_counts.tsv
sed -i '1s/^/contig\tmod_count\n/' top150_contigs_mod_counts.tsv
```




## Creating sequence logos (sinteractive session)
```
export PYTHONUSERBASE=/projappl/project_2009999/my-python-env
module load python-data
python3 src/create_logos.py HAMBI_data/bcAd1023T_matrices HAMBI_data/bcAd1023T_matrices/logos
```

## ARG Annotations
```
# Set variable
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p ../sample_names.txt)

# Run
blastn -query $sample"_contigs.fasta" \
        -subject ../../db/resfinder_db/all.fsa \
        -out $sample"_resfinder_out.txt" -outfmt 6
```

## REBASE
```
# Set variable
contig=$(sed -n ${SLURM_ARRAY_TASK_ID}p contig_names.txt)
contigShort=$(sed -n ${SLURM_ARRAY_TASK_ID}p contig_names.txt | sed 's/^[^/]*\///g' | sed 's/\.fasta//g' )

# Run
blastn -query $contig \
        -subject ../../../db/All_REBASE_Gold_Standards_DNA.fasta \
        -out REBASE_out/$contigShort"_REBASE_out.txt" -outfmt 6
```

## Analyse further the different methylation profiles of genomes of same species between different samples
### Bakta
```
export SING_IMAGE=/projappl/project_2006608/containers/bakta:1.11.0.sif
export BAKTA_DB=/scratch/project_2006608/Methylation/db/bakta_db/db-light

apptainer_wrapper exec bakta $name".fasta" \
	--prefix $name"_bakta" \
	--output $name"_bakta_out"/ \
	--db $BAKTA_DB \
	--keep-contig-headers \
	--threads $SLURM_CPUS_PER_TASK
```

### Align the reads to the reference genome
```
# Load tools
module load samtools/1.21

# Set temp dir
export TMPDIR="/scratch/project_2006608/Methylation/tmp_dir"

# Go to dir
cd /scratch/project_2006608/Methylation/HAMBI_data/MAGs/Paracoccus

# Set variable
name=$(sed -n ${SLURM_ARRAY_TASK_ID}p names.txt)

# Run alignment
apptainer exec --bind $PWD:$PWD,$TMPDIR:/scratch/project_2006608/Methylation/tmp_dir,$DATA:/scratch/project_2006608/Methylation/HAMBI_data/HiFi_bam_kinetic/ \
        /projappl/project_2006608/containers/pbmm2:1.17.0.sif pbmm2 align GCF_034627565.1_ASM3462756v1_genomic.fasta /scratch/project_2006608/Methylation/HAMBI_data/HiFi_bam_kinetic/"hifi_reads_kinetic."$name".bam" \
        $name"_Paracoccus_mapped_reads".bam --sort --preset CCS --log-level INFO --log-file $name"_mapping.log" \
        --num-threads $SLURM_CPUS_PER_TASK

# Index
samtools index $name"_Paracoccus_mapped_reads.bam"
```

# Check the taxa of those contigs that are not among our WGS data
### Remote BLAST
```
# Create folder
cd /scratch/project_2006608/Methylation/HAMBI_data/metagenomic_assembly
mkdir Missing_labels

# Split names
split -l 290 missing_IDs.txt --numeric-suffixes missing_IDs_

# Go to directory
cd /scratch/project_2006608/Methylation/HAMBI_data/metagenomic_assembly

# Set variable
name=$(sed -n ${SLURM_ARRAY_TASK_ID}p missing_IDs_00)
name=$(sed -n ${SLURM_ARRAY_TASK_ID}p missing_IDs_01)
name=$(sed -n ${SLURM_ARRAY_TASK_ID}p missing_IDs_02)

# Load tools
module load biokit

# Run
blastn -query $name".fasta" -db nt -remote -out Missing_labels/$name"_blast_out.txt" -outfmt 6 -perc_identity 90 -max_target_seqs 1
```

## Investigate the population structure more in detail
### Metaphlan4
```
# Load
module load metaphlan/4.1.1

# Set environmental variable
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_names.txt)

# Run
metaphlan /scratch/project_2006608/Methylation/HAMBI_data/HiFi_fastq/$sample".hifi_reads.fastq.gz" \
        --nproc $SLURM_CPUS_PER_TASK \
        --bowtie2out Metaphlan4_out/$sample".bowtie2.bz2" \
        --bowtie2db ../db/metaphlan_databases/ \
        --input_type fastq > Metaphlan4_out/$sample"_profile.txt"
```

### Kaiju
```
cd /scratch/project_2006608/Methylation/db
mkdir kaijudb
cd kaijudb
kaiju-makedb -s refseq

# Go to dir
cd /scratch/project_2006608/Methylation/HAMBI_data

# Load
module load kaiju/1.10.0

# Set temp dir
export TMPDIR="/scratch/project_2006608/Methylation/tmp_dir"

# Set environmental variable
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_names.txt)

# Run
kaiju -z $SLURM_CPUS_PER_TASK -t $KAIJUDB/nodes.dmp -f $KAIJUDB/refseq/kaiju_db_refseq.fmi \
	-i /scratch/project_2006608/Methylation/HAMBI_data/HiFi_fastq/$sample".hifi_reads.fastq.gz" \
	-o Kaiju_out/$sample"_kaiju_out" -v
```

## Miksi niin moni plasmidi hukataan ??