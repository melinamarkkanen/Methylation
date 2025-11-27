# Analysis steps for the synthetic community data
&nbsp;
## Preparatory analysis
### Create database ```HAMBI_genomes.fasta``` for community members with WGS data
```
# Download available HAMBI assemblies from NCBI (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1047486/)
cd HAMBI_data/WGS_data
cat *fna > HAMBI_genomes.fasta
```
### Run ```workflow/Snakefile_HAMBI_preanalysis``` for:
- generation of HiFi reads **without** kinetics tags (fastq.gz) for the assembly
- metagenomic assembly of the community
- BLASTn search between the metagenomic assemblies and WGS data
- generation of HiFi reads **with** kinetics tags (.bam) for the methylation analysis
```
module load snakemake
snakemake --profile workflow/profile --use-envmodules --use-singularity \
        --snakefile workflow/Snakefile_HAMBI_preanalysis --use-singularity \
        --singularity-args "--bind ~/Raw_subread_data/m64145_231126_001443.subreads.bam" -np
```
&nbsp;
## Preliminary analysis for the methylation detection
### Extract contigs for methylation analysis
```
# Extract metagenome assembled contigs from each samples into individual files in ```contigs/```
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
### Run ```workflow/Snakefile_HAMBI_methylation_analysis```:
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
## Generate Position Weight Matrices (PWM):
- to filter data, the methylation types that have less than 50 detected sites are filled with 0 matrices which increased the models performance
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
python3 src/scoring_matrices_HAMBI.py HAMBI_data/contigs/bcAd1023T--bcAd1023T HAMBI_data/bcAd1023T_matrices

# Clean subfolders
cd /scratch/project_2006608/Methylation/HAMBI_data/contigs
rm -r bcAd1023T--bcAd1023T/
...
```
### Create ```HAMBI_labels.txt``` for connecting taxonomical identities and methylation profiles of metagenomic contigs for method vaidation
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
#### Modify ```HAMBI_taxa_names.txt``` further in excel, save as ```HAMBI_labels.txt```:
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
#### Process missing values in ```HAMBI_labels.txt```
```
# count contigs missing values
cd HAMBI_data/metagenomic_assembly
cat *contigs.fasta | grep ">" | sed 's/>//' > contig_IDs.txt
cut -f 1 HAMBI_labels.txt > label_IDs.txt
grep -v -f label_IDs.txt contig_IDs.txt > missing_IDs.txt
less missing_IDs.txt | wc -l
# 669
```
#### Run ```./create_file.sh``` to create data entries for contigs with no labels and fill in 'NA'
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
#### Append them to the original file
```
paste -d '\t' missing_IDs.txt missing_values.txt > missing_data.txt
cat missing_data.txt >> HAMBI_labels.txt
```
#### Combine PWMs with ```HAMBI_labels.txt```
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
## Random Forest classifier
### Run ```notebooks/Random_Forest_analysis_HAMBI.ipynb```
&nbsp;
## Sequence logos
```
export PYTHONUSERBASE=/projappl/project_2009999/my-python-env
module load python-data
python3 src/create_logos.py HAMBI_data/contigs HAMBI_data/logos
```
&nbsp;
## Taxonomical composition
```
# Run Sylph 
cd HAMBI_data/
apptainer_wrapper exec sylph profile gtdb_database.syldb HiFi_fastq/*.hifi_reads.fastq.gz -t $SLURM_CPUS_PER_TASK > Sylph_out/Sylph_HAMBI.tsv

# Run Sylph-tax
apptainer_wrapper exec sylph-tax taxprof Sylph_out/Sylph_HAMBI.tsv -t GTDB_r214 -o Sylph_out/prefix_

# Merge samples
apptainer_wrapper exec sylph-tax merge *.sylphmpa --column relative_abundance -o Sylph_HAMBI_merged.txt

# Get genus
awk '$1 ~ "clade_name" || $1 ~ "g__" {print $0}' Sylph_HAMBI_merged.txt | grep -v "|t__" | grep -v "|s__"  > Sylph_HAMBI_merged_genus_full.txt

# Edit
sed -i 's/HiFi_fastq\///g' Sylph_HAMBI_merged_genus_full.txt
sed -i 's/\.hifi_reads\.fastq\.gz//g' Sylph_HAMBI_merged_genus_full.txt
sed -i 's/--bcAd10[0-9][0-9]T//g' Sylph_HAMBI_merged_genus_full.txt
```