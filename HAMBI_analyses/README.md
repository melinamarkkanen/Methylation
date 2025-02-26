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
- fill in column headers as follows:
```
contig_id       d       p	c	o	f	g	s	element	str	All	Domain to species	Domain to genus	Domain to family	Domain to order	Domain to class	Domain to phylum
---------------------------	-----------	---------------------	----------------------	-------------------	------------------	-------------	----------------------------	--------------	----------	--------------------------------------------------------------------------------------------------	----------------------------------------------------------------------------------------------------------------------------------------------------------------	---------------------------------------------------------------------------------	------------------------------------------------------------	-----------------------------------------------------	-----------------------------------	-------------------------
bcAd1023T--bcAd1023T_ptg000001l	Bacteria	Pseudomonadota	Gammaproteobacteria	Burkholderiales	Burkholderiaceae_B	Comamonas	Comamonas testosteroni_C	chromosome	HAMBI_0403	Bacteria_Pseudomonadota_Gammaproteobacteria_Burkholderiales_Burkholderiaceae_B_Comamonas_Comamonas_testosteroni_C_chromosome_HAMBI_0403	Bacteria_Pseudomonadota_Gammaproteobacteria_Burkholderiales_Burkholderiaceae_B_Comamonas_Comamonas_testosteroni_C	Bacteria_Pseudomonadota_Gammaproteobacteria_Burkholderiales_Burkholderiaceae_B_Comamonas	Bacteria_Pseudomonadota_Gammaproteobacteria_Burkholderiales_Burkholderiaceae_B	Bacteria_Pseudomonadota_Gammaproteobacteria_Burkholderiales	Bacteria_Pseudomonadota_Gammaproteobacteria	Bacteria_Pseudomonadota
bcAd1023T--bcAd1023T_ptg000003l	Bacteria	Pseudomonadota	Gammaproteobacteria	Enterobacterales	Enterobacteriaceae	Kluyvera	Kluyvera intermedia	plasmid unnamed	HAMBI_1299	Bacteria_Pseudomonadota_Gammaproteobacteria_Enterobacterales_Enterobacteriaceae_Kluyvera_Kluyvera_intermedia_plasmid_unnamed_HAMBI_1299	Bacteria_Pseudomonadota_Gammaproteobacteria_Enterobacterales_Enterobacteriaceae_Kluyvera_Kluyvera_intermedia	Bacteria_Pseudomonadota_Gammaproteobacteria_Enterobacterales_Enterobacteriaceae_Kluyvera	Bacteria_Pseudomonadota_Gammaproteobacteria_Enterobacterales_Enterobacteriaceae	Bacteria_Pseudomonadota_Gammaproteobacteria_Enterobacterales	Bacteria_Pseudomonadota_Gammaproteobacteria	Bacteria_Pseudomonadota
…																																					
```
- transfer back to Puhti for further modifications
```
# Add the contigs with no label and fill in 'NA'
cd HAMBI_data/metagenomic_assembly
cat *contigs.fasta | grep ">" | sed 's/>//' > contig_IDs.txt
cut -f 1 HAMBI_labels.txt > label_IDs.txt
grep -v -f label_IDs.txt contig_IDs.txt > missing_IDs.txt
less missing_IDs.txt | wc -l
# 669
```

#### Create file for missing contigs
#### Run ./create_file.sh
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
cat HAMBI_labels.txt missing_data.txt > tmp && mv tmp HAMBI_labels.txt
```
- bind to methylation data to create ```merged_data.tsv```
```
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
- x
- y

```
script here
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