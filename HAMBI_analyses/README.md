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
### Modify ```HAMBI_taxa_names.txt``` further in excel 
#### - modify the labels according to GTDB nomenclature
#### - fill in column headers:
```
contig	fullname	accession	element	domain	phylum	class	order	family	genus	species	secondary_accession	sample
---------------------------	------------------------	----------	-------------	--------------------	------	-------------	-----------------	------------	---------------	----------	--------------------	------------------
bcAd1023T--bcAd1023T_ptg000001l	HAMBI_0403_chromosome	HAMBI_0403	chromosome	HAMBI_0403_chromosome	Bacteria	Pseudomonadota	Gammaproteobacteria	Burkholderiales	Burkholderiaceae_B	Comamonas	Comamonas testosteroni_C	bcAd1023T--bcAd1023T
bcAd1023T--bcAd1023T_ptg000003l	HAMBI_1299_plasmid unnamed	HAMBI_1299	plasmid unnamed	HAMBI_1299_plasmid unnamed	Bacteria	Pseudomonadota	Gammaproteobacteria	Enterobacterales	Enterobacteriaceae	Kluyvera	Kluyvera intermedia	bcAd1023T--bcAd1023T
...		
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