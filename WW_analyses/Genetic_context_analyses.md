# ARG genetic contexts of selected ARG bins
- [Class D 2 beta-lactamase](#class-d-2-beta-lactamase)
- [Class C beta-lactamase](#class-c-beta-lactamase)
- [*bla*OXA-129](#blaoxa-129)
- [*sul1_9*-like](#*sul1_9*-like)
- [*erm*(F)](#ermf)
&nbsp;
## Class D 2 beta-lactamase
### Gather hits from fARGene results
```
mkdir beta_lactamase_d_2
cd beta_lactamase_d_2

# Gather ORFs
cat ../fARGene_d_2_out/predictedGenes/predicted-orfs.fasta > beta_lactamase_d_2.fasta

# Manually create lists for contig names with hits per sample to obtain the contigs
for i in $(less EFF1_lista.txt);do grep -A 1 -f <(echo "$i") ../../EFF1/EFF1_contigs.fasta > "EFF1_"$i".fasta";done
```
#### Preliminary exploration revealed that there are different types within the class D 2 beta-lactamases, we'll focus on those that were similar to the one in our bins

### Explore within all data
```
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p ID.txt)

blastn -query ../$sample/$sample"_contigs.fasta" \
        -subject beta_lactamase_d_2.fasta \
        -out $sample"_beta_lactamase_d_2_out.txt" -outfmt 6

# Filter hits
awk '$3 >= 90' EFF1_beta_lactamase_d_2_out.txt > filt_EFF1_beta_lactamase_d_2_out.txt

# Extract contigs
cut -f 1  filt_EFF1_beta_lactamase_d_2_out.txt > EFF1_beta_lactamase_d_2_lista.txt
for i in $(less EFF1_beta_lactamase_d_2_lista.txt);do grep -A 1 -f <(echo "$i") ../EFF1/EFF1_contigs.fasta > "EFF1_"$i".fasta";done
```
### Gather all sequences
```
# Concatenate sequences
cat *l.fasta > contigs_beta_lactamase_d_2.fasta
cat *c.fasta >> contigs_beta_lactamase_d_2.fasta
# Check lengths
module load seqkit/2.5.1
seqkit fx2tab --length --name --header-line contigs_beta_lactamase_d_2.fasta > contigs_beta_lactamase_d_2_lengths.txt
```
### Run Bakta for the contigs
```
ls *l.fasta | sed 's/\.fasta//g' > sample_names.txt
ls *c.fasta | sed 's/\.fasta//g' >> sample_names.txt

export SING_IMAGE=/projappl/project_2006608/containers/bakta:1.11.0.sif
export BAKTA_DB=/scratch/project_2006608/Methylation/db/bakta_db/db-light

apptainer_wrapper exec bakta $name".fasta" \
	--prefix $name \
	--output $name"_bakta_out"/ \
	--db $BAKTA_DB \
	--keep-contig-headers \
	--threads $SLURM_CPUS_PER_TASK
```
### Get ARG locations using blastn
```
module load biokit

cd WW_data/beta_lactamase_d_2/

blastn -query contigs_beta_lactamase_d_2.fasta \
        -subject beta_lactamase_d_2.fasta \
        -out blastn_out.txt -outfmt 6

# Create locations.txt according to blast hits
cut -f 1,7,8 blastn_out.txt > locations.txt

# Modify the duplicates manually
```
### Extract sequences
#### get_regions.sh
```
./get_regions.sh locations.txt
```
```
#!/bin/bash

touch startFlank.txt
touch endFlank.txt

while read -a line
do
  	contig=${line[0]}
        start=${line[1]}
        end=${line[2]}
        startFlank=$((${line[1]} - 5000))
        endFlank=$((${line[2]} + 5000))
        echo $startFlank >> startFlank.txt
        echo $endFlank >> endFlank.txt

        # Combine into one file
        paste -d '\t' locations.txt startFlank.txt endFlank.txt > locations_flanking.txt

done < $1
```
#### Polish results
```
# replace negative with one
less locations_flanking.txt | grep "-"

sed -i 's/-[0-9][0-9][0-9]/1/g' locations_flanking.txt
sed -i 's/-[0-9][0-9]/1/g' locations_flanking.txt

less locations_flanking.txt | grep "-"

rm startFlank.txt
rm endFlank.txt

mkdir extracted_data

# update sample_names.txt
cut -f 1 locations_flanking.txt > sample_names.txt
```
#### Run [extract_regions.sh](./../src/extract_regions.sh)
```
./extract_regions.sh
```
### Look for similar genes among the established ones
```
cd extracted_data

sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p ../sample_names.txt)

blastn -query $sample".fasta" \
        -subject db/resfinder_db/all.fsa \
        -out $sample"_resfinder_out.txt" -outfmt 6
```
### Bakta for extarcted sequences
```
cd extracted_data

# Set variable
name=$(sed -n ${SLURM_ARRAY_TASK_ID}p ../sample_names.txt)

# Run
export SING_IMAGE=/projappl/project_2006608/containers/bakta:1.11.0.sif
export BAKTA_DB=/scratch/project_2006608/Methylation/db/bakta_db/db-light

apptainer_wrapper exec bakta $name".fasta" \
        --prefix $name \
        --output $name"_bakta_out"/ \
        --db $BAKTA_DB \
        --keep-contig-headers \
        --threads $SLURM_CPUS_PER_TASK
```
### vsearch to dereplicate the contexts
```
cd extracted_data
cat *fasta > beta_lactamase_d_2_clustered.fasta
vsearch --cluster_fast beta_lactamase_d_2_clustered.fasta --id 0.99 --centroids beta_lactamase_d_2_clustered_vsearch_99.fasta --strand both --sizeout --threads $SLURM_CPUS_PER_TASK

# Polish
sed -i 's/;/_/g' beta_lactamase_d_2_clustered_vsearch_99.fasta
sed -i 's/=/_/g' beta_lactamase_d_2_clustered_vsearch_99.fasta

# Extract sequences into individual fasta
mkdir vsearch_99
mv beta_lactamase_d_2_clustered_vsearch_99.fasta vsearch_99
cd vsearch_99
seqretsplit beta_lactamase_d_2_clustered_vsearch_99.fasta
```
### Add bin & reference sequences
- Run Bakta
- Identify and extract the ARG locations
#### Add the bin sequences eventhought they would have be excluded by the deprelication
| ID | bin |
| ------------- | ------------- |
| INF1_s0.ctg006697l | C5f |
| INF3_s20.ctg001312l | C4f |
| INF2_s1.ctg008129l | C7f |
| INF3_s10.ctg003382l | C7f |
#### Add also the reference sequences from NCBI based in their similarity to our sequences (context/ARG)
| Accession | Description |
| ------------- | ------------- |
| KU721147 | *bla*OXA-464 for rooting the tree, root.fna |
| KU721146 | *bla*OXA-490 for rooting the tree, root.fna |
| CP183406.1 | *Arcobacter butzleri* |
| CP053835.1 | *Arcobacter defluvii* |
| CP042652.1 | *Arcobacter articola* |
| CP032100.1 | *Arcobacter suis* |
### Create the phylogenetic tree
#### The ARG sequences are gathered to ```genes_bl_d_2.fasta```
```
# Root
module load mafft/7.505
mafft ../../root.fna > ../../root_aln

# Alignment
mafft genes_bl_d_2.fasta > genes_bl_d_2_aln
mafft --add genes_bl_d_2_aln --reorder ../../root_aln > updated_aln

# Create the tree
module load raxml/8.2.12
raxmlHPC-PTHREADS -T 6 -s updated_aln -m GTRGAMMA -p 12345 -n beta_lactamase_d_2_genes_tree_vsearch_99

# Visualize the tree in iTOL
```

### Visualize the genetic contexts using pyGenomeViz
#### Some of the sequences are in reverted orientation, run [revert.py](./../src/revert.py):
```
module load python-data
module load biopythontools
python3 revert.py
```
#### (From this on locally)
#### Update the gene names so that more information on the gene annotations are available for the sequence visualization by running [replace_gene_names.py](./../src/replace_gene_names.py)
```
python3 replace_gene_names.py
```
#### Create the figures by running [pyGenomeViz_BLAST_d_2.py](./../src/pyGenomeViz_BLAST_d_2.py)
```
python3 pyGenomeViz_BLAST_d_2.py
```

### Run geNomad to check the likelyhood of the sequences being plasmid-borne
```
# Run
apptainer exec --bind $PWD:$PWD,$DB_PATH:$DB_PATH \
        /projappl/project_2006608/containers/genomad:1.8.0.sif genomad end-to-end \
        --cleanup --splits 8 $contig".fasta" geNomad_out/$contig $DB_PATH

# Check
cat */*_summary/*_plasmid_summary.tsv | grep -v "seq_name"
```
&nbsp;
## Class C beta-lactamase
```
mkdir beta_lactamase_c
cd beta_lactamase_c

# Gather ORFs
cat ../fARGene_c_out/predictedGenes/predicted-orfs.fasta > beta_lactamase_c.fasta

# Manually create lists for contig names with hits per sample to obtain the contigs
for i in $(less INF1_lista.txt);do grep -A 1 -f <(echo "$i") ../../INF1/INF1_contigs.fasta > "INF1_"$i".fasta";done
```
#### Preliminary exploration revealed that there are different types within the class C beta-lactamases, we'll focus on those that were similar to the one in our bins
### Explore within all data
```
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p ID.txt)

blastn -query ../$sample/$sample"_contigs.fasta" \
        -subject beta_lactamase_c.fasta \
        -out $sample"_beta_lactamase_c_out.txt" -outfmt 6

# Filter results
awk '$4 >= 1000' INF1_contigs_blaMCA_out.txt > filt_INF1_blaMCA_out.txt

# Extract contigs
cut -f 1  filt_INF1_beta_lactamase_c_out.txt > INF1_beta_lactamase_c_lista.txt
for i in $(less INF1_beta_lactamase_c_lista.txt);do grep -A 1 -f <(echo "$i") ../INF1/INF1_contigs.fasta > "INF1_"$i".fasta";done
```
### Gather all sequences
```
# Concatenate sequences
cat *l.fasta > contigs_beta_lactamase_c.fasta
cat *c.fasta >> contigs_beta_lactamase_c.fasta
# Check lengths
module load seqkit/2.5.1
seqkit fx2tab --length --name --header-line contigs_beta_lactamase_c.fasta > contigs_beta_lactamase_c_lengths.txt
```
### Run Bakta for the contigs
```
ls *l.fasta | sed 's/\.fasta//g' > sample_names.txt
ls *c.fasta | sed 's/\.fasta//g' >> sample_names.txt

export SING_IMAGE=/projappl/project_2006608/containers/bakta:1.11.0.sif
export BAKTA_DB=/scratch/project_2006608/Methylation/db/bakta_db/db-light

apptainer_wrapper exec bakta $name".fasta" \
	--prefix $name \
	--output $name"_bakta_out"/ \
	--db $BAKTA_DB \
	--keep-contig-headers \
	--threads $SLURM_CPUS_PER_TASK
```
### Get ARG locations using blastn
```
module load biokit

cd WW_data/beta_lactamase_c/

blastn -query contigs_beta_lactamase_c.fasta \
        -subject beta_lactamase_c.fasta \
        -out blastn_out.txt -outfmt 6

# Create locations.txt according to blast hits
cut -f 1,7,8 blastn_out.txt > locations.txt

# Modify the duplicates manually
```
### Extract sequences
#### get_regions.sh
```
./get_regions.sh locations.txt
```
```
#!/bin/bash

touch startFlank.txt
touch endFlank.txt

while read -a line
do
  	contig=${line[0]}
        start=${line[1]}
        end=${line[2]}
        startFlank=$((${line[1]} - 5000))
        endFlank=$((${line[2]} + 5000))
        echo $startFlank >> startFlank.txt
        echo $endFlank >> endFlank.txt

        # Combine into one file
        paste -d '\t' locations.txt startFlank.txt endFlank.txt > locations_flanking.txt

done < $1
```
#### Polish results
```
# replace negative with one
less locations_flanking.txt | grep "-"

sed -i 's/-[0-9][0-9][0-9]/1/g' locations_flanking.txt
sed -i 's/-[0-9][0-9]/1/g' locations_flanking.txt

less locations_flanking.txt | grep "-"

rm startFlank.txt
rm endFlank.txt

mkdir extracted_data

# update sample_names.txt
cut -f 1 locations_flanking.txt > sample_names.txt
```
#### Run [extract_regions.sh](./../src/extract_regions.sh)
```
./extract_regions.sh
```
### Bakta for extarcted sequences
```
cd extracted_data

# Set variable
name=$(sed -n ${SLURM_ARRAY_TASK_ID}p ../sample_names.txt)

# Run
export SING_IMAGE=/projappl/project_2006608/containers/bakta:1.11.0.sif
export BAKTA_DB=/scratch/project_2006608/Methylation/db/bakta_db/db-light

apptainer_wrapper exec bakta $name".fasta" \
        --prefix $name \
        --output $name"_bakta_out"/ \
        --db $BAKTA_DB \
        --keep-contig-headers \
        --threads $SLURM_CPUS_PER_TASK
```
### vsearch to dereplicate the contexts
```
cd extracted_data
cat *fasta > beta_lactamase_c_clustered.fasta
vsearch --cluster_fast beta_lactamase_c_clustered.fasta --id 0.90 --centroids beta_lactamase_c_clustered_vsearch_90.fasta --strand both --sizeout --threads $SLURM_CPUS_PER_TASK

# Polish
sed -i 's/;/_/g' beta_lactamase_c_clustered_vsearch_90.fasta
sed -i 's/=/_/g' beta_lactamase_c_clustered_vsearch_90.fasta

# Extract sequences into individual fasta
mkdir vsearch_90
mv beta_lactamase_c_clustered_vsearch_90.fasta vsearch_90
cd vsearch_90
module load biokit
seqretsplit beta_lactamase_c_clustered_vsearch_90.fasta
```
### Add reference sequences
- Run Bakta
- Identify and extract the ARG locations
| Accession | Description |
| ------------- | ------------- |
| MRSN15084 | Acinetobacter baumannii |
### Create the phylogenetic tree
#### The ARG sequences are gathered to ```genes_bl_c.fasta```
```
# Alignment
mafft genes_bl_c.fasta > genes_bl_c_aln

# Create the tree
module load raxml/8.2.12
raxmlHPC-PTHREADS -T 6 -s genes_bl_c_aln -m GTRGAMMA -p 12345 -n beta_lactamase_c_genes_tree_vsearch_90

# Visualize the tree in iTOL
```
### Visualize the genetic contexts using pyGenomeViz
#### Some of the sequences are in reverted orientation, run [revert.py](./../src/revert.py):
```
module load python-data
module load biopythontools
python3 revert.py
```
#### (From this on locally)
#### Update the gene names so that more information on the gene annotations are available for the sequence visualization by running [replace_gene_names.py](./../src/replace_gene_names.py)
```
python3 replace_gene_names.py
```
#### Create the figures by running [pyGenomeViz_BLAST_c.py](./../src/pyGenomeViz_BLAST_c.py)
```
python3 pyGenomeViz_BLAST_c.py
```
### Run geNomad to check the likelyhood of the sequences being plasmid-borne
```
# Run
apptainer exec --bind $PWD:$PWD,$DB_PATH:$DB_PATH \
        /projappl/project_2006608/containers/genomad:1.8.0.sif genomad end-to-end \
        --cleanup --splits 8 $contig".fasta" geNomad_out/$contig $DB_PATH

# Check
cat */*_summary/*_plasmid_summary.tsv | grep -v "seq_name"
```
### pdifFinder
```
# Add the program to your path
export PYTHONUSERBASE=/projappl/project_2006608/my-python-env
export PATH="/projappl/project_2006608/pdifFinder/pdifFinder_installation/bin:$PATH"

# Run
/projappl/project_2006608/my-python-env/bin/pdifFinder --inFile $name".fasta" --outdir $name"_pdifFinder_out"

# View
ls *pdifFinder_out/*
tail *_pdifFinder_out/*/pdif_site.txt
```
## *bla*OXA-129
```
mkdir blaOXA-129
cd blaOXA-129
for i in $(less INF1_lista.txt);do grep -A 1 -f <(echo "$i") ../INF1/INF1_contigs.fasta > "INF1_"$i".fasta";done
```
### Add reference sequences
| Accession | Description |
| ------------- | ------------- |
| FJWZ01000025.1 | Enterobacter hormaechei |
| NZ_AP022010.1 | Klebsiella quasipneumoniae |
| NZ_CP031449.2 | Pseudomonas aeruginosa |

### Run Bakta
```
ls I*fasta | sed 's/\.fasta//g' > contig_names.txt
ls N*fasta | sed 's/\.fasta//g' >> contig_names.txt

contig=$(sed -n ${SLURM_ARRAY_TASK_ID}p contig_names.txt)

# Run
export SING_IMAGE=/projappl/project_2006608/containers/bakta:1.11.0.sif
export BAKTA_DB=/scratch/project_2006608/Methylation/db/bakta_db/db-light

apptainer_wrapper exec bakta $contig".fasta" \
        --prefix $contig \
        --output $contig"_bakta_out" \
        --db $BAKTA_DB \
        --keep-contig-headers \
        --threads $SLURM_CPUS_PER_TASK
```
### Extract flanking regions
```
# Check locations
cat *_bakta_out/*.tsv | grep "OXA" > locations.txt

awk '{print $1,$3,$4}' locations.txt > tmp && mv tmp locations.txt
sed -i 's/ /\t/g' locations.txt
```
```
./get_regions.sh locations.txt
```
```
#!/bin/bash

touch startFlank.txt
touch endFlank.txt

while read -a line
do
  	contig=${line[0]}
        start=${line[1]}
        end=${line[2]}
        startFlank=$((${line[1]} - 10000))
        endFlank=$((${line[2]} + 10000))
        echo $startFlank >> startFlank.txt
        echo $endFlank >> endFlank.txt

        # Combine into one file
        paste -d '\t' locations.txt startFlank.txt endFlank.txt > locations_flanking.txt

done < $1
```
#### Polish results
```
# replace negative with one
less locations_flanking.txt | grep "-"

sed -i 's/-[0-9][0-9][0-9]/1/g' locations_flanking.txt
sed -i 's/-[0-9][0-9]/1/g' locations_flanking.txt

less locations_flanking.txt | grep "-"

rm startFlank.txt
rm endFlank.txt

mkdir extracted_data

# update sample_names.txt
cut -f 1 locations_flanking.txt > sample_names.txt

mkdir extracted_data
```
#### Run [extract_regions.sh](./../src/extract_regions.sh)
```
./extract_regions.sh
```
### Run bakta on the extracted squences
```
export SING_IMAGE=/projappl/project_2006608/containers/bakta:1.11.0.sif
export BAKTA_DB=/scratch/project_2006608/Methylation/db/bakta_db/db-light

apptainer_wrapper exec bakta $name".fasta" \
        --prefix $name \
        --output $name"_bakta_out"/ \
        --db $BAKTA_DB \
        --keep-contig-headers \
        --threads $SLURM_CPUS_PER_TASK
```

### Extract the genes
```
cd WW_data/blaOXA-129/extracted_data
cat *fasta > contigs_blaOXA-129.fasta

# Load the tools
module load biokit

blastn -query contigs_blaOXA-129.fasta \
        -subject ../blaOXA-129.fasta \
        -out blastn_out.txt -outfmt 6 -max_target_seqs 1
```
```
./extract_gene.sh
```
```
#!/bin/bash

# Load tools
module load seqkit/2.5.1

accessions=$(less ../sample_names.txt)

for a in $accessions;
do
        line=$(grep "$a" blastn_out.txt)
        if [[ -n "$line" ]]; then
                startFlank=$(echo $line | cut -d' ' -f 7)
                endFlank=$(echo $line | cut -d' ' -f 8)
                seqkit subseq -r $startFlank:$endFlank *$a".fasta" > $a"_gene_blaOXA-129_"$startFlank"_"$endFlank".fasta"
        fi
done
```
### Visualize the genetic contexts using pyGenomeViz
#### Some of the sequences are in reverted orientation, run [revert.py](./../src/revert.py):
```
module load python-data
module load biopythontools
python3 revert.py
```
#### (From this on locally)
#### Update the gene names so that more information on the gene annotations are available for the sequence visualization by running [replace_gene_names.py](./../src/replace_gene_names.py)
```
python3 replace_gene_names.py
```
#### Create the figures by running [pyGenomeViz_BLAST_blaOXA-129.py](./../src/pyGenomeViz_BLAST_blaOXA-129.py)
```
python3 pyGenomeViz_BLAST_blaOXA-129.py
```
## *sul1_9*-like
```
mkdir sul1_9
cd sul1_9
```
### Gather contigs with mach to sul1_9
```
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p ID.txt)

blastn -query ../$sample/$sample"_contigs.fasta" \
        -subject sul1_9.fasta \
        -out $sample"_sul1_9.txt" -outfmt 6

# Filter only the best hits to avoid getting much of the other variants of sul1

# Extract contigs
for i in $(less EFF1_lista.txt);do grep -A 1 -f <(echo "$i") ../../EFF1/EFF1_contigs.fasta > $i".fasta";done
```
### Add reference sequences
| Accession | Description |
| ------------- | ------------- |
| CP026207.1 | *Escherichia coli* |
| CP055486.1 | *Klebsiella* sp. |
| CP021775.1 | *Pseudomomas aeruginosa* |
### Get the locations of the ARGs
```
# Run
blastn -query $sample".fasta" \
        -subject sul1_9.fasta \
        -out $sample"_resfinder_out.txt" -outfmt 6 -max_target_seqs 1

# Create the locations.txt
cat *_resfinder_out.txt | cut -f 1,7-8 > locations.txt
```
```
./get_regions.sh locations.txt
```
```
#!/bin/bash

touch startFlank.txt
touch endFlank.txt

while read -a line
do
  	contig=${line[0]}
        start=${line[1]}
        end=${line[2]}
        startFlank=$((${line[1]} - 10000))
        endFlank=$((${line[2]} + 10000))
        echo $startFlank >> startFlank.txt
        echo $endFlank >> endFlank.txt

        # Combine into one file
        paste -d '\t' locations.txt startFlank.txt endFlank.txt > locations_flanking.txt

done < $1
```
#### Polish results
```
# replace negative with one
less locations_flanking.txt | grep "-"

sed -i 's/-[0-9][0-9][0-9]/1/g' locations_flanking.txt
sed -i 's/-[0-9][0-9]/1/g' locations_flanking.txt

less locations_flanking.txt | grep "-"

rm startFlank.txt
rm endFlank.txt

mkdir extracted_data

# update sample_names.txt
cut -f 1 locations_flanking.txt > sample_names.txt
```
#### Run [extract_regions.sh](./../src/extract_regions.sh)
```
./extract_regions.sh
```
### Run Bakta for the extracted data
```
cd extracted_data

apptainer_wrapper exec bakta s399.ctg000560l.fasta \
        --prefix s399.ctg000560l \
        --output s399.ctg000560l_bakta_out \
        --db $BAKTA_DB \
        --keep-contig-headers \
        --threads $SLURM_CPUS_PER_TASK
```
### vsearch to dereplicate the contexts
```
cd extracted_data
cat s*fasta > sul1_clustered.fasta

vsearch --cluster_fast sul1_clustered.fasta --id 0.90 --centroids sul1_clustered_vsearch_90.fasta --strand both --sizeout --threads $SLURM_CPUS_PER_TASK

# Polish
sed -i 's/;/_/g' sul1_clustered_vsearch_90.fasta
sed -i 's/=/_/g' sul1_clustered_vsearch_90.fasta
sed -i 's/_size_[0-9]//g' sul1_clustered_vsearch_90.fasta
sed -i 's/_size_[0-9]//g' sul1_clustered_vsearch_90.fasta
```
### Create phylogenetic tree
```
# Alignment
mafft genes_bl_c.fasta > genes_sul1_9_aln

# Create the tree
module load raxml/8.2.12
raxmlHPC-PTHREADS -T 6 -s genes_sul1_9_aln -m GTRGAMMA -p 12345 -n sul1_9_genes_tree_vsearch_90

# Visualize the tree in iTOL
```
### Visualize the genetic contexts using pyGenomeViz
#### Some of the sequences are in reverted orientation, run [revert.py](./../src/revert.py):
```
module load python-data
module load biopythontools
python3 revert.py
```
#### (From this on locally)
#### Update the gene names so that more information on the gene annotations are available for the sequence visualization by running [replace_gene_names.py](./../src/replace_gene_names.py)
```
python3 replace_gene_names.py
```
#### Create the figures by running [pyGenomeViz_BLAST_sul1_9.py](./../src/pyGenomeViz_BLAST_sul1_9.py)
```
python3 pyGenomeViz_BLAST_sul1_9.py
```
## *erm*(F)
```
mkdir erm_F
cd erm_F
```
### Gather the contigs with hits to erm(F)_1, erm(F)_3 or erm(F)_4
```
blastn -query $sample".fasta" \
        -subject erm_F.fasta \
        -out $sample"_resfinder_out.txt" -outfmt 6 -max_target_seqs 1

for i in $(less EFF1_lista.txt);do grep -A 1 -f <(echo "$i") ../EFF1/EFF1_contigs.fasta > "EFF1_"$i".fasta";done
```
### Run Bakta
```
ls *fasta | sed 's/\.fasta//g' > contig_names.txt

contig=$(sed -n ${SLURM_ARRAY_TASK_ID}p contig_names.txt)

# Run
export SING_IMAGE=/projappl/project_2006608/containers/bakta:1.11.0.sif
export BAKTA_DB=/scratch/project_2006608/Methylation/db/bakta_db/db-light

apptainer_wrapper exec bakta $contig".fasta" \
        --prefix $contig \
        --output $contig"_bakta_out" \
        --db $BAKTA_DB \
        --keep-contig-headers \
        --threads $SLURM_CPUS_PER_TASK
```

### Check locations
```
cat *_bakta_out/*l.tsv | grep "erm(F)" > locations.txt
awk '{print $1,$3,$4}' locations.txt > tmp && mv tmp locations.txt
sed -i 's/ /\t/g' locations.txt

# And and those identified only by BLAST
cat *resf* | cut -f 1,7-8 >> locations.txt

# Also some hits and their locations has to be fixed manually
```
### Extract flanking regions
#### Get the locations
```
./get_regions.sh locations.txt
```
```
#!/bin/bash

touch startFlank.txt
touch endFlank.txt

while read -a line
do
  	contig=${line[0]}
        start=${line[1]}
        end=${line[2]}
        startFlank=$((${line[1]} - 5000))
        endFlank=$((${line[2]} + 5000))
        echo $startFlank >> startFlank.txt
        echo $endFlank >> endFlank.txt

        # Combine into one file
        paste -d '\t' locations.txt startFlank.txt endFlank.txt > locations_flanking.txt

done < $1
```
#### Polish the resulys
```
# replace negative with one
less locations_flanking.txt | grep "-"

sed -i 's/-[0-9][0-9][0-9]/1/g' locations_flanking.txt
sed -i 's/-[0-9][0-9]/1/g' locations_flanking.txt

less locations_flanking.txt | grep "-"

rm startFlank.txt
rm endFlank.txt


mkdir extracted_data

# update sample_names.txt
cut -f 1 locations_flanking.txt > sample_names.txt
```
#### Run [extract_regions.sh](./../src/extract_regions.sh)
```
./extract_regions.sh
```
### Extract genes
```
cd extracted_data
cat *fasta > contigs_erm_f.fasta

# Load the tools
module load biokit

blastn -query contigs_erm_f.fasta \
        -subject ../erm_F.fasta \
        -out blastn_out.txt -outfmt 6 -max_target_seqs 1
```
```
./extract_gene.sh
```
```
#!/bin/bash

# Load tools
module load seqkit/2.5.1

accessions=$(less ../sample_names.txt)

for a in $accessions;
do
        line=$(grep "$a" blastn_out.txt)
        if [[ -n "$line" ]]; then
                startFlank=$(echo $line | cut -d' ' -f 7)
                endFlank=$(echo $line | cut -d' ' -f 8)
                seqkit subseq -r $startFlank:$endFlank *$a".fasta" > $a"_gene_erm_f_"$startFlank"_"$endFlank".fasta"
        fi
done
```
#### Wrangle
```
cat *gene_erm_f* > genes_erm_f.fasta
# Filter
awk '$2 >= 700' genes_erm_f_lengths.txt > filt_genes_erm_f_lengths.txt

# Extract contigs
cut -f 1 filt_genes_erm_f_lengths.txt > erm_f_lista.txt

# Remove extra space after gene name
sed -i 's/ //g' genes_erm_f.fasta
sed -i 's/ //g' erm_f_lista.txt
seqkit grep -f erm_f_lista.txt genes_erm_f.fasta > filt_genes_erm_f.fasta

cut -f 1 genes_erm_f_lengths.txt > ID.txt
cut -f 1 filt_genes_erm_f_lengths.txt > filt_ID.txt
sed -i 's/ //g' filt_ID.txt
grep -v -f filt_ID.txt ID.txt > remove_ID.txt

# Remove them
grep -Ff remove_ID.txt <(ls)
grep -Ff remove_ID.txt <(ls) | xargs -d '\n' rm -v

# Also from contigs_erm_f.fasta
seqkit grep -f remove_ID.txt -v contigs_erm_f.fasta > filt_contigs_erm_f.fasta
```
### Bakta for extracted
```
apptainer_wrapper exec bakta $name".fasta" \
        --prefix $name \
        --output $name"_bakta_out"/ \
        --db $BAKTA_DB \
        --keep-contig-headers \
        --threads $SLURM_CPUS_PER_TASK
```
### vsearch to reduce the number of sequences
```
cd extracted_data
cat *l.fasta > erm_f_clustered.fasta
vsearch --cluster_fast erm_f_clustered.fasta --id 0.50 --centroids erm_f_clustered_vsearch_50.fasta --strand both --sizeout --threads $SLURM_CPUS_PER_TASK

sed -i 's/;/_/g' erm_f_clustered_vsearch_50.fasta
sed -i 's/=/_/g' erm_f_clustered_vsearch_50.fasta

mkdir vsearch_50
nano vsearch_50/vsearch_50_list.txt
```
### Filter the genes in ```genes_erm_f.fasta``` according to the dereplication results
```
seqkit grep -f vsearch_50_list.txt ../genes_erm_f.fasta > genes_erm_f_vsearch_50.fasta
```
### Add reference sequences
| Accession | Description |
| ------------- | ------------- |
| M17808.1 | *Bacteroides fragilis* |
| CP054002.1 | *Bacteroides fragilis* |
### Build phylogenetic tree (no root?)
```
module load mafft/7.505
mafft genes_erm_f_vsearch_50.fasta > genes_erm_f_vsearch_50_aln
module load raxml/8.2.12
raxmlHPC-PTHREADS -T 6 -s genes_erm_f_vsearch_50_aln -m GTRGAMMA -p 12345 -n erm_f_genes_vsearch_50_tree
```
### Visualize the genetic contexts using pyGenomeViz
#### Some of the sequences are in reverted orientation, run [revert.py](./../src/revert.py):
```
module load python-data
module load biopythontools
python3 revert.py
```
#### (From this on locally)
#### Update the gene names so that more information on the gene annotations are available for the sequence visualization by running [replace_gene_names.py](./../src/replace_gene_names.py)
```
python3 replace_gene_names.py
```
#### Create the figures by running [pyGenomeViz_BLAST_ermF.py](./../src/pyGenomeViz_BLAST_ermF.py)
```
python3 pyGenomeViz_BLAST_ermF.py
```

