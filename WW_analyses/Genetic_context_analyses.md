# Analysis steps for the ARG genetic contexts of selected ARG bins
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
#### extract_regions.sh
```
./extract_regions.sh
```
```
#!/bin/bash

# Load tools
module load seqkit/2.5.1

accessions=$(less sample_names.txt)

for a in $accessions;
do
        line=$(grep "$a" locations_flanking.txt)
        if [[ -n "$line" ]]; then
                startFlank=$(echo $line | cut -d' ' -f 4)
                endFlank=$(echo $line | cut -d' ' -f 5)
                seqkit subseq -r $startFlank:$endFlank *$a"_bakta_out"/*$a".fna" > extracted_data/$a".fasta"
        fi
done
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
```
INF1_s0.ctg006697l C5f
INF3_s20.ctg001312l C4f
INF2_s1.ctg008129l C7f
INF3_s10.ctg003382l C7f
```
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
#### Some of the sequences are in reverted orientation, run ```revert.py```:
```
module load python-data
module load biopythontools
python3 revert.py
```
#### (From this on locally)
#### Update the gene names so that more information on the gene annotations are available for the sequence visualization by running ```replace_gene_names.py```
```
python3 replace_gene_names.py
```
#### Create the figures by running ```pyGenomeViz_BLAST_d_2.py```
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





## Class C beta-lactamase ### pdifFinder
## *bla*OXA-129
## *sul1*
## *erm(F)*
