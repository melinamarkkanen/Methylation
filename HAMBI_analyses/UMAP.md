# Generating attached data
## Contig lengths
## Base modification counts
## Established ARGs
## Latent ARGs
# Binning clustered contigs





### Run ```notebooks/UMAP_analysis_HAMBI.ipynb```
&nbsp;
#### For visualization, obtain contig lengths and modification data amounts
##### Contig lengths
```
cd /scratch/project_2006608/Methylation/HAMBI_data/metagenomic_assembly
module load seqkit/2.5.1
seqkit fx2tab --length --name --header-line $sample".fasta" > $sample"_lengths.txt"
```
##### Number of modified bases, run ```src/HAMBI_count_gff_lines.sh```
##### Merge contigs_lengths.tsv & contigs_mod_counts.tsv with merged_data.tsv within ```notebooks/UMAP_analysis_HAMBI.ipynb```
&nbsp;

### Build MAGs according to contigs clustered by UMAP
#### Get contigs
```
cd ../src
./HAMBI_get_contigs_for_MAGs.sh C1A
```
&nbsp;
&nbsp;
&nbsp;
&nbsp;

#### Combine small contigs sample-wise
```
cd /scratch/project_2006608/Methylation/HAMBI_data/MAGs
# Check samples
ls */*contigs/
# Check contig lengths
tail -n 200 C10*/C10*txt
# Combine
## C10
cat C10/bcAd1023T--bcAd1023T_contigs/*fasta > C10/bcAd1023T--bcAd1023T_contigs/bcAd1023T--bcAd1023T_C10.fasta
cat C10/bcAd1039T--bcAd1039T_contigs/*fasta > C10/bcAd1039T--bcAd1039T_contigs/bcAd1039T--bcAd1039T_C10.fasta

# Remove these:
bcAd1023T--bcAd1023T_ptg000007c
bcAd1023T--bcAd1023T_ptg000008l
bcAd1037T--bcAd1037T_ptg000008c
bcAd1039T--bcAd1039T_ptg000024c
bcAd1046T--bcAd1046T_ptg000003c
bcAd1063T--bcAd1063T_ptg000008c

seqkit grep -v -n -f list.txt bcAd1023T--bcAd1023T_contigs/all.fasta > tmp && mv tmp bcAd1023T--bcAd1023T_contigs/all.fasta
seqkit grep -v -n -f list.txt bcAd1037T--bcAd1037T_contigs/all.fasta > tmp && mv tmp bcAd1037T--bcAd1037T_contigs/all.fasta
seqkit grep -v -n -f list.txt bcAd1039T--bcAd1039T_contigs/all.fasta > tmp && mv tmp bcAd1039T--bcAd1039T_contigs/all.fasta
seqkit grep -v -n -f list.txt bcAd1046T--bcAd1046T_contigs/all.fasta > tmp && mv tmp bcAd1046T--bcAd1046T_contigs/all.fasta
seqkit grep -v -n -f list.txt bcAd1063T--bcAd1063T_contigs/all.fasta > tmp && mv tmp bcAd1063T--bcAd1063T_contigs/all.fasta
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
cat HAMBI_data/MAGs/C10/*_CheckM2_out/quality_report.tsv | cut -f 1-3 | grep -v "Name"
```

#### Check if plasmid or ARGs
```
cat HAMBI_labels.txt | grep -f C10_id.txt | grep -v "chromoso" | cut -f 1,9
cat *res* | grep -f C10_id.txt | cut -f 1-4
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


# Same for data with all features
# Extract those that are in HAMBI_data/metagenomic_assembly/contig_IDs.txt
cd /scratch/project_2006608/Methylation/HAMBI_data/metagenomic_assembly
cat ../../HAMBI_data/contigs/HAMBI_mod_counts.txt | grep -f contig_IDs.txt > contigs_mod_counts.txt

# Reorder to match
awk 'NR==FNR {order[$1]=NR; next} {print order[$1], $0}' contig_IDs.txt contigs_mod_counts.txt | sort -n | cut -d' ' -f2- > tmp
tr ' ' '\t' < tmp > contigs_mod_counts.tsv
sed -i '1s/^/contig\tmod_count\n/' contigs_mod_counts.tsv
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