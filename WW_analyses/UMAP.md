# Generating attached data
## Contig lengths
```
module load seqkit/2.5.1
seqkit fx2tab --length --name --header-line EFF1_contigs.fasta > EFF1_contigs_lengths.txt
```
## Base modification counts
```
cd src/
./WW_count_gff_lines.sh EFF1
sed -i '1i contig\tmod_count' EFF1_mod_counts.txt
```
## Established ARGs (WW_blastn_resfinder.sh)
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
### Filter & count ARGs / contig
```
awk '$4 >= 100' EFF1_resfinder_out.txt" > tmp && mv tmp EFF1_resfinder_out.txt

cp EFF1_resfinder_out.txt raw_EFF1_resfinder_out.txt

## First manually (set # and then remove those lines) check hits that are duplicates

sed -i '/^#/d' EFF1_resfinder_out.txt

# Then run
cd src/
./WW_count_ARGs.sh EFF1
```
### Attach ARG names
```
# Extract columns
cut -f 1-2 EFF1_resfinder_out.txt > ARG_names.txt
# Remove the accession
sed -i 's/\(.*_.*\)_.*$/\1/' ARG_names.txt
# Check those that have multiple
cut -f 1 ARG_names.txt | sort | uniq -d

## Manually fix those contigs that have multiple

# Add contig names of those where there is no ARG
cut -f 1 ARG_names.txt > contig_ARG_names.txt
grep -v -f contig_ARG_names.txt contig_names.txt > remaining_contig_names.txt
cat ARG_names.txt remaining_contig_names.txt > EFF1_ARG_names.txt
sed -i '1i contig\tARG_name' EFF1_ARG_names.txt
```
## Latent ARGs
```
fargene -i fARGene_in/*.fasta \
        --hmm-model db/fargene/fargene_analysis/models/beta_lactamase_d_2.hmm \
        -o fARGene_d_2_out -p $SLURM_CPUS_PER_TASK --score 200
```

# Binning clustered contigs
## Extract clusters
### ./WW_get_contigs_for_MAGs.sh <sample> <cluster> <above>
```
./WW_get_contigs_for_MAGs.sh EFF1 C1
```