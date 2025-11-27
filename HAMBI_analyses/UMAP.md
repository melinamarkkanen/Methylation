# Uniform Manifold Approximation and Projection visualization for the synthetic community data
&nbsp;
## Generating attached data
### Contig lengths
```
cd HAMBI_data/metagenomic_assembly
module load seqkit/2.5.1
seqkit fx2tab --length --name --header-line $sample".fasta" > $sample"_lengths.txt"
```
### Base modification counts
#### Run ```src/HAMBI_count_gff_lines.sh```
##### Combine contigs_lengths.tsv & contigs_mod_counts.tsv & ARG names & counts to merged_data.tsv within ```notebooks/UMAP_analysis_HAMBI.ipynb```
&nbsp;
## Binning clustered contigs
### Gather contigs based n coordinates in UMAP (```notebooks/UMAP_analysis_HAMBI.ipynb```)
```
cd ../src
./HAMBI_get_contigs_for_MAGs.sh C1
```
### Combine contigs sample-wise
```
cd /HAMBI_data/MAGs
# Check samples
ls */*contigs/
# Check contig lengths
tail -n 200 C1*/C1*txt
# Combine
## C1
cat C1/bcAd1023T--bcAd1023T_contigs/*fasta > C1/bcAd1023T--bcAd1023T_contigs/bcAd1023T--bcAd1023T_C1.fasta
```
### Run CheckM2
```
# Create sample list
cd HAMBI_data/MAGs/C1
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
cat HAMBI_data/MAGs/C1/*_CheckM2_out/quality_report.tsv | cut -f 1-3 | grep -v "Name"
```
### Run GTDB-Tk
```
# Go to folder
cd HAMBI_data/MAGs/C1

# Set variable
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_names.txt)

# Load the environment and variables
export PATH="/projappl/project_2006608/GTDB-Tk/bin:$PATH"
export GTDBTK_DATA_PATH=/scratch/project_2006608/GTDB-Tk/release220/

# Run
gtdbtk classify_wf --genome_dir $sample"_contigs" -x fasta \
        --out_dir $sample"_GTDB_out" --skip_ani_screen --cpus $SLURM_CPUS_PER_TASK

# Summarize
cat *GTDB_out/g*tsv | cut -f 1,5
```