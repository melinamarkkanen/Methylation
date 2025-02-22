#!/bin/bash
#SBATCH --job-name=snakemake_hq
#SBATCH --account=project_2002265
#SBATCH --error=slurm_logs/snakemake_hq_%j_err.txt
#SBATCH --output=slurm_logs/snakemake_hq_%j_out.txt
#SBATCH --time=16:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks-per-node=1
#SBATCH --partition=small

module load hyperqueue
module load snakemake/8.4.6

# Create a per job directory
export HQ_SERVER_DIR=$PWD/.hq-server-$SLURM_JOB_ID
mkdir -p $HQ_SERVER_DIR

hq server start &
srun --cpu-bind=none --hint=nomultithread --mpi=none -N $SLURM_NNODES -n $SLURM_NNODES -c 40 hq worker start --cpus=40 &

num_up=$(hq worker list | grep RUNNING | wc -l)
while true; do

    echo "Checking if workers have started"
    if [[ $num_up -eq $SLURM_NNODES ]];then
        echo "Workers started"
        break
    fi
    echo "$num_up/$SLURM_NNODES workers have started"
    sleep 1
    num_up=$(hq worker list | grep RUNNING | wc -l)

done

snakemake -s workflow/Snakefile_WW_methylation_analysis -j 100 --rerun-incomplete \
        --use-singularity --use-envmodules --executor cluster-generic --cluster-generic-submit-cmd "hq submit --cpus 2"

hq worker stop all
hq server stop