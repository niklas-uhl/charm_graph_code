#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --output=weighted_no_projections.out
#SBATCH --partition=cpu
#SBATCH --account=mzu-delta-cpu
#SBATCH --job-name=myjobtest
#SBATCH --exclusive
#SBATCH --time=00:10:00      # hh:mm:ss for the job

srun -n $SLURM_NTASKS ./weighted_htram_nonsmp $SLURM_NTASKS 100000 ../facebook_clean_data/big_graph.csv 100 1
