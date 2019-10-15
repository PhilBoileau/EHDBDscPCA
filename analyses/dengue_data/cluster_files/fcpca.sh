#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --ntasks=1
#SBATCH --nodes=1

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
R CMD BATCH --no-save fcpca.R fcpca.Rout