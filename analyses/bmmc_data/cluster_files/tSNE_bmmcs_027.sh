#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1

R CMD BATCH --no-save tSNE_bmmcs_027.R tSNE_bmmcs_027.Rout
