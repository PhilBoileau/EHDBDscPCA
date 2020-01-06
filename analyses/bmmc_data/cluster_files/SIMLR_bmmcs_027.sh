#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1

R CMD BATCH --no-save SIMLR_bmmcs_027.R SIMLR_bmmcs_027.Rout
