#!/bin/bash
#SBATCH --cpus-per-task=30
#SBATCH --ntasks=1
#SBATCH --nodes=1

R CMD BATCH --no-save scPCA_bmmcs_035.R scPCA_bmmcs_035.Rout
