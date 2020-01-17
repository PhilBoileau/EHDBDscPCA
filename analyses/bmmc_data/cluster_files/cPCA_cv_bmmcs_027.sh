#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --nodes=1

R CMD BATCH --no-save cPCA_cv_bmmcs_027.R cPCA_cv_bmmcs_027.Rout
