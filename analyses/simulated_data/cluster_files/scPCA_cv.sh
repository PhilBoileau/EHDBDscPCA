#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1

R CMD BATCH --no-save scPCA_cv.R scPCA_cv.Rout
