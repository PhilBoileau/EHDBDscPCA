#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1

R CMD BATCH --no-save cpca_cv.R cpca_cv.Rout