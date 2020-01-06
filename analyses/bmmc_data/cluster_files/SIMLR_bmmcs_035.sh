#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --nodes=1

R CMD BATCH --no-save SIMLR_bmmcs_035.R SIMLR_bmmcs_035.Rout
