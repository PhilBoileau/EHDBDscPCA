#!/bin/bash
#SBATCH --cpus-per-task=30
#SBATCH --ntasks=1
#SBATCH --nodes=1

R CMD BATCH --no-save zinbwave_035.R zinbwave_035.Rout
