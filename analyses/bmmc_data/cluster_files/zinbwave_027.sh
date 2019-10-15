#!/bin/bash
#SBATCH --cpus-per-task=30
#SBATCH --ntasks=1
#SBATCH --nodes=1

R CMD BATCH --no-save zinbwave_027.R zinbwave_027.Rout
