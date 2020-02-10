#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --nodes=1

R CMD BATCH --no-save rt_analysis.R rt_analysis.Rout
