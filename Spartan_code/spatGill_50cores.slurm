#!/bin/bash
#SBATCH --job-name=spatGill
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=50

module load r/4.0.0

Rscript --vanilla spatGill.R
