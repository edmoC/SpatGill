#!/bin/bash
#SBATCH --job-name=mvmt
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=50

module load r/4.0.0

Rscript --vanilla mvmt.R
