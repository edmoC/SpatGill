#!/bin/bash
#SBATCH --job-name=mvmt_pairs
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=50

module load r/4.0.0

Rscript --vanilla mvmt_pairs.R
