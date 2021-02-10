#!/bin/bash
#SBATCH --job-name=h_mvmt
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=50
#SBATCH --array=3-20

module load r/4.0.0

Rscript --vanilla h_mvmt.R ${SLURM_ARRAY_TASK_ID}