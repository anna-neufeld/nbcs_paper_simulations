#!/bin/bash
#SBATCH --partition=gilad
#SBATCH --time=05:00:00
#SBATCH --mem=450G
#SBATCH --output=logs/split.out
#SBATCH --error=logs/split.err
#SBATCH --job-name=split

module load R/4.2.0

Rscript ../code/countsplit.R "$@"
