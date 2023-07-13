#!/bin/bash
#SBATCH --partition=gilad
#SBATCH --time=02:00:00
#SBATCH --mem=250G
#SBATCH --output=logs/traintest_comparison_ari.out
#SBATCH --error=logs/traintest_comparison_ari.err
#SBATCH --job-name=ariall

module load -f gdal/2.4.1
module load R/4.2.0

Rscript ../code/traintest_comparison_ari.R "$@"
