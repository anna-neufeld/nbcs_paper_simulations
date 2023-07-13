#!/bin/bash
#SBATCH --partition=gilad
#SBATCH --time=10:00
#SBATCH --mem=250G
#SBATCH --output=logs/subset_counts.out
#SBATCH --error=logs/subset_counts.err
#SBATCH --job-name=subset

module load -f gdal/2.4.1
module load R/4.2.0

Rscript ../code/subset_counts.R
