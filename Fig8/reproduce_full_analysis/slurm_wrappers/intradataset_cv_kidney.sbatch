#!/bin/bash
#SBATCH --partition=gilad
#SBATCH --time=25:00:00
#SBATCH --mem=250G
#SBATCH --output=logs/cross_val_kid.out
#SBATCH --error=logs/cross_val_kid.err
#SBATCH --job-name=svkid

module load -f gdal/2.4.1
module load R/4.2.0

Rscript ../code/intradataset_cv_kidney.R
