#!/bin/bash
#SBATCH --partition=gilad
#SBATCH --time=05:00:00
#SBATCH --mem=250G
#SBATCH --output=logs/cross_val_meta.out
#SBATCH --error=logs/cross_val_meta.err
#SBATCH --job-name=svmeta

module load -f gdal/2.4.1
module load R/4.2.0

Rscript ../code/intradataset_cv_metanephric.R

