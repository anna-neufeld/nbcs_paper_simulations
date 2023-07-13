#!/bin/bash
#SBATCH --partition=gilad
#SBATCH --time=05:00:00
#SBATCH --mem=250G
#SBATCH --output=logs/create_cds.out
#SBATCH --error=logs/create_cds.err
#SBATCH --job-name=cdsmk

module load -f gdal/2.4.1
module load R/4.2.0

Rscript ../code/create_cds.R "$@"

