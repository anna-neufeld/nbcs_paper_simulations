#!/bin/bash
#SBATCH --partition=gilad
#SBATCH --time=02:00:00
#SBATCH --mem=250G
#SBATCH --output=logs/cluster_cds.out
#SBATCH --error=logs/cluster_cds.err
#SBATCH --job-name=clust

module load -f gdal/2.4.1
module load R/4.2.0

Rscript ../code/cluster_cds.R "$@"
