#!/bin/bash
#$ -cwd

Rscript run_roleEps_diffExp.R --simname $1 --nreps $2

