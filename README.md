#nbcs_paper

This repository contains all code needed to reproruce the analyses from "Negative binomial count splitting for single-cell RNA sequencing data", by Neufeld, Popp, Gao, Battle, and Witten (2023). 

Each figure from the paper has a corresponding folder, with its own README and R files needed to create that figure. 

The file "simulation_functions.R" is important for figures 3-7, and will need to be sourced before running the code in any of those folders. 

The res folder contains output from our own simulation studies, which can be used to more quickly recreate figures 3,4,5 and 6. The res folder only contains results from 500 replicates of our own simulation studies, due to github repository size constraints. The remaining 1500 replicates can be reproduced by running the simulation code from scratch. 

