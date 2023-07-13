library(tidyverse)
library(sctransform)
library(Matrix)
library(Matrix.utils)
library(stats)
options(future.globals.maxSize=1*1000^3)
source("/project2/gilad/jpopp/count-splitting-nb/code/helpers_countsplit.R")

args = commandArgs(trailingOnly=TRUE)
dataset = as.character(args[[1]])
splitter = as.character(args[[2]])
seed = as.integer(args[[3]])
set.seed(seed)

counts_loc = paste0("../data/fca_", dataset, "_full_counts.rds")
counts = readRDS(counts_loc)

if (splitter %in% c("poisson", "Poisson", "pois", "Pois")) {
  counts_split <- countsplit(counts)
  
  counts_train <- counts_split$train
  saveRDS(counts_train, paste0("../data/fca_", dataset, "_poistrain_counts_", seed, ".rds"))
  
  counts_test <- counts_split$test
  saveRDS(counts_test, paste0("../data/fca_", dataset, "_poistest_counts_", seed, ".rds"))
} else if (splitter %in% c("nb", "NB", "negativebinomial", "NegativeBinomial", "negbin", "NegBin")) {
  sctransform_fit <- sctransform::vst(counts, min_cells=10, residual_type='none')
  overdisp_estimates <- sctransform_fit$model_pars_fit[,"theta"]
  saveRDS(sctransform_fit, paste0("../data/fca_", dataset, "_sctransform_fit_", seed, ".rds"))
  
  counts_split <- nb.countsplit(counts, overdisps=overdisp_estimates)
  
  counts_train <- counts_split$train
  saveRDS(counts_train, paste0("../data/fca_", dataset, "_nbtrain_counts_", seed, ".rds"))
  
  counts_test <- counts_split$test
  saveRDS(counts_test, paste0("../data/fca_", dataset, "_nbtest_counts_", seed, ".rds"))
} else {
  print(paste0("Error: \"", splitter, "\" is not a valid splitter for this data"))
}
