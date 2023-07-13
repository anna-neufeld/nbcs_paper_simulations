library(caret)
library(corrplot)
library(mclust)
library(scmap)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
dataset = as.character(args[[1]])
split = as.character(args[[2]])

if (dataset == "metanephric") {
  cluster_res <- "2e-5"
} else {
  cluster_res <- "1e-6"
}

seeds <- c(75, 803, 498, 184, 463, 552, 145, 191, 401, 254)

ari_list <- c()

for (s in seeds) {
  clusts_train_loc = paste0("/project2/gilad/jpopp/count-splitting-nb/data/fca_", 
                            dataset, "_", split, "train_clusters_", cluster_res, "_", s, ".tsv")
  clusts_test_loc = paste0("/project2/gilad/jpopp/count-splitting-nb/data/fca_", 
                           dataset, "_", split, "test_clusters_", cluster_res, "_", s, ".tsv")
  
  clusts_train <- read_tsv(clusts_train_loc) %>%
    rename(train_cluster=value)
  clusts_test <- read_tsv(clusts_test_loc) %>%
    rename(test_cluster=value)
  
  both_clusts <- inner_join(clusts_train, clusts_test, by="cell")
  
  ari_list <- c(ari_list, adjustedRandIndex(both_clusts$train_cluster, both_clusts$test_cluster))
}

save(ari_list, file=paste0("/project2/gilad/jpopp/count-splitting-nb/data/fca_", 
                           dataset, "_", split, "train_clusters_", cluster_res, "_ari_all.Rdata"))
