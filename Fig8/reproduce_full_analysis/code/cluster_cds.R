library(tidyverse)
library(monocle3)
options(future.globals.maxSize=1*1000^3)
set.seed(1234)

args = commandArgs(trailingOnly=TRUE)
dataset = as.character(args[[1]])
split = as.character(args[[2]])
seed = as.character(args[[3]])

cds_loc <- paste0("/project2/gilad/jpopp/count-splitting-nb/data/fca_", dataset, "_", split, "_cds_", seed, ".rds")

if (dataset == "kidney") {
  resolution <- "1e-6"
} else if (dataset == "metanephric") {
  resolution <- "2e-5"
}

cluster_assignment_loc <- paste0("/project2/gilad/jpopp/count-splitting-nb/data/fca_", dataset, "_", split, "_clusters_", resolution, "_", seed, ".tsv")
cds_clustered_loc <- paste0("/project2/gilad/jpopp/count-splitting-nb/data/fca_", dataset, "_", split, "_clusteredcds_", resolution, "_", seed, ".rds")

cds <- readRDS(cds_loc)

cds <- cluster_cells(cds, resolution=as.numeric(resolution))
cluster_assignments <- as_tibble(clusters(cds), rownames="cell") %>%
  write_tsv(cluster_assignment_loc)

saveRDS(cds, cds_clustered_loc)
