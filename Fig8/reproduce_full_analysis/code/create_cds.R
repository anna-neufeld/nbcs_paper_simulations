library(tidyverse)
library(monocle3)
options(future.globals.maxSize=1*1000^3)
set.seed(1234)

args = commandArgs(trailingOnly=TRUE)
dataset = as.character(args[[1]])
split = as.character(args[[2]])
seed = as.character(args[[3]])

counts_loc <- paste0("/project2/gilad/jpopp/count-splitting-nb/data/fca_", dataset, "_", split, "_counts_", seed, ".rds")
cell_metadata_loc <- "/project2/gilad/jpopp/count-splitting-nb/data/fca_cell_metadata.rds"
gene_metadata_loc <- "/project2/gilad/jpopp/count-splitting-nb/data/fca_gene_metadata.rds"

pc_embedding_loc <- paste0("/project2/gilad/jpopp/count-splitting-nb/data/fca_", dataset, "_", split, "_pcs_", seed, ".tsv")
umap_embedding_loc <- paste0("/project2/gilad/jpopp/count-splitting-nb/data/fca_", dataset, "_", split, "_umap_", seed, ".tsv")
cds_loc <- paste0("/project2/gilad/jpopp/count-splitting-nb/data/fca_", dataset, "_", split, "_cds_", seed, ".rds")

counts <- readRDS(counts_loc)
cell_metadata <- readRDS(cell_metadata_loc) %>%
  dplyr::filter(sample %in% colnames(counts)) %>%
  column_to_rownames("sample")
gene_metadata <- readRDS(gene_metadata_loc) %>%
  dplyr::filter(gene_id %in% rownames(counts)) 

cds <- new_cell_data_set(counts,
                         cell_metadata = cell_metadata[colnames(counts),], #reorder if needed
                         gene_metadata = gene_metadata[rownames(counts),]) #reorder if needed

cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds)

pc_embedding <- as_tibble(reducedDim(cds), rownames="cell") %>% write_tsv(pc_embedding_loc)
umap_embedding <- as_tibble(reducedDim(cds, "UMAP"), rownames="cell") %>% write_tsv(umap_embedding_loc)

saveRDS(cds, cds_loc)