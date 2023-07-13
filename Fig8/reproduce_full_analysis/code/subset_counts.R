library(tidyverse)
library(monocle3)
options(future.globals.maxSize=1*1000^3)
set.seed(1234)

kidney_loc <- "/project2/gilad/jpopp/count-splitting-nb/data/fca_kidney_raw_counts.rds"
cell_metadata_loc <- "/project2/gilad/jpopp/count-splitting-nb/data/fca_cell_metadata.rds"
gene_metadata_loc <- "/project2/gilad/jpopp/count-splitting-nb/data/fca_gene_metadata.rds"

kidney_all <- readRDS(kidney_loc)
cell_metadata <- readRDS(cell_metadata_loc) %>%
  filter(sample %in% colnames(kidney_all)) %>%
  column_to_rownames("sample")
gene_metadata <- readRDS(gene_metadata_loc) %>%
  filter(gene_id %in% rownames(kidney_all))

# filter full dataset to genes expressed in at least 10 cells
gene_keepers_kidney <- rownames(kidney_all)[rowSums(kidney_all > 0) >= 10]
kidney <- kidney_all[gene_keepers_kidney,]

# filter to metanephric cells
cell_keepers_metanephric <- cell_metadata %>%
  filter(Organ_cell_lineage=="Kidney-Metanephric cells") %>%
  rownames
metanephric <- kidney_all[,cell_keepers_metanephric]

# filter metanephric subset to genes expressed in at least 10 metanephric cells
gene_keepers_metanephric <- rownames(metanephric)[rowSums(metanephric > 0) >= 10]
metanephric <- metanephric[gene_keepers_metanephric,]

# save both objects
saveRDS(kidney, "/project2/gilad/jpopp/count-splitting-nb/data/fca_kidney_full_counts.rds")
saveRDS(metanephric, "/project2/gilad/jpopp/count-splitting-nb/data/fca_metanephric_full_counts.rds")