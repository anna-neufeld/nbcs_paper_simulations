## Real Data Analysis
### Overview
The code to run the full analysis is contained in the `code` folder. The `slurm_wrappers` folder 
contains scripts to implement the count splitting pipeline efficiently under different distributional assumptions (Poisson vs NB), across multiple subsets of the data (under different random seeds used for count splitting), and for different datasets (full fetal kidney vs metanephric cell type)

### Reproducing the analysis
*Note:* in order to reproduce the analysis, restructure the paths at the top and bottom of each file to match your own directory structure.

1. Data Acquisition
- `code/fca_kidney_download.sh` allows you to download the fetal cell atlas data necessary for the analysis

2. Filter Count Matrices
- `code/subset_counts.R` filters each count matrix (all kidney cells, and metanephric cells only) to genes expressed in at least 10 cells 

3. Apply Count Splitting
- `code/countsplit.R` applies count splitting (Poisson or NB, depending on the user input), with help 
from some functions contained in `code/helpers_countsplit.R`

4. Create CDS
- `code/create_cds.R` creates a cell_data_set object from a counts matrix (or a train/test split after count splitting) for subsequent use in a Monocle pipeline. It also applies normalization and generates multiple embeddings of the data (PCA and UMAP) as recommended in the Monocle3 tutorial. 

5. Apply clustering
- `code/cluster_cds.R` clusters a cell_data_set object as in the Monocle3 tutorial, using a resolution parameter specific to kidney or metanephric cells which partitions cells into a similar number of cell types as used by Cao et al. 2020

6. Apply intra-dataset cross-validation
- `code/intradataset_cv_*.R` applies intra-dataset cross-validation to evaluate clusters in kidney or metanephric cells

7. Compute adjusted rand index
- `code/traintest_comparison_ari.R` takes two sets of clustering assignments for a dataset and computes the adjusted rand index between them
