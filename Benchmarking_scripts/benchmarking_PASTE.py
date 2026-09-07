import scanpy as sc
import paste as pst
import anndata as ad
import pandas as pd
import numpy as np
import random


#### SeqFISH mouse embryo data with 3 samples

# setting seed for reproducibility
seed = 567
random.seed(seed)
np.random.seed(seed)
# reading the anndata files
me_adata1 = sc.read_h5ad("./data_outputs/annDataFiles/mouse_embryo/Sample_embryo1.h5ad")
me_adata2 = sc.read_h5ad("./data_outputs/annDataFiles/mouse_embryo/Sample_embryo2.h5ad")
me_adata3 = sc.read_h5ad("./data_outputs/annDataFiles/mouse_embryo/Sample_embryo3.h5ad")
# adding to list
slices_me = [me_adata1, me_adata2, me_adata3]
# aligning slices
pis_me = []
for i in range(len(slices_me) - 1):
    pi_me = pst.pairwise_align(slices_me[i], slices_me[i+1])
    pis_me.append(pi_me)
# generating aligned coordinates stored in adata.obsm['spatial']
aligned_slices_me = pst.stack_slices_pairwise(slices_me, pis_me)
# extracting aligned coordinates
aligned_list_me = []
for i, adata in enumerate(aligned_slices_me):
    df = pd.DataFrame(adata.obsm['spatial'], index = adata.obs_names, columns = ['x_aligned', 'y_aligned'])
    df['Sample'] = f"Sample{i+1}"
    aligned_list_me.append(df)
me_aligned_coords = pd.concat(aligned_list_me)
# saving common coordinates
me_aligned_coords.to_csv("./data_outputs/benchmarkOut/Aligned_coords_ME.csv")
# merging all sample data and saving
adata_merged_me = ad.concat(aligned_slices_me)
adata_merged_me.write_h5ad("./data_outputs/annDataFiles/mouse_embryo/Embryo_aligned_merged.h5ad")
