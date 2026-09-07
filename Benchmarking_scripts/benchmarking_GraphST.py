# IMPORTANT: requires python 3.8
import os
import sys
import glob
import random
import torch
import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd
from GraphST import GraphST
from GraphST.utils import clustering 
from memory_profiler import memory_usage
import time

torch.set_num_threads(1)

# function for running GraphST
def run_graphst(in_file_path, seed, n_clusters, start, end, increment):
    # setting seed
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    
    # setting the device to cpu
    device = torch.device('cpu')
    
    adata_list = []
    # reading the anndata file names
    filenames = [os.path.basename(f) for f in glob.glob(in_file_path + "*.h5ad")]
    for in_file in filenames:
        # reading anndata file
        adata = sc.read_h5ad(in_file_path + in_file)
        # running graphST
        model = GraphST.GraphST(adata, datatype = 'Stereo', device = device)
        # running model training
        adata = model.train()
        # adding output to list
        adata_list.append(adata)
    # joining all sample adata
    adata_comb = ad.concat(adata_list)
    # louvain clustering
    clustering(adata_comb, n_clusters, method='louvain', start = start, end = end, increment = increment)
    return adata_comb




#### SeqFISH mouse embryo data
seed = 567
n_clusters = 23
start, end, increment = 0.1, 3.0, 0.01 # default
in_file_path = "./data_outputs/annDataFiles/mouse_embryo/"

start_time = time.perf_counter()
mem_usage, adata_out = memory_usage((run_graphst, (in_file_path, seed, n_clusters, start, end, increment)), retval = True, max_usage = True)
end_time = time.perf_counter()
end_time - start_time

# saving data
labels_df = adata_out.obs[['domain', 'louvain']].copy()
labels_df.to_csv("./data_outputs/benchmarkOut/graphst_me.csv", index = True, index_label = "cell_id")

metrics_df = pd.DataFrame({
    "peak_memory_MB": [mem_usage],    # Update the unit label if yours is GB/KB
    "runtime_sec": [end_time - start_time]  # Update the unit label if yours is minutes/hours
})
metrics_df.to_csv("./data_outputs/benchmarkOut/mem_me_graphst.csv", index=False)

print(f"Maximum memory usage: {mem_usage} MiB")
print(f"Total running time: {end_time - start_time} seconds")




#### Xenium breast cancer data
seed = 567
n_clusters = 19
start, end, increment = 0.1, 3.0, 0.01 # default
in_file_path = "./data_outputs/annDataFiles/breast_cancer/"

start_time = time.perf_counter()
mem_usage, adata_out = memory_usage((run_graphst, (in_file_path, seed, n_clusters, start, end, increment)), retval = True, max_usage = True)
end_time = time.perf_counter()
end_time - start_time
# Run ended with memory usage >400 GB




#### MERFISH mouse hypothalamus
seed = 567
n_clusters = 16
start, end, increment = 0.1, 3.0, 0.01 # default
in_file_path = "./data_outputs/annDataFiles/mouse_brain/"

start_time = time.perf_counter()
mem_usage, adata_out = memory_usage((run_graphst, (in_file_path, seed, n_clusters, start, end, increment)), retval = True, max_usage = True)
end_time = time.perf_counter()
end_time - start_time
# Runtime exceeded 4 days.




#### CosMx lung cancer data
seed = 567
n_clusters = 22
start, end, increment = 0.1, 3.0, 0.01 # default
in_file_path = "./data_outputs/annDataFiles/lung_cancer/"

start_time = time.perf_counter()
mem_usage, adata_out = memory_usage((run_graphst, (in_file_path, seed, n_clusters, start, end, increment)), retval = True, max_usage = True)
end_time = time.perf_counter()
end_time - start_time
# Run ended with memory usage >400 GB




#### Stereo-seq human ovarian cancer
seed = 567
n_clusters = 13
start, end, increment = 0.1, 3.0, 0.01 # default
in_file_path = "./data_outputs/annDataFiles/ovarian_cancer/"

start_time = time.perf_counter()
mem_usage, adata_out = memory_usage((run_graphst, (in_file_path, seed, n_clusters, start, end, increment)), retval = True, max_usage = True)
end_time = time.perf_counter()
end_time - start_time
# Run ended with memory usage >400 GB




#### VisiumHD human COAD
seed = 567
n_clusters = 16
start, end, increment = 0.1, 3.0, 0.01 # default
in_file_path = "./data_outputs/annDataFiles/colon_cancer/"

start_time = time.perf_counter()
mem_usage, adata_out = memory_usage((run_graphst, (in_file_path, seed, n_clusters, start, end, increment)), retval = True, max_usage = True)
end_time = time.perf_counter()
end_time - start_time
# Run ended with memory usage >400 GB




#### SeqFISH mouse embryo data with PASTE alignment
# setting seed
random.seed(567)
np.random.seed(567)
torch.manual_seed(567)
# setting the device to cpu
device = torch.device('cpu')
# reading merged anndata file with aligned coordinates
adata = sc.read_h5ad("./data_outputs/annDataFiles/mouse_embryo/Embryo_aligned_merged.h5ad")
# running graphST
model = GraphST.GraphST(adata, datatype = 'Stereo', device = device)
# running model training
adata = model.train()
# louvain clustering
clustering(adata, n_clusters = 23, method = 'louvain', start = 0.1, end = 3.0, increment = 0.01)
# saving data
labels_df = adata.obs[['domain', 'louvain']].copy()
labels_df.to_csv("./data_outputs/benchmarkOut/paste_graphst_me.csv", index = True, index_label = "cell_id")
adata.write_h5ad("./data_outputs/annDataFiles/mouse_embryo/Embryo_aligned_merged_graphst.h5ad")
