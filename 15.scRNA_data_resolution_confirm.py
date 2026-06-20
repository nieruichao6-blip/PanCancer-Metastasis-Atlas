import os
import gc
import pandas as pd
import numpy as np
import scanpy as sc
import diopy
import bbknn
import scrublet as scr
import matplotlib.pyplot as plt
import scanpy.external as sce
import argparse

sc.settings.set_figure_params(dpi=500,figsize=(5, 5))
sc.logging.print_header()

# data dir
input_dir = " "
output_dir = " "
os.listdir(input_dir)

## pre work
cell_type_list = ["Myeloid_cell"]
# cell_type_list = ["B_cell","Plasma_cell","Myeloid_cell","Epithelia","T&NK_cell"]
# cell_type_list = ["Neutrophils","Endothelia","Fibroblast","B_cell","Plasma_cell","Myeloid_cell","Epithelia","T&NK_cell"]
cluster_method = "leiden"
bbknn_ridge = True
batch_method = "bbknn"
n_iterations = -1
resolution_list = [0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0]
random_state = 123

for cell_type in cell_type_list:

    ## data read
    scRNA_current = sc.read_h5ad(f"{input_dir}/{cell_type}/scRNA_bbknn.h5ad")

    ## dir create
    output_file = f"{output_dir}/{cell_type}"
    os.makedirs(output_file, exist_ok=True)

    ## umap
    sc.tl.umap(scRNA_current)

    ## diff res
    for resolution in resolution_list:
        print(f"############################# cell type : {cell_type} - resolution : {resolution} process ##########################################")

        ### cluster
        if batch_method == "harmony":
            sc.pp.neighbors(scRNA_current, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=40)
        elif batch_method == "scvi":
            sc.pp.neighbors(scRNA_current, use_rep="X_scVI")

        if cluster_method == "leiden":
            sc.tl.leiden(scRNA_current, resolution=resolution, random_state=random_state, key_added=f"leiden_{resolution}",n_iterations=n_iterations)
            gc.collect()
        else:
            sc.tl.louvain(scRNA_current, resolution=resolution, random_state=random_state)

        ### data plot
        sc.pl.umap(scRNA_current, color=f"leiden_{resolution}",show = False, legend_loc='on data')
        plt.savefig(f'{output_file}/scRNA_umap_cluster_resolution_{resolution}.png', dpi=500)
        plt.close()

        gc.collect

    ## sample plot
    sc.pl.umap(scRNA_current, color='study_ID',show = False)
    plt.savefig(f'{output_file}/scRNA_umap_sample.png', dpi=500)
    sc.pl.umap(scRNA_current, color='tissue_type',show = False)
    plt.savefig(f'{output_file}/scRNA_umap_tissue.png', dpi=500)

    ## data save
    diopy.output.write_h5(scRNA_current, file = f"{output_file}/scRNA_cell_type.h5",save_X=False)

    ## space free
    del scRNA_current
    gc.collect








