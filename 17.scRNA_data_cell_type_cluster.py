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

sc.settings.set_figure_params(dpi=1000,figsize=(5, 5))
sc.logging.print_header()

parser = argparse.ArgumentParser(
    description="multi tumor code cell type cluster",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

# parameter set
parser.add_argument("-i", "--input_dir", required=True, help="文件所在文件夹")
parser.add_argument("-ct", "--cell_type", required=True, help="分析的细胞类型")
parser.add_argument("-cm", "--cluster_method", required=True, help="聚类方式")
# parser.add_argument("-br", "--bbknn_ridge", required=True, default=True, help="聚类方式")
parser.add_argument("-r", "--resolution", type=float, default=2.0,  help="分辨率")
parser.add_argument("-ni", "--n_iterations", type=int, default=-1, help="迭代参数")
parser.add_argument("-rs", "--random_state", type=int, default=123, help="种子设置")
parser.add_argument("-bm", "--batch_method", required=True,  help="去批次的方式")
parser.add_argument("-o", "--output_dir", default="infercnv_results", help="输出目录")

args = parser.parse_args()

# data dir
# input_dir = " "
# output_dir = " "
# os.listdir(input_dir)

## pre work
# cell_type = "Epithelia_tumor"
# cluster_method = "leiden"
bbknn_ridge = True
# batch_method = "bbknn"
# resolution = 2
# random_state = 123
# n_iterations = -1

## data read
# scRNA_current = diopy.input.read_h5(file = f"{input_dir}/scRNA_{cell_type}.h5")
scRNA_current = sc.read_h5ad(f"{args.input_dir}/{args.cell_type}/scRNA_bbknn.h5ad")

## dir create
output_file = f"{args.output_dir}/{args.cell_type}"
os.makedirs(output_file, exist_ok=True)
gc.collect()

print(f"############################# cell type : {args.cell_type} - resolution : {args.resolution} process ##########################################")

## cluster
if args.batch_method == "harmony":
    sc.pp.neighbors(scRNA_current, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=40)
elif args.batch_method == "scvi":
    sc.pp.neighbors(scRNA_current, use_rep="X_scVI")

if args.cluster_method == "leiden":
    sc.tl.leiden(scRNA_current, resolution=args.resolution, random_state=args.random_state, key_added="leiden",n_iterations=args.n_iterations)
    gc.collect()
else:
    sc.tl.louvain(scRNA_current, resolution=args.resolution, random_state=args.random_state)
    
if bbknn_ridge == False:
    sc.tl.umap(scRNA_current)

print(f"############################# start bbknn ridge ##########################################")
scRNA_current.write_h5ad(f"{output_file}/scRNA_cell_type_cluster_noridge.h5ad", compression="gzip")
diopy.output.write_h5(scRNA_current, file = f"{output_file}/scRNA_cell_type_cluster_noridge.h5",save_X=False)

## bbknn ridge
if bbknn_ridge == True:
    bbknn.ridge_regression(scRNA_current, batch_key=["sample_ID","study_ID","tissue_type"], confounder_key = ["leiden"])
    sc.tl.pca(scRNA_current, svd_solver='arpack', use_highly_variable=True, n_comps=50)
    sc.external.pp.bbknn(scRNA_current, batch_key= "sample_ID", n_pcs = 50)
    scRNA_current.write_h5ad(f"{output_file}/scRNA_cell_type_cluster_bbknn_ridge.h5ad", compression="gzip")
    gc.collect()
    sc.tl.leiden(scRNA_current, resolution=args.resolution, random_state=args.random_state,n_iterations=args.n_iterations)
    sc.tl.umap(scRNA_current)

sc.pl.umap(scRNA_current, color='leiden',show = False, legend_loc='on data')
plt.savefig(f'{output_file}/scRNA_umap_cluster.png', dpi=1000)

sc.pl.umap(scRNA_current, color='study_ID',show = False)
plt.savefig(f'{output_file}/scRNA_umap_sample.png', dpi=1000)

sc.pl.umap(scRNA_current, color='tissue_type',show = False)
plt.savefig(f'{output_file}/scRNA_umap_tissue.png', dpi=1000)

# data save
if args.cell_type == "Epithelia_tumor" :
    scRNA_current.write_h5ad(f"{output_file}/scRNA_cell_type_cluster.h5ad", compression="gzip")
    diopy.output.write_h5(scRNA_current, file = f"{output_file}/scRNA_cell_type_cluster.h5",save_X=True)
else: 
    diopy.output.write_h5(scRNA_current, file = f"{output_file}/scRNA_cell_type_cluster.h5",save_X=False)
