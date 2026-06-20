# environment -------------------------------------------------------------
rm(list = ls())
options(java.parameters = "-Xmx16384m")
options(stringsAsFactors = F)
options(warn=-1)
options(pkgType="source")
# options(future.globals.maxSize = 30000)
Sys.setenv(VROOM_CONNECTION_SIZE=500072)
path.sep <- .Platform$file.sep
# Sys.setenv(JAGS_HOME=" ")

library("idmap3")
library("GEOquery")
library("tools")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("optparse")
library("reshape2")
library("ggfortify")
library("ggpubr")
library("scales")
library("rJava")
library("xlsxjars")
library('glue')
library('funr')
library('xlsx')
library('limma')
library('edgeR')
library("ggrepel")
library("clusterProfiler")
library("ggridges")
library("statmod")
library("S4Vectors")
library("stats4")
library("BiocGenerics")
library("AnnotationDbi")
library("FactoMineR")
library("factoextra")
library("TCGAbiolinks")
library("SummarizedExperiment")
library("devtools")
library("Seurat")
library("monocle")
library("CCA")
library("clustree")
library("SCpubr")
library("harmony")
library("randomcoloR")
library("patchwork")
library("harmony")
library("ggprism")
library("scRNAtoolVis")
library("rjags")
library("AUCell")
library("GSVA")
library("GSEABase")
library("plyr")
library("randomcoloR")
library("future")
library("ggforce")
library("ggsci")
library("doFuture")
library("infercnv")
library("png")
library("dior")
library("sceasy")
library("reticulate")
library("Cottrazm")
(.packages())
setwd(' ')
display.brewer.all()

function_plot <- function(filename_prefix, width, height){
  ggsave(filename = str_c(filename_prefix, '.pdf'), plot = last_plot(), device = 'pdf', width = width, height = height, units = 'mm', dpi = 500)
  ggsave(filename = str_c(filename_prefix, '.png'), plot = last_plot(), device = 'png', width = width, height = height, units = 'mm', dpi = 2000)
  ggsave(filename = str_c(filename_prefix, '.svg'), plot = last_plot(), device = 'svg', width = width, height = height, units = 'mm', dpi = 500)
  try(ggsave(filename = str_c(filename_prefix, '.wmf'), plot = last_plot(), device = 'wmf', width = width, height = height, units = 'mm', dpi = 500),
      silent = T)
}


# data dir ----------------------------------------------------------------
cancer_dir <- " "
output_dir <- " "

# read seurat object ------------------------------------------------------
cancer_file <- list.files(cancer_dir)
# unlink(output_dir, recursive = T, force = T, expand = T)
# dir.create(output_dir, recursive = T)

## run pre
for(file in 1:length(cancer_file)){
  file <- 4
  print(str_glue("################################ {cancer_file[file]} Run Cluster ###################################"))
  set.seed(1234)
  scRNA_list <- list()
  input_file <- str_glue("{cancer_dir}/{cancer_file[file]}")
  file_name <- list.files(input_file)
  
  ## dir create
  output_file <- str_glue("{output_dir}/{cancer_file[file]}")
  unlink(output_file, recursive = T, force = T, expand = T)
  dir.create(output_file, recursive = T)
  
  for(i in 6:length(file_name)){
    # i<- 4
    print(str_glue("           -------------- {file_name[i]} Sample Run"))
    ## data read
    stRNA_list <- readRDS(str_glue("{cancer_dir}/{cancer_file[file]}/{file_name[i]}"))
    
    # stRNA_list <- readRDS(" ")
    ## SCT nor
    stRNA_list <- SCTransform(stRNA_list, assay = "spatial", verbose = FALSE, variable.features.n = 5000)
    
    ## cluster
    stRNA_list <- RunPCA(stRNA_list, assay = "SCT", verbose = FALSE)
    stRNA_list <- FindNeighbors(stRNA_list, reduction = "pca", dims = 1:30)
    stRNA_list <- FindClusters(stRNA_list, verbose = FALSE ,resolution = 1)
    stRNA_list <- RunUMAP(stRNA_list, reduction = "pca", dims = 1:30)
    stRNA_list <- RunTSNE(stRNA_list, reduction = "pca", dims = 1:30)
    
    # saveRDS(stRNA_list," ")
    ## cluster plot
    # SpatialDimPlot(stRNA_list, label = TRUE, label.size = 3,group.by = 'seurat_clusters')
    
    
    ## data save
    saveRDS(stRNA_list,str_glue("{output_dir}/{cancer_file[file]}/{file_name[i]}"))
    
  }
}








