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
library("SCP")
(.packages())
setwd(' ')
display.brewer.all()


# data dir ----------------------------------------------------------------
input_list_dir <- " "
output_dir <- " "
data_format <- "10x_format"


# scRNA data load ---------------------------------------------------------

## run step
dir_name <- list.files(input_list_dir)
for(dir in 1:length(dir_name)){
  print(str_glue("############################## Create {dir_name[dir]} Seurat Object ############################"))
  # dir <- 1
  scRNA_list <- list() 
  cccc
  output_file <- str_glue("{output_dir}/{dir_name[dir]}")
  file_name <- list.files(input_file)
  
  ## output file create
  unlink(output_file, recursive = T, force = T, expand = T)
  dir.create(output_file, recursive = T)
  
  ## seurat object create
  for(i in 1:length(file_name)){
    # i <- 1
    print(str_glue("     ----------- {file_name[i]} read"))
    
    counts <- Read10X(data.dir = str_glue("{input_file}/{file_name[i]}/"))
    scRNA_list[[i]] <- CreateSeuratObject(counts,project = file_name[i],
                                          min.cells = 1, min.features = 200)
    scRNA_list[[i]]$orig.ident <- file_name[i]
  }
  
  ## seurat object output
  for(i in 1:length(scRNA_list)){
    # i <- 1
    print("already output !!!")
    output_data <- scRNA_list[[i]]
    saveRDS(output_data,str_glue("{output_file}/{str_sub(output_data@project.name,1,10)}.RData"))
  }
  
  ## env clean
  remove(scRNA_list)
  gc()
  
}




# # scRNA data load ---------------------------------------------------------
# ## attention data format
# input_list_dir <- " "
# data_format <- "exp_format"
# dir_name <- list.files(input_list_dir)
# scRNA_list <- list()
# 
# if(data_format == '10x_format'){
#   for(i in 1:length(dir_name)){
#     # i <- 19
#     counts <- Read10X(data.dir = str_glue("{input_list_dir}/{dir_name[i]}/"))
#     scRNA_list[[i]] <- CreateSeuratObject(counts,project = dir_name[i],
#                                           min.cells = 3, min.features = 200)
#   }
# }else if(data_format == 'exp_format'){
#   for(i in 1:length(dir_name)){
#     # i <- 1
#     exp_data <- read_tsv(gzfile(str_glue("{input_list_dir}/{dir_name[i]}")))
#     exp_data <- column_to_rownames(exp_data,var = 'gene')
#     scRNA_list[[i]] <- CreateSeuratObject(exp_data,project = str_sub(dir_name[i],1,10),
#                                           min.cells = 3, min.features = 200)
#     
#   }
# }else if(data_format == 'rds_format'){
#   scRNA_list <- read_rds(str_glue("{input_list_dir}/{dir_name}"))
# }else if(data_format == 'delim_format'){
#   for(i in 1:length(dir_name)){
#     # i <- 12
#     exp_data <- read_delim(gzfile(str_glue("{input_list_dir}/{dir_name[i]}")),delim = '\t')
#     # exp_data <- exp_data[,-1]
#     exp_data$GENE <- make.unique(exp_data$GENE)
#     exp_data <- column_to_rownames(exp_data,var = 'GENE')
#     scRNA_list[[i]] <- CreateSeuratObject(exp_data,project = dir_name[i],
#                                           min.cells = 3, min.features = 200)
#   }
# }else if(data_format == 'csv_format'){
#   for(i in 1:length(dir_name)){
#     i <- 1
#     exp_data <- read_csv(gzfile(str_glue("{input_list_dir}/{dir_name[i]}")))
#     colnames(exp_data)[1] <- "Name"
#     exp_data$Name <- make.unique(exp_data$Name)
#     exp_data <- column_to_rownames(exp_data,var = 'Name')
#     scRNA_list[[i]] <- CreateSeuratObject(exp_data,project = dir_name[i],
#                                           min.cells = 3, min.features = 200)
#     }
# }else if(data_format == 'h5_format'){
#   for(i in 1:length(dir_name)){
#     # i <- 1
#     exp_data <- Read10X_h5(str_glue("{input_list_dir}/{dir_name[i]}"))
#     scRNA_list[[i]] <- CreateSeuratObject(exp_data,project = dir_name[i],
#                                           min.cells = 3, min.features = 200)
#   }
# }
# 
# 
# 
# # data output rds format --------------------------------------------------
# for(i in 1:length(scRNA_list)){
#   # i <- 1
#   output_data <- scRNA_list[[i]]
#   saveRDS(output_data,str_glue("{output_dir}/{str_sub(output_data@project.name,1,10)}.RData"))
# }
# 
# 











