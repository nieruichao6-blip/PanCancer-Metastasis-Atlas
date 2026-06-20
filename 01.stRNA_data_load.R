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
(.packages())
setwd(' ')
display.brewer.all()


# data dir ----------------------------------------------------------------
input_list_dir <- " "
output_dir <- " "


# stRNA data load ---------------------------------------------------------

## run step
dir_name <- list.files(input_list_dir)
for(dir in 1:length(dir_name)){
  dir <- 4
  print(str_glue("############################## Create {dir_name[dir]} Seurat Object ############################"))

  input_file <- str_glue("{input_list_dir}/{dir_name[dir]}")
  output_file <- str_glue("{output_dir}/{dir_name[dir]}")
  file_name <- list.files(input_file)
  
  ## output file create
  unlink(output_file, recursive = T, force = T, expand = T)
  dir.create(output_file, recursive = T)
  
  ## seurat object create
  for(i in 6:length(file_name)){
    # i <- 2
    print(str_glue("          ------------ {file_name[i]} process "))
    if(file.exists(str_glue("{input_file}/{file_name[i]}/spatial/tissue_lowres_image.png")) & file.exists(str_glue("{input_file}/{file_name[i]}/filtered_feature_bc_matrix.h5"))){
      stRNA_object <- Seurat::Load10X_Spatial(data.dir = str_glue("{input_file}/{file_name[i]}/"),
                                                   filename = "filtered_feature_bc_matrix.h5",
                                                   assay='spatial',
                                                   slice = file_name[i])
      stRNA_object@meta.data$orig.ident <- file_name[i]
      stRNA_object@project.name <- file_name[i]
      
    }else if(file.exists(str_glue("{input_file}/{file_name[i]}/filtered_feature_bc_matrix/barcodes.tsv.gz"))){
      stRNA_list <- Read10X(str_glue("{input_file}/{file_name[i]}/filtered_feature_bc_matrix/"))
      image2 <- Read10X_Image(image.dir = file.path(str_glue("{input_file}/{file_name[i]}/"),"spatial"), filter.matrix = TRUE)
      stRNA_list <- CreateSeuratObject(counts = stRNA_list, assay = "spatial")
      image2 <- image2[Cells(x = stRNA_list)]
      DefaultAssay(stRNA_list = image2) <- "spatial"
      stRNA_list[["slice1"]] <- image2
      stRNA_object <- stRNA_list
      
      meta_data <- stRNA_object@meta.data
      colnames(meta_data)[2:3] <- c("nCount_spatial","nFeature_spatial")
      stRNA_object@meta.data  <- meta_data
      stRNA_object@meta.data$orig.ident <- file_name[i]
      stRNA_object@project.name <- file_name[i]
    }else{
      # img <- read_image(str_glue("{input_file}/{file_name[i]}/spatial/"))
      img <- Read10X_Image(str_glue("{input_file}/{file_name[i]}/spatial/"),image.name = "tissue_hires_image.png")
      stRNA_list <- Seurat::Load10X_Spatial(data.dir = str_glue("{input_file}/{file_name[i]}/"),
                                            assay = "spatial",
                                            filename = "filtered_feature_bc_matrix.h5",
                                            image = img)
      stRNA_list@images$slice1@scale.factors$lowres <- stRNA_list@images$slice1@scale.factors$hires
      stRNA_object <- stRNA_list
      
      meta_data <- stRNA_object@meta.data
      colnames(meta_data)[2:3] <- c("nCount_spatial","nFeature_spatial")
      stRNA_object@meta.data  <- meta_data
      stRNA_object@meta.data$orig.ident <- file_name[i]
      stRNA_object@project.name <- file_name[i]
      
      
    }
    
    ### output data
    saveRDS(stRNA_object,str_glue("{output_file}/{file_name[i]}.RData"))
    # dior::write_h5(stRNA_object, file=str_glue("{output_file}/{file_name[i]}.h5"),assay.name = "spatial",object.type = 'seurat')
    # sceasy::convertFormat(stRNA_object, from="seurat", to="anndata",assay = "spatial",outFile = str_glue("{output_file}/{file_name[i]}.h5ad"),main_layer = "counts")
    
  }
}



