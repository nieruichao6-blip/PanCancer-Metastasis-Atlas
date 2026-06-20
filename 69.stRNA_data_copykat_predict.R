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
library("ape")
library("magrittr")
library("DropletUtils")
library("copykat")
(.packages())
# setwd(' ')
display.brewer.all()

function_plot <- function(filename_prefix, width, height){
  ggsave(filename = str_c(filename_prefix, '.pdf'), plot = last_plot(), device = 'pdf', width = width, height = height, units = 'mm', dpi = 500)
  ggsave(filename = str_c(filename_prefix, '.png'), plot = last_plot(), device = 'png', width = width, height = height, units = 'mm', dpi = 1500)
  # ggsave(filename = str_c(filename_prefix, '.svg'), plot = last_plot(), device = 'svg', width = width, height = height, units = 'mm', dpi = 500)
}

# data color --------------------------------------------------------------
spot_color <- c("malignant"="#E41A1C","normal"="#48dc88")

# data dir ----------------------------------------------------------------
input_dir <- " "
output_dir <- " "
output_plot <- " "

# data process ------------------------------------------------------------
dir_name <- list.files(input_dir)
for(dir in 1:length(dir_name)){
  dir <- 6
  print(str_glue("############################## process {dir_name[dir]} copykat ############################"))
  
  
  input_file <- str_glue("{input_dir}/{dir_name[dir]}")
  output_file <- str_glue("{output_dir}/{dir_name[dir]}")
  output_plot_file <- str_glue("{output_plot}/{dir_name[dir]}")
  file_name <- list.files(input_file)
  
  ## output file create
  unlink(output_file, recursive = T, force = T, expand = T)
  dir.create(output_file, recursive = T)
  unlink(output_plot_file, recursive = T, force = T, expand = T)
  dir.create(output_plot_file, recursive = T)
  
  for(i in 7:length(file_name)){
    i <- 12
    
    ## data prepare
    current_stRNA <- readRDS(str_glue("{input_file}/{file_name[i]}"))
    exp.rawdata <- as.matrix(current_stRNA@assays$spatial@counts)

    ## copykat process
    copykat_process <- copykat(rawmat = exp.rawdata,id.type = "S",ngene.chr = 5,win.size = 25,KS.cut = 0.1,
                               distance = "euclidean",output.seg = "FALSE",plot.genes = "TRUE",
                               genome = "hg20",n.cores = 40)

    ## result transform
    copykat_predict <- as.data.frame(copykat_process$prediction)
    table(copykat_predict$copykat.pred)

    ## add to object
    colnames(copykat_predict)[1] <- "cell_name"
    meta_data <- current_stRNA@meta.data %>%
      rownames_to_column(var = "cell_name") %>%
      left_join(.,copykat_predict,by = "cell_name") %>%
      mutate(copykat_result = ifelse(copykat.pred == "aneuploid","malignant","normal")) %>%
      column_to_rownames(var = "cell_name")
    current_stRNA@meta.data <- meta_data
    current_stRNA <- subset(current_stRNA,copykat_result %in% c("malignant","normal"))
    
    current_stRNA <- readRDS(str_glue("{output_file}/{file_name[i]}"))
    
    ## data plot
    # SpatialFeaturePlot(current_stRNA,features = "KRT19") + NoLegend()
    SpatialDimPlot(current_stRNA,group.by = "copykat_result",cols = spot_color,stroke = NA,pt.size.factor = 1.5) + NoLegend()
    function_plot(filename_prefix = str_glue("{output_plot_file}/{file_name[i]}_plot"), width = 200, height = 200)
    
    ## data save
    saveRDS(current_stRNA,str_glue("{output_file}/{file_name[i]}"))
    
    
    
  }
  
}

