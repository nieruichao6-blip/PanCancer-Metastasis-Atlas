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
library("AnnoProbe")

(.packages())
setwd(' ')
display.brewer.all()

function_plot <- function(filename_prefix, width, height){
  ggsave(filename = str_c(filename_prefix, '.pdf'), plot = last_plot(), device = 'pdf', width = width, height = height, units = 'mm', dpi = 500)
  ggsave(filename = str_c(filename_prefix, '.png'), plot = last_plot(), device = 'png', width = width, height = height, units = 'mm', dpi = 1500)
  # ggsave(filename = str_c(filename_prefix, '.svg'), plot = last_plot(), device = 'svg', width = width, height = height, units = 'mm', dpi = 500)
}


# data dir ----------------------------------------------------------------
input_dir <- " "
output_dir <- " "

# data process ------------------------------------------------------------
dir_name <- list.files(input_dir)
for(dir in 1:length(dir_name)){
  dir <- 4
  print(str_glue("############################## process {dir_name[dir]} infercnv ############################"))
  
  input_file <- str_glue("{input_dir}/{dir_name[dir]}")
  output_file <- str_glue("{output_dir}/{dir_name[dir]}")
  file_name <- list.files(input_file)
  
  ## output file create
  unlink(output_file, recursive = T, force = T, expand = T)
  dir.create(output_file, recursive = T)
  
  for(i in 6:length(file_name)){
    # i <- 1
    
    output_sample <- str_glue("{output_dir}/{dir_name[dir]}/{file_name[i]}/")
    unlink(output_sample, recursive = T, force = T, expand = T)
    dir.create(output_sample, recursive = T)
    
    ## data read
    current_stRNA <- readRDS(str_glue("{input_file}/{file_name[i]}"))
    dat <- GetAssayData(current_stRNA,layer = 'counts',assay = 'spatial')
    groupinfo <- data.frame(v1 = colnames(dat),
                            v2 = current_stRNA@meta.data$copykat_result)
    geneInfor <- annoGene(rownames(dat),"SYMBOL","human")
    geneInfor$chr_num <- as.numeric(sub("chr", "", geneInfor$chr))
    geneInfor <- geneInfor[with(geneInfor,order(chr_num,start)),c(1,4:6)]
    geneInfor <- geneInfor[!duplicated(geneInfor[,1]),]
    
    dat <- dat[rownames(dat) %in% geneInfor[,1],]
    dat <- dat[match(geneInfor[,1],rownames(dat)),]
    
    ## data prepare
    expFile <- 'expFile.txt'
    colnames(dat) <- gsub("-", "_", colnames(dat))
    write.table(dat, file = str_glue("{output_sample}/{expFile}"), sep = '\t', quote = F)
    groupFiles <- 'groupFiles.txt'
    groupinfo$v1 <- gsub("-", "_", groupinfo$v1)
    write.table(groupinfo,file = str_glue("{output_sample}/{groupFiles}"), sep = '\t',
            quote = F, col.names = F, row.names = F)
    geneFile <- 'geneFile.txt'
    write.table(geneInfor, file = str_glue("{output_sample}/{geneFile}"), sep = '\t',
            quote = F, col.names = F, row.names = F)
            
    ## infercnv run
    infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = str_glue("{output_sample}/{expFile}"),
                                     annotations_file = str_glue("{output_sample}/{groupFiles}"),
                                     delim = "\t",
                                     gene_order_file = str_glue("{output_sample}/{geneFile}"),
                                     ref_group_names = "normal"
                                                         )

    infercnv_obj2 <- infercnv::run(infercnv_obj,
                               cutoff = 0.1,
                               out_dir = str_glue("{output_sample}/infercnv_run_result"),
                               cluster_by_groups = T,
                               hclust_method = "ward.D2",
                               analysis_mode = "subclusters",
                               denoise = TRUE,
                               HMM = F,
                               plot_steps = F,
                               leiden_resolution = "auto",
                               num_threads = 20,
                               # tumor_subcluster_partition_method = "random_trees",
                               output_format = "pdf"
                               )
    
    gc()
    
  }
  
  
  
}



