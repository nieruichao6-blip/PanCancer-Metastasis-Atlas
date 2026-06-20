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
library("CellChat")
library("infercnv")
library("copykat")
library("dior")
library("sceasy")
library("msigdbr")
library("fgsea")
library("GeneNMF")
library("UCell")
library("SCP")
(.packages())
setwd(' ')

function_plot <- function(filename_prefix, width, height){
  ggsave(filename = str_c(filename_prefix, '.pdf'), plot = last_plot(), device = 'pdf', width = width, height = height, units = 'mm', dpi = 500)
  ggsave(filename = str_c(filename_prefix, '.svg'), plot = last_plot(), device = 'svg', width = width, height = height, units = 'mm', dpi = 500)
  ggsave(filename = str_c(filename_prefix, '.png'), plot = last_plot(), device = 'png', width = width, height = height, units = 'mm', dpi = 500)
}

nrc_theme <- theme(panel.background = element_rect(fill = "transparent")) + 
  theme(panel.border = element_rect(color = "transparent", fill='transparent')) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(axis.title = element_text(size = 15, face = 'bold')) +
  theme(axis.line = element_line(color = 'black', linewidth = 1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1 )) +
  theme(axis.ticks = element_line(linewidth = 1), axis.ticks.length.y = unit("2", "mm")) +
  theme(axis.text = element_text(color="black", size = 12, face = 'bold')) +
  theme(legend.text = element_text(color="black", size = 12, face = 'bold'),
        legend.title = element_text(size = 15, face = 'bold')) 

cluster_theme <- theme(plot.title = element_text(face = "bold",size = 12,color="black",hjust = 0.5),
                       axis.title = element_text(face = "bold",size = 12,color ="black"), 
                       axis.text = element_text(face = "bold",size= 12,color = "black"),
                       panel.grid.minor.y = element_blank(),
                       panel.grid.minor.x = element_blank(),
                       axis.text.x = element_text(angle = 0, hjust = 0.5 ),
                       panel.grid=element_blank(),
                       legend.text = element_text(face = "bold",size= 12),
                       legend.title= element_text(face = "bold",size= 12))


# data dir ----------------------------------------------------------------
stRNA_dir <- " "
scRNA_dir <- " "
output_dir <- " "


# data read ---------------------------------------------------------------
module_gene_list <- readRDS(str_glue("{scRNA_dir}/meta_program_top50.RData"))


# data process ------------------------------------------------------------
st_file <- list.files(stRNA_dir)

## AUcell run
for(file in 1:length(st_file)){
  # file <- 2
  print(str_glue("################################ {st_file[file]} Run AUcell ###################################"))
  set.seed(1234)
  scRNA_list <- list()
  input_file <- str_glue("{stRNA_dir}/{st_file[file]}")
  file_name <- list.files(input_file)
  
  ## dir create
  output_file <- str_glue("{output_dir}/{st_file[file]}")
  unlink(output_file, recursive = T, force = T, expand = T)
  dir.create(output_file, recursive = T)
  
  for(i in 1:length(file_name)){
    # i<- 1
    print(str_glue("           -------------- {file_name[i]} Sample Run"))
    ## data read
    stRNA_list <- readRDS(str_glue("{stRNA_dir}/{st_file[file]}/{file_name[i]}"))
    
    ## AUcell
    cells_rankings <- AUCell_buildRankings(stRNA_list@assays$SCT@data,splitByBlocks=TRUE)
    cells_AUC <- AUCell_calcAUC(module_gene_list, cells_rankings,
                                aucMaxRank=nrow(cells_rankings)*0.1)
    auc_score <- as.data.frame(t(getAUC(cells_AUC)[names(module_gene_list)[names(module_gene_list) %in% rownames(cells_AUC)],]))
    stRNA_list <- AddMetaData(stRNA_list,auc_score)
    
    ## data save
    saveRDS(stRNA_list,str_glue("{output_file}/{file_name[i]}"))
    
  }
}


