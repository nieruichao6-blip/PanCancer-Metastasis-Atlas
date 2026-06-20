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
(.packages())
# setwd(' ')
display.brewer.all()

function_plot <- function(filename_prefix, width, height){
  ggsave(filename = str_c(filename_prefix, '.pdf'), plot = last_plot(), device = 'pdf', width = width, height = height, units = 'mm', dpi = 500)
  ggsave(filename = str_c(filename_prefix, '.svg'), plot = last_plot(), device = 'svg', width = width, height = height, units = 'mm', dpi = 500)
  ggsave(filename = str_c(filename_prefix, '.png'), plot = last_plot(), device = 'png', width = width, height = height, units = 'mm', dpi = 2000)
}

nrc_theme <- theme(panel.background = element_rect(fill = "transparent")) + 
  theme(panel.border = element_rect(color = "transparent", fill='transparent')) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(axis.title = element_text(size = 15, face = 'bold')) +
  theme(axis.line = element_line(color = 'black', linewidth = 1)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1 )) +
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

# source(" ")

# data dir ----------------------------------------------------------------
origin_scRNA_dir <- " "
new_scRNA_celltype <- " "
output_dir <- " "


# hub ref -----------------------------------------------------------------
tumor_list <- c("ANS","BRCA","CESC","CRC","ESCC","GC","KIRC","LSCC","NSCLC","LUAD","LUSC","MA","NPC","PAAD","PNET",
                "SARC","TGCT","THCA","ACC","Liver_healthy","Lymph_node_healthy","Brain_healthy","HNSC")
# tumor_list <- c("Liver_healthy","Lymph_node_healthy","Brain_healthy","HNSC")
hub_cell_type <- c("Epithelia","Fibroblast","Endothelia","T_cell","NK_cell",
                   "B_cell","Plasma_cell","Myeloid_cell","Neutrophils")
merge_cell_type <- c("Epithelia","Fibroblast","Endothelia","Acinar_cell","T_cell","NK_cell",
                     "B_cell","Plasma_cell","Myeloid_cell","Neutrophils")

# all_new_cell_type <- data.frame(bind_id = as.character())

# mapping back to origin data ---------------------------------------------
for(tumor in tumor_list){
  # tumor <- "CRC"
  
  print(str_glue("############################## {tumor} Pre Processing ############################"))
  all_new_cell_type <- data.frame(bind_id = as.character())
  ## dir create
  output_file <- str_glue("{output_dir}/{tumor}/")
  unlink(output_file, recursive = T, force = T, expand = T)
  dir.create(output_file, recursive = T)
  
  ## scRNA origin data read
  scRNA_origin <- readRDS(str_glue("{origin_scRNA_dir}/{tumor}/scRNA_infercnv_combine.RData"))
  
  ## origin meta data
  origin_meta_data <- scRNA_origin@meta.data %>% 
    rownames_to_column(var = "bind_id")
  
  for(celltype in hub_cell_type){
    # celltype <- "T_cell"
    
    print(str_glue("          -------------- {celltype} back merge"))
    
    ## new cell type read
    new_meta_data <- read_csv(str_glue("{new_scRNA_celltype}/{celltype}/cell_type_meta.csv"))
    # scRNA_cell_type <- readRDS(str_glue("{new_scRNA_celltype}/{celltype}/scRNA_annotation.RData"))
    
    ## new meta data
    new_meta_data <- new_meta_data %>% 
      filter(tumor_code == tumor) %>% 
      dplyr::select(cell_id,new_cell_type)
    colnames(new_meta_data) <- c("cell_id","final_cell_type")
    
    ## combine all new cell type
    all_new_cell_type <- rbind(all_new_cell_type,new_meta_data)
    

  }
  
  ## combine and meta
  origin_meta_data <- origin_meta_data %>% 
    left_join(.,all_new_cell_type,by = "cell_id") %>% 
    column_to_rownames(var = "bind_id")
  origin_meta_data[is.na(origin_meta_data$final_cell_type),]$final_cell_type <- origin_meta_data[is.na(origin_meta_data$final_cell_type),]$new_cell_type
  
  ## add meta
  scRNA_origin@meta.data <- origin_meta_data
  
  ## remove undeifne cell
  origin_cell_type <- c("Epithelia","Fibroblast","Endothelia","T_cell","NK_cell","B_cell","Plasma_cell","Myeloid_cell","Neutrophils")
  keep_cell_type <- unique(scRNA_origin@meta.data$final_cell_type)[!unique(scRNA_origin@meta.data$final_cell_type) %in% origin_cell_type]
  scRNA_origin_hub <- subset(scRNA_origin,final_cell_type %in% keep_cell_type)
  
  ## data save
  saveRDS(scRNA_origin_hub,str_glue("{output_file}/scRNA_annotation_new.RData"))
  dior::write_h5(scRNA_origin_hub, file=str_glue("{output_file}/scRNA_annotation_new.h5"),
                 assay.name = "RNA",object.type = 'seurat')
  
  ## space free
  remove(scRNA_origin)
  gc()
  
}




