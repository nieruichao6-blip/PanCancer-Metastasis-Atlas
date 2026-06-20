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
# library('xlsx')
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
# library("scPDtools")
library("liana")
library("ComplexHeatmap")
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

# data color --------------------------------------------------------------
CM_color <- c("CM1"="#ea5c6f","CM2"="#f7905a","CM3"="#3187cb","CM4"="#fb948d","CM5"="#e2b159",
              "CM6"="#ebed6f","CM7"="#b2db87","CM8"="#7ee7bb","CM9"="#64cccf")
tumor_type_color <- c("ACC"="#5A7B8F","ANS"="#EE4C97","BRCA"="#3D806F","CESC"="#A08634","CRC"="#F37C95","ESCC"="#608541",
                      "HNSC"="#7D4E57","KIRC"="#BC3C29","LUAD"="#958056","LUSC"="#9FAFA3","NPC"="#6F99AD","PAAD"="#0072B5",
                      "PNET"="#CFC59A","SARC"="#E18727","SKCM"="#FFDC91","GC"="#718DAE","TGCT"="#7876B1","THCA"="#F9AC93","MA"="#FDC44D",
                      "NSCLC"="#E3BC06","LSCC"="#2285FE","Liver_healthy"="#3E6086","Lymph_node_healthy"="#20854E","Brain_healthy"="#4A7985")
tissue_type_color <- c("Bone_marrow"="#CE1261","Kidney"="#5E4FA2","Lymph_node"="#8CA77B","Brain"="#00441B","Lung"="#ed3d30",
                       "Colon"="#DEDC00","Liver"="#B3DE69","Peritoneal"="#999999","Stomach"="#725e82","Ovary"="#8c4356",
                       "Head_and_neck"="#db5a6b","Pancreas"="#21a675","Vaginal"="#d9b611","Appendix"="#e29c45",
                       "Bone"="#4c8dae","Nasopharynx"="#003371","Salivary_gland"="#a1afc9","Testis"="#88ada6","Cervix"="#c93756",
                       "Skin"="#204B75","Breast"="#588257","Thyroid"="#B6DB7B","Esophagus"="#4285BF","Larynx"="#1CB232")
tumor_status_color <- c("Tumor" = "#5BC0EB","Metastasis"="#8DD3C7","NAT" = "#9BC53D","Healthy" = "#C3423F")
spot_color <- c("malignant"="#E41A1C","normal"="#48dc88")
cellchat_color <- c("Secreted Signaling"="#FF7F0E","ECM-Receptor"="#D62728","Cell-Cell Contact"="#2CA02C")

# data dir ----------------------------------------------------------------
metastasis_tumor_dir <- " "
CM_dir <- " "
output_dir <- " "


# data process ------------------------------------------------------------
# for(CM_type in c("CM1","CM2","CM3","CM4","CM5","CM6","CM7","CM8","CM9")){
for(CM_type in c("CM7","CM8","CM9")){
  # CM_type <- "CM5"
  
  ## CM read
  CM_data <- read_h5(file = str_glue("{CM_dir}/scRNA_{CM_type}_origin.h5"),
                     assay.name = 'RNA', 
                     target.object = 'seurat')
  CM_data@meta.data$type <- CM_data@meta.data$final_cell_type
  
  for(metastasis_tissue in c("Brain","Liver","Lymph_node","Other")){
    # metastasis_tissue <- "Lymph_node"
    
    print(str_glue("   ------------- {CM_type} -- {metastasis_tissue}"))
    
    ## metastasis tumor read
    tumor_data <- readRDS(str_glue("{metastasis_tumor_dir}/scRNA_{metastasis_tissue}.RData"))
    tumor_data@meta.data$type <- str_glue("{metastasis_tissue}_metastasis")
    
    ## data merge
    scRNA_list <- merge(x = CM_data,y=tumor_data)
    scRNA_list <- NormalizeData(scRNA_list, normalization.method = "LogNormalize", scale.factor = 10000)
    gc()
    
    ## liana process
    liana_result <- liana_wrap(scRNA_list, idents_col = 'type',
                               method  = 'cellphonedb', resource = 'OmniPath',
                               min_cells = 5)
    
    ## liana result filter
    liana_result_filter <- liana_result %>%
      filter(source != target) %>%
      filter(pvalue <= 0.05 & lr.mean >= 0.1)
    
    ## data save
    write_csv(liana_result_filter,str_glue("{output_dir}/interaction_{CM_type}_{metastasis_tissue}.csv"))
    
    ## space free
    remove(scRNA_list)
    remove(tumor_data)
    gc()
    
  }
  
}


 






