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
library("ROGUE")
library("ggrastr")
library("slingshot")
library("CytoTRACE2")
library("spacexr")
library("SPOTlight")
library("SpatialExperiment")

(.packages())
# setwd(' ')
display.brewer.all()

function_plot <- function(filename_prefix, width, height){
  ggsave(filename = str_c(filename_prefix, '.pdf'), plot = last_plot(), device = 'pdf', width = width, height = height, units = 'mm', dpi = 500)
  # ggsave(filename = str_c(filename_prefix, '.svg'), plot = last_plot(), device = 'svg', width = width, height = height, units = 'mm', dpi = 500)
  ggsave(filename = str_c(filename_prefix, '.png'), plot = last_plot(), device = 'png', width = width, height = height, units = 'mm', dpi = 1000)
}

nrc_theme <- theme(panel.background = element_rect(fill = "transparent")) + 
  theme(panel.border = element_rect(color = "transparent", fill='transparent')) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(axis.title = element_text(size = 15, face = 'bold')) +
  theme(axis.line = element_line(color = 'black', linewidth = 1)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1 )) +
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


# color -------------------------------------------------------------------
cell_type_color <- c("Tumor_cell"="#E41A1C",
                     "Epithelia"="#A6CEE3",
                     "Melanoma_cell"="#CAB2D6",
                     "Fibroblast"="#FDBF6F",
                     "SMC&Pericyte"="#33A02C",
                     "Endothelia"="#B2DF8A",
                     "T_CD4_cell"="#984EA3",
                     "T_CD8_cell"="#1F78B4",
                     "NK_cell"="#64abc0",
                     "B&Plasma_cell"="#FF7F00",
                     # "Neuron"="#21a675",
                     "Macrophage"="#9ec9e1",
                     "Monocyte"="#E66F00",
                     "Dendritic_cell"="#1d92c0",
                     "Neutrophils"="#FFFF99",
                     "Mast_cell"="#6A3D9A",
                     "Neuron"="#DEDC00",
                     "Glial_cell"="#FB9A99"
)
spot_color <- c("malignant"="#E41A1C","normal"="#48dc88")



# data dir ----------------------------------------------------------------
stRNA_dir <- " "
sample_dir <- " "
# output_dir <- " "
output_dir <- " "

# data read ---------------------------------------------------------------
sample_data <- read.xlsx(str_glue("{sample_dir} "),sheetName = "stRNA_seq") 
tumor_list <- list.files(stRNA_dir)

# data process ------------------------------------------------------------
for(tumor in tumor_list){
  tumor <- "PAAD"
  print(str_glue("############################### {tumor} ###############################"))
  
  ## sample list
  input_file <- str_glue("{stRNA_dir}/{tumor}")
  st_file_name <- list.files(input_file)[grepl("seurat.RData",list.files(input_file))]
  RCTD_file_name <- list.files(input_file)[grepl("RCTD.RData",list.files(input_file))]
  
  ## output file create
  output_file <- str_glue("{output_dir}/{tumor}")
  # unlink(output_file, recursive = T, force = T, expand = T)
  # dir.create(output_file, recursive = T)
  
  for(file in 1:length(st_file_name)){
    file <- 3
    print(str_glue("     ----------------- {st_file_name[file]}"))
    
    ### data read
    stRNA_data <- readRDS(str_glue("{stRNA_dir}/{tumor}/{st_file_name[file]}"))
    RCTD_data <- readRDS(str_glue("{stRNA_dir}/{tumor}/{RCTD_file_name[file]}"))
    sample_name <- gsub("_seurat.RData$","",st_file_name[file])
    
    ### percenage plot
    decon_mtrx <- as.matrix(RCTD_data@results$weights)
    decon_mtrx <- normalize_weights(decon_mtrx)
    local_data <- RCTD_data@spatialRNA@coords
    plotSpatialScatterpie(local_data, as.data.frame(decon_mtrx), scatterpie_alpha = 1, pie_scale = 0.4, img = F) + 
      scale_fill_manual(values = cell_type_color) + scale_y_reverse() + NoLegend()
    function_plot(filename_prefix = str_glue("{output_file}/{sample_name}_cell_type_plot"), width = 230, 230)
    
    
    ### copykat result plot
    SpatialDimPlot(stRNA_data,group.by = "copykat_result",cols = spot_color,pt.size.factor = 1.45,stroke = NA) + NoLegend()
    function_plot(filename_prefix = str_glue("{output_file}/{sample_name}_copykat_plot"), width = 230, height = 230)
    
    ### cell type plot
    for(hub_cell_type in colnames(decon_mtrx)){
      # hub_cell_type <- "B&Plasma_cell"
      SpatialFeaturePlot(stRNA_data,features = hub_cell_type,alpha = c(1, 1),pt.size.factor = 1.45,stroke = NA) + NoLegend()
      function_plot(filename_prefix = str_glue("{output_file}/{sample_name}_{hub_cell_type}_nolegend_plot"), width = 230, 230)
    }
    
  }
  
}








