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
library("SCP")
(.packages())
setwd(' ')
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

# data color --------------------------------------------------------------
cell_type_color <- c("Epithelia_tumor"="#E41A1C",
                     "Melanoma_tumor"="#E41A1C",
                     "Osteoblastic_tumor"="#E41A1C",
                     "Epithelia"="#A6CEE3",
                     "Fibroblast"="#FDBF6F",
                     "Endothelia"="#B2DF8A",
                     "Acinar_cell"="#984EA3",
                     "SMC&Pericyte"="#33A02C",
                     "T_cell"="#1F78B4",
                     "NK_cell"="#64abc0",
                     "B_cell"="#FF7F00",
                     "Plasma_cell"="#FB9A99",
                     "Myeloid_cell"="#CAB2D6",
                     "Mast_cell"="#6A3D9A",
                     "Neutrophils"="#FFFF99",
                     "Hepatocyte"="#9ec9e1",
                     "Melanoma_cell"="#E66F00",
                     "Neuron"="#1d92c0",
                     "Glial_cell"="#B15928",
                     "Osteoblastic_cell"="#fcc5c1",
                     "Alpha_cell"="#42aa5e",
                     "Bela_cell"="#c22b86"
)

# data dir ----------------------------------------------------------------
origin_data_dir <- " "
infercnv_result_dir <- " "
output_dir <- " "


# data process ------------------------------------------------------------
cancer_file <- list.files(origin_data_dir)

for(file in cancer_file){
  # file <- "BRCA"
  print(str_glue("############################## {file} data combine ############################"))
  
  ## dir create
  output_file <- str_glue("{output_dir}/{file}")
  unlink(output_file, recursive = T, force = T, expand = T)
  dir.create(output_file, recursive = T)
  
  ## origin data read
  scRNA_origin <- read_h5(file = str_glue("{origin_data_dir}/{file}/scRNA_annotation_transform.h5"),
                         assay.name = 'RNA', 
                         target.object = 'seurat')
  
  ## infercnv result read
  if(file.exists(str_glue("{infercnv_result_dir}/{file}/scRNA_infercnv_meta.csv"))){
    infercnv_result <- read_csv(str_glue("{infercnv_result_dir}/{file}/scRNA_infercnv_meta.csv"))
    infercnv_result_filter <- infercnv_result %>% 
      dplyr::select(cell_id,cnv_leiden,cnv_score,cnv_status)
    meta_data <- scRNA_origin@meta.data %>% 
      rownames_to_column(var = "bind_id") %>% 
      left_join(.,infercnv_result_filter,by = "cell_id") %>% 
      column_to_rownames(var = "bind_id")
    meta_data$cnv_status[is.na(meta_data$cnv_status)] <- "other"
    # scRNA_origin@meta.data <- meta_data
    
  }else{
    meta_data <- scRNA_origin@meta.data %>% 
      rownames_to_column(var = "bind_id") %>% 
      column_to_rownames(var = "bind_id")
    meta_data$cnv_leiden <- 0
    meta_data$cnv_score <- 0
    meta_data$cnv_status <- "other"
    # scRNA_origin@meta.data <- meta_data
  }
  
  ## new cell type combine
  meta_data$cell_type <- as.character(meta_data$cell_type)
  meta_data_new <- meta_data %>% 
    mutate(new_cell_type = ifelse(tumor_status == "Metastasis" & cell_type == "Epithelia","Epithelia_tumor",
                           ifelse(tumor_status == "Metastasis" & cell_type == "Melanoma_cell","Melanoma_tumor",
                           ifelse(tumor_status == "Metastasis" & cell_type == "Osteoblastic_cell","Osteoblastic_tumor",
                           ifelse(tumor_status == "Tumor" & cnv_status == "tumor" & cell_type == "Epithelia","Epithelia_tumor",
                           ifelse(tumor_status == "Tumor" & cnv_status == "tumor" & cell_type == "Melanoma_cell","Melanoma_tumor",
                           ifelse(tumor_status == "Tumor" & cnv_status == "tumor" & cell_type == "Osteoblastic_cell","Osteoblastic_tumor",
                                  cell_type)))))))
  scRNA_origin@meta.data <- meta_data_new
  
  ## data plot
  DimPlot(scRNA_origin, reduction = "umap",group.by = "new_cell_type",raster=FALSE,cols = cell_type_color,pt.size = nrow(meta_data_new)/900000) +
    cluster_theme + theme_void() + NoLegend() +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),axis.ticks.y = element_blank()) + labs(title = NULL)
  function_plot(filename_prefix = str_glue("{output_file}/UMAP_infercnv_combine"), width = 300, height = 300)
  
  ## data save
  saveRDS(scRNA_origin,str_glue("{output_file}/scRNA_infercnv_combine.RData"))
  dior::write_h5(scRNA_origin, file=str_glue("{output_file}/scRNA_infercnv_combine.h5"),
                 assay.name = "RNA",object.type = 'seurat',save.scale = T)
}



