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
                     "Neuron"="#FB9A99",
                     "Macrophage"="#9ec9e1",
                     "Monocyte"="#E66F00",
                     "Dendritic_cell"="#1d92c0",
                     "Neutrophils"="#FFFF99",
                     "Mast_cell"="#6A3D9A",
                     "Neuron"="",
                     "Glial_cell"=""
)



# data dir ----------------------------------------------------------------
scRNA_dir <- " "
stRNA_dir <- " "
sample_dir <- " "
output_dir <- " "

# data read ---------------------------------------------------------------
sample_data <- read.xlsx(str_glue("{sample_dir} "),sheetName = "stRNA_seq") 
tumor_list <- list.files(stRNA_dir)

# data process ------------------------------------------------------------
for(tumor in tumor_list){
  # tumor <- "MA"
  print(str_glue("############################### {tumor} ###############################"))
  
  ## sample list
  input_file <- str_glue("{stRNA_dir}/{tumor}")
  file_name <- list.files(input_file)
  
  ## output file create
  output_file <- str_glue("{output_dir}/{tumor}")
  unlink(output_file, recursive = T, force = T, expand = T)
  dir.create(output_file, recursive = T)
  
  ## scRNA read
  if(tumor == "LUAD"){
    scRNA_data1 <- read_h5(file = str_glue("{scRNA_dir}/LUAD/scRNA_annotation_new.h5"),
                          assay.name = 'RNA', 
                          target.object = 'seurat')
    scRNA_data2 <- read_h5(file = str_glue("{scRNA_dir}/NSCLC/scRNA_annotation_new.h5"),
                          assay.name = 'RNA', 
                          target.object = 'seurat')
    scRNA_data <- merge(scRNA_data1,scRNA_data2)
    remove(scRNA_data1)
    remove(scRNA_data2)
    gc()
    
  }else{
    scRNA_data <- read_h5(file = str_glue("{scRNA_dir}/{tumor}/scRNA_annotation_new.h5"),
                          assay.name = 'RNA', 
                          target.object = 'seurat')
  }
  
  ## data transform
  scRNA_data@meta.data$final_cell_type <- gsub("^Ep_.*", "Epithelia", scRNA_data@meta.data$final_cell_type)
  
  ## cell type transform
  scRNA_data@meta.data <- scRNA_data@meta.data %>%
    mutate(most_cell_type = case_when(
      grepl("tumor", final_cell_type) ~ "Tumor_cell",
      grepl("^B_|^Plasma_", final_cell_type) ~ "B&Plasma_cell",
      grepl("^EC", final_cell_type) ~ "Endothelia",
      # grepl("^Plasma_", final_cell_type) ~ "Plasma_cell",
      grepl("^Fib_", final_cell_type) ~ "Fibroblast",
      grepl("^Epithelia", final_cell_type) ~ "Epithelia",
      grepl("^Melanoma", final_cell_type) ~ "Melanoma_cell",
      grepl("^Pericyte|^SMC", final_cell_type) ~ "SMC&Pericyte",
      grepl("^T_CD8", final_cell_type) ~ "T_CD8_cell",
      grepl("^T_CD4", final_cell_type) ~ "T_CD4_cell",
      grepl("^NK_", final_cell_type) ~ "NK_cell",
      grepl("^Mono_", final_cell_type) ~ "Monocyte",
      grepl("^Mph_", final_cell_type) ~ "Macrophage",
      grepl("^DC_", final_cell_type) ~ "Dendritic_cell",
      grepl("^Neu_", final_cell_type) ~ "Neutrophils",
      grepl("^Mast", final_cell_type) ~ "Mast_cell",
      grepl("^Melanoma_", final_cell_type) ~ "Melanoma_cell",
      grepl("^Neuron", final_cell_type) ~ "Neuron",
      grepl("^Glial_", final_cell_type) ~ "Glial_cell",
      # grepl("^Mast|^Myeloid|^Neu_", final_cell_type) ~ "Myeloid_Other",
      TRUE ~ "Other"
    ))
  scRNA_data <- subset(scRNA_data,most_cell_type != "Other")
  gc()
  
  
  ## RCTD process
  for(file in file_name[5:7]){
    # file <- "GSM5420750.RData"
    print(str_glue("     ----------------- {file}"))
    ## stRNA read
    if(gsub("\\.RData$","",file) %in% sample_data$sample_ID){
      stRNA_data <- readRDS(str_glue("{stRNA_dir}/{tumor}/{file}"))
    }
    
    ## scRNA data filter
    status <- as.character(sample_data[sample_data$sample_ID == gsub("\\.RData$","",file),]$tumor_status)
    meta_tissue <- as.character(sample_data[sample_data$sample_ID == gsub("\\.RData$","",file),]$Metastasis)
    
    # if(status == "Tumor"){
    #   scRNA_data_filter <- subset(scRNA_data,tumor_status == status)
    # }else if(status == "Metastasis" & tumor == "KIRC"){
    #   scRNA_data_filter <- subset(scRNA_data,tumor_status == status)
    # }else{
    #   scRNA_data_filter <- subset(scRNA_data,tumor_status == status)
    #   scRNA_data_filter <- subset(scRNA_data_filter,metastasis_tissue == meta_tissue)
    # }
    
    ## data filter
    if(status == "Metastasis"){
      hub_cell_type <- unique(scRNA_data@meta.data$most_cell_type)[!unique(scRNA_data@meta.data$most_cell_type) %in% c("Epithelia","Melanoma_cell")]
      scRNA_data_filter <- subset(scRNA_data, most_cell_type %in% hub_cell_type)
    }else if(status == "Tumor"){
      hub_cell_type <- unique(scRNA_data@meta.data$most_cell_type)[!unique(scRNA_data@meta.data$most_cell_type) %in% c("Glial_cell","Neuron")]
      scRNA_data_filter <- subset(scRNA_data, most_cell_type %in% hub_cell_type)
    }
    
    current_meta <- scRNA_data_filter@meta.data
    cell_num_stat <- as.data.frame(table(current_meta$most_cell_type))
    colnames(cell_num_stat) <- c("cell_type","num")
    cell_num_stat <- cell_num_stat[cell_num_stat$num >= 30,]
    scRNA_data_filter <- subset(scRNA_data_filter,most_cell_type %in% cell_num_stat$cell_type)
    gc()
    
    ## reference create
    counts <- scRNA_data_filter[["RNA"]]@counts
    cluster <- as.factor(scRNA_data_filter$most_cell_type)
    names(cluster) <- colnames(scRNA_data_filter)
    nUMI <- scRNA_data_filter$nCount_RNA
    names(nUMI) <- colnames(scRNA_data_filter)
    # reference <- Reference(counts, cluster, nUMI, n_max_cells = 5000)
    reference <- Reference(counts, cluster, n_max_cells = 10000)
    
    ## spatialRNA create
    counts <- stRNA_data[["spatial"]]@counts
    coords <- GetTissueCoordinates(stRNA_data,cols = c("imagerow", "imagecol"),scale = NULL)
    colnames(coords) <- c("x", "y")
    coords[is.na(colnames(coords))] <- NULL
    query <- SpatialRNA(coords = coords, counts = counts, colSums(counts))
    
    ## run RCTD
    RCTD <- create.RCTD(query, reference, max_cores = 80, CELL_MIN_INSTANCE = 20)
    gc()
    RCTD <- run.RCTD(RCTD, doublet_mode = "full")
    decon_mtrx <- as.matrix(RCTD@results$weights)
    decon_mtrx <- normalize_weights(decon_mtrx)
    
    stRNA_data$cell <- names(stRNA_data$orig.ident)
    stRNA_data <- subset(stRNA_data, cell %in% rownames(decon_mtrx))
    
    ## add RCTD result to stRNA
    stRNA_data <- AddMetaData(stRNA_data,as.data.frame(decon_mtrx))
    
    ## data plot
    local_data <- RCTD@spatialRNA@coords
    # pie_scale_size <- 650/nrow(local_data)
    plotSpatialScatterpie(local_data, as.data.frame(decon_mtrx), scatterpie_alpha = 1, pie_scale = 0.38, img = F, axis = "h") + 
      scale_fill_manual(values = cell_type_color)
    # plotSpatialScatterpie(local_data, as.data.frame(decon_mtrx), scatterpie_alpha = 1, pie_scale = pie_scale_size, img = F, axis = "h") + 
    #   scale_fill_manual(values = cell_type_color) + NoLegend()
    function_plot(filename_prefix = str_glue("{output_file}/{file}_plot"), width = 250, 230)
    
    ## data save
    sample_name <- gsub("\\.RData$","",file)
    saveRDS(stRNA_data,str_glue("{output_file}/{sample_name}_seurat.RData"))
    saveRDS(RCTD,str_glue("{output_file}/{sample_name}_RCTD.RData"))
    
  }
  ## space free
  remove(scRNA_data)
  gc()
}


