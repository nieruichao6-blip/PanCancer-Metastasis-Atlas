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



# data dir ----------------------------------------------------------------
scRNA_dir <- " "
stRNA_dir <- " "
sample_dir <- " "
output_dir <- " "

RCTD_save <- data.frame(cell_type = as.character())

# data read ---------------------------------------------------------------
sample_data <- read.xlsx(str_glue("{sample_dir} "),sheetName = "stRNA_seq") 
tumor_list <- list.files(stRNA_dir)

# data process ------------------------------------------------------------
for(tumor in tumor_list[5:6]){
  # tumor <- "PAAD"
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
  
  
  ## RCTD process
  for(file in file_name){
    # file <- "P07.RData"
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
    
    if(status == "Metastasis" & meta_tissue == "Brain"){
      hub_cell_type <- unique(scRNA_data@meta.data$final_cell_type)[!unique(scRNA_data@meta.data$final_cell_type) %in% c("Epithelia","Melanoma_cell","Hepatocyte")]
      scRNA_data_filter <- subset(scRNA_data, final_cell_type %in% hub_cell_type)
    }else if(status == "Metastasis" & meta_tissue == "Liver"){
      hub_cell_type <- unique(scRNA_data@meta.data$final_cell_type)[!unique(scRNA_data@meta.data$final_cell_type) %in% c("Epithelia","Melanoma_cell","Glial_cell","Neuron")]
      scRNA_data_filter <- subset(scRNA_data, final_cell_type %in% hub_cell_type)
    }else if(status == "Metastasis" & meta_tissue == "Lymph_node"){
      hub_cell_type <- unique(scRNA_data@meta.data$final_cell_type)[!unique(scRNA_data@meta.data$final_cell_type) %in% c("Epithelia","Melanoma_cell","Glial_cell","Neuron","Hepatocyte")]
      scRNA_data_filter <- subset(scRNA_data, final_cell_type %in% hub_cell_type)
    }else if(status == "Tumor"){
      hub_cell_type <- unique(scRNA_data@meta.data$final_cell_type)[!unique(scRNA_data@meta.data$final_cell_type) %in% c("Glial_cell","Neuron","Hepatocyte")]
      scRNA_data_filter <- subset(scRNA_data, final_cell_type %in% hub_cell_type)
    }
    
    current_meta <- scRNA_data_filter@meta.data
    cell_num_stat <- as.data.frame(table(current_meta$final_cell_type))
    colnames(cell_num_stat) <- c("cell_type","num")
    cell_num_stat <- cell_num_stat[cell_num_stat$num >= 5,]
    scRNA_data_filter <- subset(scRNA_data_filter,final_cell_type %in% cell_num_stat$cell_type)
    gc()
    
    ## reference create
    counts <- scRNA_data_filter[["RNA"]]@counts
    cluster <- as.factor(scRNA_data_filter$final_cell_type)
    names(cluster) <- colnames(scRNA_data_filter)
    nUMI <- scRNA_data_filter$nCount_RNA
    names(nUMI) <- colnames(scRNA_data_filter)
    reference <- Reference(counts, cluster, nUMI, n_max_cells = 10000)
    
    ## spatialRNA create
    counts <- stRNA_data[["spatial"]]@counts
    coords <- GetTissueCoordinates(stRNA_data,cols = c("imagerow", "imagecol"),scale = NULL)
    colnames(coords) <- c("x", "y")
    coords[is.na(colnames(coords))] <- NULL
    query <- SpatialRNA(coords = coords, counts = counts, colSums(counts))
    
    ## run RCTD
    RCTD <- create.RCTD(query, reference, max_cores = 80, CELL_MIN_INSTANCE = 5)
    gc()
    RCTD <- run.RCTD(RCTD, doublet_mode = "full")
    decon_mtrx <- as.matrix(RCTD@results$weights)
    decon_mtrx <- as.data.frame(normalize_weights(decon_mtrx))
    decon_mtrx$sample_ID <- gsub("\\.RData$","",file)
    decon_mtrx$tumor_status <- status
    current_save <- decon_mtrx %>% 
      t(.) %>% 
      as.data.frame(.) %>% 
      rownames_to_column(var = "cell_type")
    
    stRNA_data$cell <- names(stRNA_data$orig.ident)
    stRNA_data <- subset(stRNA_data, cell %in% rownames(decon_mtrx))
    
    ## add RCTD result to stRNA
    stRNA_data <- AddMetaData(stRNA_data,as.data.frame(decon_mtrx))
    
    # SpatialFeaturePlot(stRNA_data, features = c("Epithelia_tumor"))
    
    ## all meta data collect
    RCTD_save <- full_join(RCTD_save,current_save,by = "cell_type")
    
    ## data save
    saveRDS(stRNA_data,str_glue("{output_file}/{file}"))
    
    ## space free
    remove(scRNA_data_filter)
    remove(reference)
    
  }
  ## space free
  remove(scRNA_data)
  gc()
}

## all data save
RCTD_save_transform <- RCTD_save %>% 
  column_to_rownames(var = "cell_type") %>% 
  t(.) %>% 
  as.data.frame(.) %>% 
  rownames_to_column(var = "cell_id")
write_csv(RCTD_save_transform,str_glue("{output_dir}/RCTD_save_transform.csv"))



