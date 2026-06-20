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
library("liana")
library("ComplexHeatmap")
library("Matrix")

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
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5 )) +
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
tumor_type_color <- c("ACC"="#5A7B8F","ANS"="#EE4C97","BRCA"="#3D806F","CESC"="#A08634","CRC"="#F37C95","ESCC"="#608541",
                      "HNSC"="#7D4E57","KIRC"="#BC3C29","LUAD"="#958056","LUSC"="#9FAFA3","NPC"="#6F99AD","PAAD"="#0072B5",
                      "PNET"="#CFC59A","SARC"="#E18727","SKCM"="#FFDC91","GC"="#718DAE","TGCT"="#7876B1","THCA"="#F9AC93",
                      "NSCLC"="#E3BC06","LSCC"="#2285FE","Liver_healthy"="#3E6086","Lymph_node_healthy"="#20854E","Brain_healthy"="#4A7985")
tissue_type_color <- c("Bone_marrow"="#CE1261","Kidney"="#5E4FA2","Lymph_node"="#8CA77B","Brain"="#00441B","Lung"="#ed3d30",
                       "Colon"="#DEDC00","Liver"="#B3DE69","Peritoneal"="#999999","Stomach"="#725e82","Ovary"="#8c4356",
                       "Head_and_neck"="#db5a6b","Pancreas"="#21a675","Vaginal"="#d9b611","Appendix"="#e29c45",
                       "Bone"="#4c8dae","Nasopharynx"="#003371","Salivary_gland"="#a1afc9","Testis"="#88ada6","Cervix"="#c93756",
                       "Skin"="#204B75","Breast"="#588257","Thyroid"="#B6DB7B","Esophagus"="#4285BF","Larynx"="#1CB232")
tumor_status_color <- c("Tumor" = "#5BC0EB","Metastasis"="#8DD3C7","NAT" = "#9BC53D","Healthy" = "#C3423F")
study_database_color <- c("GSE202813"="#55d151","GSE237425"="#0e09b7","GSE171306"="#f69a95","GSE152938"="#aa329a","GSE210038"="#1d96a9",
                          "GSE277742"="#EDE816","GSE143423"="#48dc88","GSE234832"="#b15b0a","GSE256490"="#ac66d2","GSE164789"="#dc1c6f",
                          "HRA006468"="#CB00FC","GSE178318"="#af3292","GSE183916"="#5da3ac","GSE231559"="#c6ab48","GSE261388"="#0c9986",
                          "GSE245552"="#2f52d1","GSE271913"="#aa152e","GSE201347"="#729e0c","GSE163558"="#be5793","GSE246662"="#3cc0c0",
                          "GSE181919"="#9ebe71","GSE173468"="#9a221c","GSE188737"="#FD9A97","GSE256136"="#0c9c54","GSE162708"="#C0AFE0",
                          "GSE156405"="#D283FA","GSE197177"="#A30097","GSE281288"="#FDC44D","GSE205013"="#737D58","GSE212966"="#9e2d71",
                          "GSE244031"="#5A2AB4","GSE270231"="#0D6094","GSE237070"="#9F4700","GSE249361"="#16A795","GSE280127"="#BF658B",
                          "GSE217084"="#59d185","GSE197778"="#FED8A3","GSE228501"="#7E9800","GSE208653"="#F82263","GSE174401"="#A35FFD",
                          "GSE169147"="#2285FE","GSE280584"="#F65FFE","GSE200218"="#75ADCD","GSE222446"="#959099","GSE225600"="#FB5A8A",
                          "GSE190870"="#BD856E","GSE184362"="#A61C79","GSE232237"="#FC7216","GSE199619"="#76AE86","GSE196756"="#F87AA3",
                          "GSE188900"="#E3DDEF","GSE206332"="#CEE769","GSE252490"="#E7ABFD","GSE243977"="#53475A","GSE185477"="#847A74",
                          "GSE243639"="#6FCDDC","GSE203067"="#1CB232")            



# data dir ----------------------------------------------------------------
input_dir <- " "
output_dir <- " "
CM_list <- c("CM1","CM2","CM3","CM4","CM5","CM6","CM7","CM8","CM9")


# data process ------------------------------------------------------------
for(CM_type in CM_list){
  print(str_glue("############################### {CM_type} process ###############################"))
  # CM_type <- "CM3"
  
  ## output file create
  output_file <- str_glue("{output_dir}/{CM_type}")
  # unlink(output_file, recursive = T, force = T, expand = T)
  # dir.create(output_file, recursive = T)
  
  ## data read
  scRNA_CM_data <- read_h5(file = str_glue("{input_dir}/scRNA_{CM_type}_origin.h5"),
                           assay.name = 'RNA', 
                           target.object = 'seurat')
  
  ## data nor
  scRNA_CM_data <- NormalizeData(scRNA_CM_data, normalization.method = "LogNormalize", scale.factor = 10000)
  
  for(status in c("all","NAT","Tumor","Metastasis")){
  # for(status in c("all")){
    # status <- "all"
    print(str_glue("     --------------- {status} process"))
    
    ## data filter
    if(status == "all"){
      scRNA_CM_data_filter <- scRNA_CM_data
    }else{
      scRNA_CM_data_filter <- subset(scRNA_CM_data, tumor_status == status)
    }
    
    
    
    ## create cellchat object
    cellchat <- createCellChat(object = scRNA_CM_data_filter, group.by = "final_cell_type")
    
    ## database
    CellChatDB <- CellChatDB.human
    # CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
    # cellchat@DB <- CellChatDB.use
    cellchat@DB <- CellChatDB.human
    
    ## pre analysis
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat,adj = PPI.human)
    # cellchat <- smoothData(cellchat, smooth_batch_size = 1000, parallelize_smoothing = TRUE, seed = 123, ncores = 20)
    
    cellchat@idents <- factor(cellchat@idents,levels = as.character(unique(cellchat@idents)))
    cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = FALSE, population.size = FALSE,nboot = 100)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    groupSize <- as.numeric(table(cellchat@idents))
    Cairo::CairoSVG(str_glue("{output_file}/{status}_netVisual_circle_count.svg"))
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                     weight.scale = T, label.edge= F, title.name = "Number of interactions")
    dev.off()
    Cairo::CairoSVG(str_glue("{output_file}/{status}_netVisual_circle_weight.svg"))
    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                     weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    dev.off()
    
    ## data save
    # cellchat <- readRDS(str_glue("{output_file}/all_cellchat_obj.RData"))
    saveRDS(cellchat,str_glue("{output_file}/{status}_cellchat_obj.RData"))

    
  }
  
  
}
netVisual_chord_gene(cellchat,signaling = "SPP1")

# diff communication strength ---------------------------------------------
tumor_status <- c("NAT","Tumor","Metastasis")
for(CM_type in CM_list){
  CM_type <- "CM3"
  print(str_glue("############################### {CM_type} process ###############################"))
  
  ## data read
  scRNA_cellchat <- list()
  for(status in tumor_status){
    scRNA_cellchat[[status]] <- readRDS(str_glue("{output_dir}/{CM_type}/{status}_cellchat_obj.RData"))
  }
  
  ## data merge
  scRNA_cellchat <- mergeCellChat(scRNA_cellchat,add.names = names(scRNA_cellchat))
  
  ## diff plot
  compareInteractions(scRNA_cellchat ,show.legend = F, measure = "weight") + nrc_theme
  function_plot(filename_prefix = str_glue("{output_dir}/{CM_type}/diff_plot_{CM_type}_bar_weight"), width = 70, height = 90)
  compareInteractions(scRNA_cellchat ,show.legend = F, measure = "count") + nrc_theme
  function_plot(filename_prefix = str_glue("{output_dir}/{CM_type}/diff_plot_{CM_type}_bar_count"), width = 70, height = 90)
  
  Cairo::CairoSVG(str_glue("{output_dir}/{CM_type}/diff_plot_{CM_type}_dotplot_count.svg"))
  netVisual_diffInteraction(scRNA_cellchat,weight.scale = T, measure = "count")
  dev.off()
  
  Cairo::CairoSVG(str_glue("{output_dir}/{CM_type}/diff_plot_{CM_type}_dotplot_weight.svg"))
  netVisual_diffInteraction(scRNA_cellchat,weight.scale = T, measure = "weight")
  dev.off()
  
  levels(scRNA_cellchat@idents$joint)
  Cairo::CairoSVG(str_glue("{output_dir}/diff_plot_dotplot_CM_to_tumor.svg"))
  netVisual_bubble(scRNA_cellchat,sources.use = c(), targets.use = c(), comparison = c(1,2),
                   max.dataset = 2, angle.x = 45, remove.isolate = T)
  dev.off()
  Cairo::CairoSVG(str_glue("{output_dir}/diff_plot_dotplot_tumor_to_CM.svg"))
  netVisual_bubble(scRNA_cellchat,sources.use = c(), targets.use = c(), comparison = c(1,2),
                   max.dataset = 2, angle.x = 45, remove.isolate = T)
  dev.off()
  
}



