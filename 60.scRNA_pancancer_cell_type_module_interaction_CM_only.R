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


# data dir ----------------------------------------------------------------
input_dir <- " "
tumor_cell_dir <- " "
output_dir <- " "

read_list <- c("CM1","CM2","CM3","CM4","CM5","CM6","CM7","CM8","CM9")


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

# LIANA analysis ----------------------------------------------------------
scRNA_list <- list()
all_commuication_stat <- data.frame(
  source = as.character(),
  target = as.character(),
  count = as.numeric()
)

## data read
for(i in 1:length(read_list)){
  # i <- 1
  scRNA_list[[i]] <- read_h5(file = str_glue("{input_dir}/scRNA_{read_list[i]}_origin.h5"),
                             assay.name = 'RNA', 
                             target.object = 'seurat')
  Idents(scRNA_list[[i]]) <- scRNA_list[[i]]$final_cell_type
  # scRNA_list[[i]] <- subset(scRNA_list[[i]], tumor_status == "Metastasis")
  scRNA_list[[i]] <- subset(scRNA_list[[i]], downsample = 5000)
  scRNA_list[[i]]@meta.data$type <- read_list[i]
}

## data merge
scRNA_list <- merge(x = scRNA_list[[1]],y=scRNA_list[-1])

## data nor
scRNA_list <- NormalizeData(scRNA_list, normalization.method = "LogNormalize", scale.factor = 10000)

## liana process
liana_result <- liana_wrap(scRNA_list, idents_col = 'type',
                           method  = 'cellphonedb', resource = 'OmniPath',
                           min_cells = 5)

## liana result filter
liana_result_filter <- liana_result %>%
  filter(source != target) %>%
  filter(pvalue <= 0.05 & lr.mean >= 0.3)

## connection stat
for(source_cell in unique(liana_result_filter$source)){
  
  for(target_cell in unique(liana_result_filter$target)){
    
    current_result_filter <- liana_result_filter %>% 
      filter(source == source_cell & target == target_cell)
    communication_num <- nrow(current_result_filter)
    
    ### communication stat collect
    current_stat <- data.frame(
      source = source_cell,
      target = target_cell,
      count = communication_num
    )
    all_commuication_stat <- rbind(all_commuication_stat,current_stat)
    
  }
}

## connection plot
all_commuication_stat_transform <- all_commuication_stat %>%
  pivot_wider(names_from = source,
              values_from = count,
              values_fill = 0) %>%
  column_to_rownames(var = "target")

## other para
breaks <- c(seq(0,250, length = 100))
col.heat <- colorRampPalette(brewer.pal(n = 7, name = "RdBu")[4:7])
cols <- col.heat(100)
col1 <- colorRamp3(breaks, cols)

svg(str_glue("{output_dir}/liana_cell_type_communication_CM.svg"),width = 6,height = 5)
Heatmap(
  as.matrix(all_commuication_stat_transform),
  show_column_names = TRUE,
  show_row_names = TRUE,
  cluster_rows = T,
  cluster_columns = T,
  col = col1,
  color = cols,
  clustering_distance_rows = "euclidean",
  cluster_column_slices = "euclidean"
  
)
dev.off()


# cellchat analysis -------------------------------------------------------

## create cellchat object
cellchat <- createCellChat(object = scRNA_list, group.by = "type")

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
cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = FALSE, population.size = FALSE,nboot = 1000)
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
Cairo::CairoSVG(str_glue("{output_dir}/netVisual_circle_count_CM.svg"))
netVisual_circle(cellchat@net$count, vertex.weight = 8000,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
Cairo::CairoSVG(str_glue("{output_dir}/netVisual_circle_weight_CM.svg"))
netVisual_circle(cellchat@net$weight, vertex.weight = 8000, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()


# diff l-g stat -----------------------------------------------------------
for(type in c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact")){
  # type <- "Secreted Signaling"
  
  CellChatDB <- CellChatDB.human
  CellChatDB.use <- subsetDB(CellChatDB, search = type)
  cellchat@DB <- CellChatDB.use
  # cellchat@DB <- CellChatDB.human
  
  ## pre analysis
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat,adj = PPI.human)
  # cellchat <- smoothData(cellchat, smooth_batch_size = 1000, parallelize_smoothing = TRUE, seed = 123, ncores = 20)
  
  cellchat@idents <- factor(cellchat@idents,levels = as.character(unique(cellchat@idents)))
  cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = FALSE, population.size = FALSE,nboot = 1000)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  ## count result
  interaction_count <- as.data.frame(rowSums(cellchat@net$count))
  interaction_count <- rownames_to_column(interaction_count,var = "CM")
  colnames(interaction_count)[2] <- "count"
  interaction_count <- arrange(interaction_count,desc(count))
  interaction_count$count <- interaction_count$count / 9
  
  ## data plot
  usecolor <- cellchat_color[[type]]
  ggplot(interaction_count,aes(x = reorder(CM, count), y = count,fill = usecolor)) +
    geom_bar(stat = "identity", width = 0.81) +
    # geom_text(y = 0.2, fontface = "bold",size = 4, hjust = 0) +
    nrc_theme +
    xlab("") + ylab("") +
    coord_flip()+ NoLegend()
  function_plot(filename_prefix = str_glue("{output_dir}/Cellchat_{type}"), width = 150, height = 160)
  
}


# cellchat result analsys ------------------------------------------------
cellchat <- readRDS(str_glue("{output_dir}/cellchat_all_CM.RData"))

## connect plot
breaks <- c(seq(0,0.8, length = 300))
col.heat <- colorRampPalette(brewer.pal(n = 7, name = "RdBu")[3:1])
cols <- col.heat(300)
col1 <- colorRamp3(breaks, cols)

svg(str_glue("{output_dir}/cellchat_CM_communication_weight.svg"),width = 6,height = 5)
Heatmap(
  as.matrix(cellchat@net$weight),
  show_column_names = TRUE,
  show_row_names = TRUE,
  cluster_rows = T,
  cluster_columns = T,
  col = col1,
  color = cols,
  clustering_distance_rows = "euclidean",
  cluster_column_slices = "euclidean"
  
)
dev.off()
# "#B2182B" "#EF8A62" "#FDDBC7" "#F7F7F7" "#D1E5F0" "#67A9CF" "#2166AC"
## data save
# cellchat <- readRDS(str_glue("{output_dir}/cellchat_all_CM.RData"))
saveRDS(cellchat,str_glue("{output_dir}/cellchat_all_CM.RData"))









