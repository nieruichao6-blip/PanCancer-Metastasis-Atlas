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
library("ggrastr")
(.packages())
setwd(' ')
show_palettes(palette_names = c("Paired","nejm","simspec","Spectral"))

function_plot <- function(filename_prefix, width, height){
  ggsave(filename = str_c(filename_prefix, '.pdf'), plot = last_plot(), device = 'pdf', width = width, height = height, units = 'mm', dpi = 500)
  ggsave(filename = str_c(filename_prefix, '.svg'), plot = last_plot(), device = 'svg', width = width, height = height, units = 'mm', dpi = 500)
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
cell_type_color <- c("tumor cell"="#E41A1C",
                     "Epithelia"="#00CDD1",
                     "Fibroblast"="#377EB8",
                     "Endothelia"="#4DAF4A",
                     "Acinar_cell"="#984EA3",
                     "Pericyte&SMC"="#F29403",
                     "T_cell"="#F781BF",
                     "NK_cell"="#FFDC91B2",
                     "B_cell"="#BC9DCC",
                     "Plasma_cell"="#A65628",
                     "Myeloid_cell"="#54B0E4",
                     "Mast_cell"="#222F75",
                     "Neutrophils"="#B2DF8A",
                     "Hepatocyte"="#1B9E77",
                     "Melanoma_cell"="#E3BE00",
                     "Neuron"="#FB9A99",
                     "Astrocyte"="#E7298A",
                     "Glial_cell"="#910241",
                     "Osteoblastic_cell"="#48dc88",
                     "Alpha_cell"="#C0AFE0",
                     "Bela_cell"="#FDC44D",
                     "undefine"="#A6CEE3"
)

tumor_type_color <- c("ACC"="#5A7B8F","ANS"="#EE4C97","BRCA"="#3D806F","CESC"="#A08634","CRC"="#F37C95","ESCC"="#608541",
                      "HNSC"="#7D4E57","KIRC"="#BC3C29","LUAD"="#958056","LUSC"="#9FAFA3","NPC"="#6F99AD","PAAD"="#0072B5",
                      "PNET"="#CFC59A","SARC"="#E18727","SKCM"="#FFDC91","GC"="#718DAE","TGCT"="#7876B1","THCA"="#F9AC93",
                      "NSCLC"="#E3BC06","LSCC"="#2285FE","Liver_healthy"="#3E6086","Lymph_node_healthy"="#20854E","Brain_healthy"="#4A7985")
tissue_type_color <- c("Bone_marrow"="#CE1261","Kidney"="#5E4FA2","Lymph_node"="#8CA77B","Brain"="#00441B","Lung"="#ed3d30",
                       "Colon"="#DEDC00","Liver"="#B3DE69","Peritoneal"="#999999","Stomach"="#725e82","Ovary"="#8c4356",
                       "Head_and_neck"="#db5a6b","Pancreas"="#21a675","Vaginal"="#d9b611","Appendix"="#e29c45",
                       "Bone"="#4c8dae","Nasopharynx"="#003371","Salivary_gland"="#a1afc9","Testis"="#88ada6","Cervix"="#c93756",
                       "Skin"="#204B75","Breast"="#588257","Thyroid"="#B6DB7B","Esophagus"="#4285BF","larynx"="#1CB232")
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
                          "GSE243639"="#6FCDDC")

"#E64B35B2" "#4DBBD5B2" "#00A087B2" "#3C5488B2" "#F39B7FB2" "#8491B4B2" "#91D1C2B2" "#DC0000B2" "#7E6148B2" "#B09C85B2"
"#3B4992B2" "#EE0000B2" "#008B45B2" "#631879B2" "#008280B2" "#BB0021B2" "#5F559BB2" "#A20056B2" "#808180B2" "#1B1919B2"
"#BC3C29B2" "#0072B5B2" "#E18727B2" "#20854EB2" "#7876B1B2" "#6F99ADB2" "#FFDC91B2" "#EE4C97B2"
"#00468BB2" "#ED0000B2" "#42B540B2" "#0099B4B2" "#925E9FB2" "#FDAF91B2" "#AD002AB2" "#ADB6B6B2" "#1B1919B2"
"#374E55B2" "#DF8F44B2" "#00A1D5B2" "#B24745B2" "#79AF97B2" "#6A6599B2" "#80796BB2"
"#2A6EBBB2" "#F0AB00B2" "#C50084B2" "#7D5CC6B2" "#E37222B2" "#69BE28B2" "#00B2A9B2" "#CD202CB2" "#747678B2"



# data dir ----------------------------------------------------------------
sample_idr <- " "
cancer_dir <- " "
output_dir <- " "


"#55d151", "#0e09b7", "#f69a95", "#aa329a", "#1d96a9", "#EDE816", "#48dc88", "#b15b0a", "#ac66d2", "#dc1c6f",
"#CB00FC", "#af3292", "#5da3ac", "#c6ab48", "#0c9986", "#2f52d1", "#aa152e", "#729e0c", "#be5793", "#3cc0c0", "#9ebe71", "#9a221c",
"#FD9A97", "#0c9c54", "#C0AFE0", "#D283FA", "#A30097", "#FDC44D", "#737D58", "#9e2d71", "#5A2AB4", "#0D6094", "#9F4700", "#16A795",
"#BF658B", "#59d185", "#FED8A3", "#7E9800", "#F82263", "#A35FFD", "#F65FFE", "#75ADCD", "#959099", "#FB5A8A", "#BD856E",
"#92628B", "#C5D0B7", "#C4B43B", "#E89B62", "#0DA8E7", "#822E35", "#5DE9EB", "#FEC4EC", "#8F00B3", "#1CB232", "#F668DE", "#9878CE",
"#A61C79", "#FC7216", "#76AE86", "#F87AA3", "#E3DDEF", "#CEE769", "#E7ABFD", "#53475A", "#0D56B2" 
"#E41B1B", "#4376AC", "#48A75A", "#87638F", "#D87F32", "#737690", "#D690C6", "#B17A7D", "#847A74",
 "#FA9B93", "#E9358B", "#A0094E", "#999999", "#6FCDDC", "#BD5E95"

all_cell <- c("Astrocyte","Oligodendrocyte","Neuron","Epithelia","Fibroblast","Endothelia",
              "T cell","B cell","Plasma cell","NK cell","Macrophages","Monocyte","Dendritic cell",
              "Mast cell","Smooth muscle cell","Neutrophils","Pericyte")

## stat collect
all_cell_type_num <- data.frame(
  cell_type = character(),
  cell_num = numeric(),
  tumor_type = character()
)

all_meta <- data.frame(
  sample_ID = character()
)

# data read ----------------------------------------------------------------
sample_data <- read_csv(str_glue("{sample_idr}/scRNA_data_sample_filter.csv"))
hub_cell_type <- c("Neutrophils","Endothelia","Fibroblast","B_cell","Plasma_cell","Myeloid_cell","Epithelia","T&NK_cell")

# run plot ----------------------------------------------------------------
cancer_file <- list.files(cancer_dir)

for(file in 1:length(cancer_file)){
  # file <- 4
  print(str_glue("############################## {cancer_file[file]} data plot ############################"))
  
  ## data read
  scRNA_annotation <- readRDS(str_glue("{cancer_dir}/{cancer_file[file]}/scRNA_annotation.RData"))
  cell_type_data <- read_csv(str_glue("{cancer_dir}/{cancer_file[file]}/cell_type.csv"))
  
  ## stat data
  current_meta <- scRNA_annotation@meta.data
  
  ## figure size
  if(nrow(current_meta) <= 80000){
    size = 300
  }else if(nrow(current_meta) >= 200000){
    size = 500
  }else{
    size = 400
  }

  ## dir create
  output_file <- str_glue("{output_dir}/{cancer_file[file]}")
  # unlink(output_file, recursive = T, force = T, expand = T)
  # dir.create(output_file, recursive = T)
  
  ## umap plot
  # if(T){
  #   DimPlot(scRNA_annotation, reduction = "umap",label = TRUE,raster=FALSE,cols = cell_type_color) +
  #     cluster_theme + NoLegend() +
  #     theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
  #           axis.ticks.x = element_blank(),axis.ticks.y = element_blank())
  #   function_plot(filename_prefix = str_glue("{output_file}/UMAP_annotation"), width = size, height = size)
  #   DimPlot(scRNA_annotation, reduction = "umap",label = FALSE,raster=FALSE,cols = cell_type_color) +
  #     cluster_theme + NoLegend() +
  #     theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
  #           axis.ticks.x = element_blank(),axis.ticks.y = element_blank())
  #   function_plot(filename_prefix = str_glue("{output_file}/UMAP_annotation_nolabel"), width = size, height = size)
  # }
 

  ### dot plot
  # current_cell_type <- unique(scRNA_annotation$cell_type)
  # y_direct <- subset(all_cell, all_cell %in% current_cell_type)
  # if(str_glue("{cancer_file[file]}") == "GBM"){
  #   DotPlot(scRNA_annotation, features = cell_marker_GBM) + nrc_theme + ylab("") + xlab("") +
  #     scale_y_discrete(limits=y_direct) +
  #     scale_size_continuous(range = c(0.1, 8),breaks = c(0,20,40,60,80,100)) +
  #     scale_color_gradient2(high = "#D7301F",mid = "#A1D99B",low = "#6BAED6",limits=c(-2,3))
  #   function_plot(filename_prefix = str_glue("{output_file}/markgene_exp"), width = 460, height = 20 + 12 * length(current_cell_type))
  # }else{
  #   DotPlot(scRNA_annotation, features = cell_marker) + nrc_theme + ylab("") + xlab("") +
  #     scale_y_discrete(limits=y_direct) +
  #     scale_size_continuous(range = c(0.1, 8),breaks = c(0,20,40,60,80,100)) +
  #     scale_color_gradient2(high = "#D7301F",mid = "#A1D99B",low = "#6BAED6",limits=c(-2,3))
  #   function_plot(filename_prefix = str_glue("{output_file}/markgene_exp"), width = 460, height = 20 + 12 * length(current_cell_type))
  # }

  ## tumor status plot
  # DimPlot(scRNA_annotation, reduction = "umap",group.by = "tumor_status",raster=FALSE,cols = tumor_status_color) +
  #   cluster_theme + NoLegend() +
  #   theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
  #         axis.ticks.x = element_blank(),axis.ticks.y = element_blank()) + labs(title = NULL)
  # function_plot(filename_prefix = str_glue("{output_file}/UMAP_tumor_status"), width = size, height = size)
  # 
  # ## all sample plot
  # DimPlot(scRNA_annotation, reduction = "umap",group.by = "sample_ID",raster=FALSE) +
  #   cluster_theme + NoLegend() +
  #   theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
  #         axis.ticks.x = element_blank(),axis.ticks.y = element_blank()) + labs(title = NULL)
  # function_plot(filename_prefix = str_glue("{output_file}/UMAP_sample_distribution"), width = size, height = size)
  
  ## cell type combine
  all_cell_type_num <- rbind(all_cell_type_num,cell_type_data)
  all_meta <- rbind(all_meta,current_meta)
}


# cell QC plot ------------------------------------------------------------
violin_data <- all_meta[,c("nCount_RNA","nFeature_RNA","percent.mt")]
violin_data_melt <- melt(violin_data)
colnames(violin_data_melt) <- c("QC_type","value")
# violin_data_melt <- violin_data_melt[violin_data_melt$QC_type == "nCount_RNA",]
nCount_RNA_plot <- ggplot(violin_data_melt[violin_data_melt$QC_type == "nCount_RNA",],aes(x = QC_type,y = value,fill = "red")) +
  geom_violin() + xlab("") + ylab("") + NoLegend() +
  theme(panel.background = element_rect(fill = "transparent")) + 
  theme(panel.border = element_rect(color = "transparent", fill='transparent')) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(axis.title = element_text(size = 15, face = 'bold')) +
  theme(axis.line = element_line(color = 'black', linewidth = 1)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5 )) +
  theme(axis.ticks = element_line(linewidth = 1), axis.ticks.length.y = unit("2", "mm")) +
  theme(axis.text = element_text(color="black", size = 12, face = 'bold')) +
  theme(legend.text = element_text(color="black", size = 12, face = 'bold'),
        legend.title = element_text(size = 15, face = 'bold'))
nFeature_RNA_plot <- ggplot(violin_data_melt[violin_data_melt$QC_type == "nFeature_RNA",],aes(x = QC_type,y = value,fill = "red")) +
  geom_violin() + xlab("") + ylab("") + NoLegend() +
  theme(panel.background = element_rect(fill = "transparent")) + 
  theme(panel.border = element_rect(color = "transparent", fill='transparent')) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(axis.title = element_text(size = 15, face = 'bold')) +
  theme(axis.line = element_line(color = 'black', linewidth = 1)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5 )) +
  theme(axis.ticks = element_line(linewidth = 1), axis.ticks.length.y = unit("2", "mm")) +
  theme(axis.text = element_text(color="black", size = 12, face = 'bold')) +
  theme(legend.text = element_text(color="black", size = 12, face = 'bold'),
        legend.title = element_text(size = 15, face = 'bold')) 
percent_mt_plot <- ggplot(violin_data_melt[violin_data_melt$QC_type == "percent.mt",],aes(x = QC_type,y = value,fill = "red")) +
  geom_violin() + xlab("") + ylab("") + NoLegend() +
  theme(panel.background = element_rect(fill = "transparent")) + 
  theme(panel.border = element_rect(color = "transparent", fill='transparent')) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(axis.title = element_text(size = 15, face = 'bold')) +
  theme(axis.line = element_line(color = 'black', linewidth = 1)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5 )) +
  theme(axis.ticks = element_line(linewidth = 1), axis.ticks.length.y = unit("2", "mm")) +
  theme(axis.text = element_text(color="black", size = 12, face = 'bold')) +
  theme(legend.text = element_text(color="black", size = 12, face = 'bold'),
        legend.title = element_text(size = 15, face = 'bold')) 
nCount_RNA_plot + nFeature_RNA_plot + percent_mt_plot
function_plot(filename_prefix = str_glue("{output_dir}/QC_after_plot"), width = 250, height = 150)

# cell type dis -----------------------------------------------------------
all_cell_type_num$exist_code <- 1
all_cell_type_num <- all_cell_type_num %>% 
  filter(cell_type != "undefine")

ggplot() +
  geom_tile(data=all_cell_type_num, aes(x=tumor_type, y=cell_type, fill=exist_code),color = "white",fill = "#DC0000B2") + 
  theme(panel.background = element_rect(fill = "transparent", colour = NA)) +
  theme(axis.title = element_text(size = 15, face = 'bold')) +
  theme(axis.text.x = element_text(color="black", size = 12, face = 'bold')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1 )) +
  theme(axis.ticks.x = element_line(size = 1), axis.ticks.length.x = unit("2", "mm")) + 
  theme(axis.text.y = element_text(color="black", size = 12, face = 'bold')) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
  scale_y_discrete(limits=c("T_cell","Fibroblast","Myeloid_cell","B_cell","NK_cell","Endothelia","Plasma_cell","Epithelia",
                            "SMC&Pericyte","Mast_cell","Neutrophils","Glial_cell","Hepatocyte","Neuron",
                            "Osteoblastic_cell","Acinar_cell","Bela_cell","Melanoma_cell")) +
  xlab('') + ylab('')
ggsave(filename = str_c(str_glue("{output_dir}/cell_type_exist_result"), '.pdf'), plot = last_plot(), device = 'pdf', width = 230, height = 210, units = 'mm', dpi = 500)
ggsave(filename = str_c(str_glue("{output_dir}/cell_type_exist_result"), '.svg'), plot = last_plot(), device = 'svg', width = 230, height = 210, units = 'mm', dpi = 500)
ggsave(filename = str_c(str_glue("{output_dir}/cell_type_exist_result"), '.png'), plot = last_plot(), device = 'png', width = 230, height = 210, units = 'mm', dpi = 1000)


# all meta tissue stat ----------------------------------------------------
write_csv(all_meta,str_glue("{output_dir}/all_meta.csv"))
# all_meta <- read_csv(str_glue("{output_dir}/all_meta.csv"))

## data filter
all_meta_filter <- all_meta %>% 
  filter(cell_type != "undefine")

all_tissue_cell_type_stat <- all_meta_filter %>%
  group_by(tissue_type, cell_type) %>%
  dplyr::summarise(count = n(), .groups = "drop")

all_tissue_type_stat <- as.data.frame(table(all_meta_filter$tissue_type))
colnames(all_tissue_type_stat) <- c("tissue_type","num")
all_tumor_status_stat <- as.data.frame(table(all_meta_filter$tumor_status))
colnames(all_tumor_status_stat) <- c("tumor_status","num")
all_cell_type_stat <- as.data.frame(table(all_meta_filter$cell_type))
colnames(all_cell_type_stat) <- c("cell_type","num")

## data save
if(T){
  write_csv(all_tissue_type_stat,str_glue("{output_dir}/all_tissue_type_stat.csv"))
  write_csv(all_tumor_status_stat,str_glue("{output_dir}/all_tumor_status_stat.csv"))
  write_csv(all_cell_type_stat,str_glue("{output_dir}/all_cell_type_stat.csv"))
}

# cell type plot ----------------------------------------------------------
write_csv(all_cell_type_num,str_glue("{output_dir}/all_cell_type_num.csv"))
# all_cell_type_num <- read_csv(str_glue("{output_dir}/all_cell_type_num.csv"))

ggplot(all_cell_type_num, aes(x = tumor_type, y = cell_num, fill = cell_type)) + 
  geom_bar(stat = "identity",position = "fill") +
  theme_bw() +
  labs(x = "", y = "Proportion of Cell Type(%)") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = cell_type_color) +
  nrc_theme
ggsave(filename = str_c(str_glue("{output_dir}/cell_proportion_result"), '.pdf'), plot = last_plot(), device = 'pdf', width = 250, height = 200, units = 'mm', dpi = 500)
ggsave(filename = str_c(str_glue("{output_dir}/cell_proportion_result"), '.svg'), plot = last_plot(), device = 'svg', width = 250, height = 200, units = 'mm', dpi = 500)
ggsave(filename = str_c(str_glue("{output_dir}/cell_proportion_result"), '.png'), plot = last_plot(), device = 'png', width = 250, height = 200, units = 'mm', dpi = 1000)


ggplot(all_tissue_cell_type_stat, aes(x = tissue_type, y = count, fill = cell_type)) + 
  geom_bar(stat = "identity",position = "fill") +
  theme_bw() +
  labs(x = "", y = "Proportion of Cell Type(%)") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = cell_type_color) +
  nrc_theme
ggsave(filename = str_c(str_glue("{output_dir}/tissue_cell_proportion_result"), '.pdf'), plot = last_plot(), device = 'pdf', width = 250, height = 200, units = 'mm', dpi = 500)
ggsave(filename = str_c(str_glue("{output_dir}/tissue_cell_proportion_result"), '.svg'), plot = last_plot(), device = 'svg', width = 250, height = 200, units = 'mm', dpi = 500)
ggsave(filename = str_c(str_glue("{output_dir}/tissue_cell_proportion_result"), '.png'), plot = last_plot(), device = 'png', width = 250, height = 200, units = 'mm', dpi = 1000)

# cell type proportion diff tissue --------------------------------------------
diff_status_counts <- all_meta %>%
  group_by(tumor_status, tumor_code, cell_type) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  arrange(tumor_status, tumor_code, desc(count))

metastasis_cell_type_num <- diff_status_counts %>% 
  filter(tumor_status == "Metastasis")
primary_cell_type_num <- diff_status_counts %>% 
  filter(tumor_status == "Tumor")
normal_cell_type_num <- diff_status_counts %>% 
  filter(tumor_status == "NAT")

## data plot
ggplot(metastasis_cell_type_num, aes(x = tumor_code, y = count, fill = cell_type)) + 
  geom_bar(stat = "identity",position = "fill") +
  theme_bw() +
  labs(x = "", y = "Proportion of Cell Type(%)") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = cell_type_color) +
  nrc_theme
ggsave(filename = str_c(str_glue("{output_dir}/metastasis_cell_proportion_result"), '.pdf'), plot = last_plot(), device = 'pdf', width = 250, height = 200, units = 'mm', dpi = 500)
ggsave(filename = str_c(str_glue("{output_dir}/metastasis_cell_proportion_result"), '.svg'), plot = last_plot(), device = 'svg', width = 250, height = 200, units = 'mm', dpi = 500)
ggsave(filename = str_c(str_glue("{output_dir}/metastasis_cell_proportion_result"), '.png'), plot = last_plot(), device = 'png', width = 250, height = 200, units = 'mm', dpi = 1000)

ggplot(primary_cell_type_num, aes(x = tumor_code, y = count, fill = cell_type)) + 
  geom_bar(stat = "identity",position = "fill") +
  theme_bw() +
  labs(x = "", y = "Proportion of Cell Type(%)") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = cell_type_color) +
  nrc_theme
ggsave(filename = str_c(str_glue("{output_dir}/primary_cell_proportion_result"), '.pdf'), plot = last_plot(), device = 'pdf', width = 250, height = 200, units = 'mm', dpi = 500)
ggsave(filename = str_c(str_glue("{output_dir}/primary_cell_proportion_result"), '.svg'), plot = last_plot(), device = 'svg', width = 250, height = 200, units = 'mm', dpi = 500)
ggsave(filename = str_c(str_glue("{output_dir}/primary_cell_proportion_result"), '.png'), plot = last_plot(), device = 'png', width = 250, height = 200, units = 'mm', dpi = 1000)

ggplot(normal_cell_type_num, aes(x = tumor_code, y = count, fill = cell_type)) + 
  geom_bar(stat = "identity",position = "fill") +
  theme_bw() +
  labs(x = "", y = "Proportion of Cell Type(%)") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = cell_type_color) +
  nrc_theme
ggsave(filename = str_c(str_glue("{output_dir}/normal_cell_proportion_result"), '.pdf'), plot = last_plot(), device = 'pdf', width = 250, height = 200, units = 'mm', dpi = 500)
ggsave(filename = str_c(str_glue("{output_dir}/normal_cell_proportion_result"), '.svg'), plot = last_plot(), device = 'svg', width = 250, height = 200, units = 'mm', dpi = 500)
ggsave(filename = str_c(str_glue("{output_dir}/normal_cell_proportion_result"), '.png'), plot = last_plot(), device = 'png', width = 250, height = 200, units = 'mm', dpi = 1000)





