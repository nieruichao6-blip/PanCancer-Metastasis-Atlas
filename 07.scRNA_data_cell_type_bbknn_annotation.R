# environment -------------------------------------------------------------
rm(list = ls())
options(java.parameters = "-Xmx16384m")
options(stringsAsFactors = F)
options(warn=-1)
options(pkgType="source")
Sys.setenv(VROOM_CONNECTION_SIZE=500072)
path.sep <- .Platform$file.sep

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
library("dior")
library("sceasy")
library("SCP")

(.packages())
setwd(' ')
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
input_dir <- " "
output_dir <- " "
list.files(input_dir)
tumor_code <- "THCA"

[1] "ACC"                "ANS"                "brain_healthy"      "BRCA"               "CESC"               "CRC"               
[7] "ESCC"               "GC"                 "HNSC"               "KIRC"               "liver_healthy"      "LSCC"              
[13] "LUAD"               "LUSC"               "lymph_node_healthy" "MA"                 "NPC"                "NSCLC"             
[19] "PAAD"               "PNET"               "SARC"               "TGCT"               "THCA"


# data read and pre -------------------------------------------------------
scRNA_bbknn <- read_h5(file = str_glue("{input_dir}/{tumor_code}/scRNA_cluster.h5"),
                         assay.name = 'RNA', 
                         target.object = 'seurat')
# scRNA_bbknn$n_genes <- 0
Idents(scRNA_bbknn) <- scRNA_bbknn$leiden

## join cell id
scRNA_bbknn@meta.data$cell_id <- c(1:nrow(scRNA_bbknn@meta.data))
scRNA_bbknn@meta.data$cell_id <- str_c(tumor_code,scRNA_bbknn@meta.data$cell_id,sep = "_")

# cluster annotation ------------------------------------------------------

## dir create
output_file <- str_glue("{output_dir}/{tumor_code}")
unlink(output_file, recursive = T, force = T, expand = T)
dir.create(output_file, recursive = T)

## gene marker
table(scRNA_bbknn@meta.data$leiden)

DimPlot(scRNA_bbknn, reduction = "umap",label = TRUE,raster=FALSE) + 
  cluster_theme + NoLegend() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),axis.ticks.y = element_blank())

cell_marker <- c("EPCAM","KRT18","KRT19","KRT5","PRSS1","CPA1","FOXJ1","TG","AMBP",   # Epithelia # Acinar_cell
                 "DCN","COL1A1","COL1A2","COL3A1",   # Fibroblast
                 "PECAM1","VWF","PLVAP","RAMP2",  # Endothelia
                 "MYH11","CNN1","MYL9","TAGLN", # SMC
                 "RGS5","CD36","ACTA2","MT1M", # Pericyte
                 "CD3D","CD3E","CD2","NKG7","GNLY","KLRF1",     # T&NK_cell
                 "CD79A","CD79B","MS4A1",  # B_cell
                 "IGHA1","MZB1","IGHG1", # Plasma_cell
                 # "CD68","CD14","C1QC","C1QA","FCN1","LILR2","CD1C","CD1E",  # Myeloid cell
                 # "CD68","CD163","CD14","C1QC","C1QA","FCN1","FCGR3A","LST1","CD1C","CD1E",
                 "LYZ","CD14","CD16","CD68","CSF1R","FCER1G","C1QC","FCN1", # Myeloid_cell
                 "TPSAB1","CPA3","TPSB2",  # Mast_cell
                 "CSF3R","IL1R2","CXCL8", # Neutrophils
                 "ALB","FABP1","APOC3","CYP2E1", # Hepatocyte
                 "ALPL","RUNX2","CLEC11A",
                 "MAP2","STMN2","SRRM4",      # Neuron
                 "GFAP","AQP4","S100B",   # Astrocyte  Glial_cell
                 "MBP","PLP1","MOG","TF",   #  Oligodendrocyte
                 "MLANA","PMEL","DCT", # Melanoma_cell
                 "INS","MAFA","PDX1","GCG"
)

DotPlot(scRNA_bbknn, features = cell_marker,scale = T) + RotatedAxis()
DEGs_marker <- FindMarkers(scRNA_bbknn,ident.1 = 66, min.pct = 0.25, logfc.threshold = 0.25)
DEGs_marker_filter <- DEGs_marker %>% 
  rownames_to_column(var = "gene_name") %>% 
  arrange(desc(avg_log2FC))

## cell type annotation
new.cluster.ids <- c(
  "Epithelia", # 0
  "Epithelia", # 1
  "B_cell", # 2
  "Epithelia", # 3
  "Epithelia", # 4
  "B_cell", # 5
  "NK_cell", # 6
  "T_cell", # 7
  "Epithelia", # 8
  "T_cell", # 9
  "Epithelia", # 10
  "Myeloid_cell", # 11
  "Myeloid_cell", # 12
  "NK_cell", # 13
  "NK_cell", # 14
  "undefine", # 15
  "T_cell", # 16
  "T_cell", # 17
  "NK_cell", # 18
  "T_cell", # 19
  "T_cell", # 20
  "Myeloid_cell", # 21
  "Endothelia", # 22
  "Epithelia", # 23
  "SMC&Pericyte", # 24
  "Epithelia", # 25
  "Epithelia", # 26
  "T_cell", # 27
  "Myeloid_cell", # 28
  "NK_cell", # 29
  "T_cell", # 30
  "T_cell", # 31
  "T_cell", # 32
  "B_cell", # 33
  "Epithelia", # 34
  "T_cell", # 35
  "NK_cell", # 36
  "T_cell", # 37
  "Myeloid_cell", # 38
  "B_cell", # 39
  "T_cell", # 40
  "Myeloid_cell", # 41
  "NK_cell", # 42
  "Fibroblast", # 43
  "undefine", # 44
  "SMC&Pericyte", # 45
  "T_cell", # 46
  "T_cell", # 47
  "T_cell", # 48
  "Plasma_cell", # 49
  "T_cell", # 50
  "Endothelia", # 51
  "Fibroblast", # 52
  "undefine", # 53
  "Endothelia", # 54
  "undefine", # 55
  "Myeloid_cell", # 56
  "Myeloid_cell", # 57
  "undefine", # 58
  "T_cell", # 59
  "Endothelia", # 60
  "B_cell", # 61
  "undefine", # 62
  "undefine", # 63
  "undefine", # 64
  "Mast_cell" # 65
  
)

scRNA_ano <- scRNA_bbknn
names(new.cluster.ids) <- levels(scRNA_ano)
scRNA_ano <- RenameIdents(scRNA_ano, new.cluster.ids)

## join cell type in meta data
new_cell <- rownames_to_column(as.data.frame(new.cluster.ids),var = "leiden")
colnames(new_cell) <- c("leiden","cell_type")
meta_data <- scRNA_ano@meta.data %>% 
  rownames_to_column(var = 'id') %>% 
  left_join(.,new_cell,by = 'leiden') %>% 
  column_to_rownames(var = 'id')
scRNA_ano@meta.data <- meta_data
cell_type <- as.data.frame(table(scRNA_ano@meta.data$cell_type))
colnames(cell_type) <- c("cell_type","cell_num")
cell_type$tumor_type <- tumor_code

## remove undefine cell
scRNA_ano <- subset(scRNA_ano,cell_type != "undefine")

## cell type plot
### umap plot
if(T){
  DimPlot(scRNA_ano, reduction = "umap",label = TRUE,raster=FALSE) + 
    cluster_theme + NoLegend() +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),axis.ticks.y = element_blank())
  function_plot(filename_prefix = str_glue("{output_file}/UMAP_annotation"), width = 200, height = 200)
  DimPlot(scRNA_ano, reduction = "umap",label = FALSE,raster=FALSE) + 
    cluster_theme + NoLegend() +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),axis.ticks.y = element_blank())
  function_plot(filename_prefix = str_glue("{output_file}/UMAP_annotation_nolabel"), width = 200, height = 200)
}

## filter tumor sample cell
scRNA_tumor_ano <- subset(scRNA_ano, tumor_status == "Tumor",)

CellDimPlot(scRNA_ano,group.by = "cell_type",reduction = "UMAP",raster=FALSE,theme_use = ggplot2::theme_classic,
            label = TRUE, label_insitu = TRUE)

## data save
if(T){
  write_csv(cell_type,str_glue("{output_dir}/{tumor_code}/cell_type.csv"))
  saveRDS(scRNA_ano,str_glue("{output_dir}/{tumor_code}/scRNA_annotation.RData"))
  dior::write_h5(scRNA_ano, file=str_glue("{output_dir}/{tumor_code}/scRNA_annotation.h5"),
                 assay.name = "RNA",object.type = 'seurat',save.scale = T)
  saveRDS(scRNA_tumor_ano,str_glue("{output_dir}/{tumor_code}/scRNA_tumor_annotation.RData"))
  dior::write_h5(scRNA_tumor_ano, file=str_glue("{output_dir}/{tumor_code}/scRNA_tumor_annotation.h5"),
  assay.name = "RNA",object.type = 'seurat')
}

