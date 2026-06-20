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
# library("scPDtools")
library("scRNAtoolVis")
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
cluster_color <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#F29403", "#F781BF", # 1-6
                   "#BC9DCC", "#A65628", "#54B0E4", "#222F75", "#1B9E77", "#B2DF8A", # 7-12
                   "#E3BE00", "#FB9A99", "#E7298A", "#910241", "#00CDD1", "#A6CEE3", # 13-18
                   "#CE1261", "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00", "#B3DE69", # 19-24
                   "#8DD3C7", "#999999", "#D87F32", "#737690", "#B17A7D", "#BD5E95"  # 25-30
                   )

# data dir ----------------------------------------------------------------
input_dir <- " "
output_dir <- " "
celltype <- "Endothelia"


# data read ---------------------------------------------------------------
scRNA_cell_type_merge <- read_h5(file = str_glue("{input_dir}/{celltype}/scRNA_cell_type_cluster.h5"),
                                 assay.name = 'RNA', 
                                 target.object = 'seurat')
Idents(scRNA_cell_type_merge) <- scRNA_cell_type_merge$leiden

DimPlot(scRNA_cell_type_merge, reduction = "umap",label = TRUE,raster=FALSE) + 
  cluster_theme + NoLegend() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),axis.ticks.y = element_blank())

FeaturePlot(scRNA_cell_type_merge,"CXCL12")

# new cell type annotation ------------------------------------------------

## dir create
output_file <- str_glue("{output_dir}/{celltype}")
unlink(output_file, recursive = T, force = T, expand = T)
dir.create(output_file, recursive = T)

## find marker
DEGs_marker <- FindAllMarkers(scRNA_cell_type_merge,only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
# DEGs_marker <- FindMarkers(scRNA_cell_type_merge,ident.1 = 0, min.pct = 0.25, logfc.threshold = 0.25)

## DEG filter
DEGs_marker_filter <- DEGs_marker %>%
  rownames_to_column(var = "gene_name") %>% 
  filter(pct.1 >= 0.3) %>% 
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 50) %>%
  ungroup()
# DEGs_marker_filter <- DEGs_marker_filter[order(match(DEGs_marker_filter$cluster, cluster.order)), ]
DEGs_marker_filter <- read_csv(str_glue("{output_dir}/{celltype}/DEGs_marker_again_filter.csv"))
DEGs_marker_again_filter <- DEGs_marker_filter


cluster.order <- c("0","1","27","2","12","5","11","15","3","8","4","6","24","25",
                  "23","9","10","17","13","14","18","19","20","21","22","26","28","29","7","16","30")

## marker plot
Cairo::CairoSVG(str_glue("{output_dir}/cluster_marker.pdf"),width = 6, height = 12)
averageHeatmap(object = scRNA_cell_type_merge, markerGene = DEGs_marker_filter$gene,
               gene.order = DEGs_marker_filter$gene,use_raster=F,
              #  cluster.order = cluster.order,
               myanCol = cluster_color)
dev.off()

## cell type annotation
new.cluster.ids <- c(
  "EC_01_ACKR1", # 0
  "EC_02_CXCL12", # 1
  "EC_01_ACKR1", # 2
  "EC_03_CD36", # 3
  "EC_04_PDGFD", # 4
  "EC_05_RGCC", # 5
  "EC_05_RGCC", # 6
  "undefine", # 7
  "EC_06_PROX1", # 8
  "EC_07_COL1A1", # 9
  "EC_08_DNAJB1", # 10
  "EC_09_ADAMTSL1", # 11
  "EC_10_FCN3", # 12
  "EC_11_SLC6A4", # 13
  "undefine", # 14
  "EC_12_PGF", # 15
  "EC_13_HPGD", # 16
  "undefine", # 17
  "EC_14_SLC7A5", # 18
  "EC_15_HSPA1B", # 19
  "EC_16_DKK2", # 20
  "EC_17_STMN1", # 21
  "EC_18_C7" # 22
  
)

scRNA_ano <- scRNA_cell_type_merge
names(new.cluster.ids) <- levels(scRNA_ano)
scRNA_ano <- RenameIdents(scRNA_ano, new.cluster.ids)
DimPlot(scRNA_ano, reduction = "umap",label = TRUE,raster=FALSE) + 
  cluster_theme + 
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),axis.ticks.y = element_blank())


## marker DEGs
DEGs_marker_again <- FindAllMarkers(scRNA_ano,only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
DEGs_marker_again_filter <- DEGs_marker_again %>% 
  rownames_to_column(var = "gene_name") %>% 
  filter(pct.1 >= 0.4) %>% 
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 50) %>%
  ungroup()

new.cluster.ids <- c(
  "Pericyte_01_HIGD1B", # Fib_01
  "Fib_01_RUNX2", # Fib_02
  "Fib_02_CFD", # Fib_03
  "Pericyte_02_HTR1F", # Fib_04
  "Fib_03_MMP11", # Fib_05
  "Fib_04_PI16", # Fib_06
  "undefine", # Fib_07
  "Fib_02_CFD", # Fib_08
  "SMC", # Fib_09
  "Fib_05_KIF26B", # Fib_10
  "Fib_02_CFD", # Fib_11
  "SMC", # Fib_12
  "Fib_06_POSTN", # Fib_13
  "undefine", # Fib_14
  "Fib_07_C3", # Fib_15
  "Fib_03_MMP11", # Fib_16
  "Pericyte_03_CCL21", # Fib_17
  "Fib_08_SOX6", # Fib_18
  "Fib_09_PTGDS", # Fib_19
  "undefine", # Fib_20
  "Fib_10_CLU", # Fib_21
  "Fib_11_STMN1", # Fib_22
  "undefine", # Fib_23
  "undefine", # Fib_24
  "Fib_09_PTGDS", # Fib_25
  "Fib_12_ROBO2", # Fib_26
  "Fib_13_NRXN1", # Fib_27
  "undefine" # Fib_28
)


names(new.cluster.ids) <- levels(scRNA_ano)
scRNA_ano <- RenameIdents(scRNA_ano, new.cluster.ids)

## join cell type in meta data
scRNA_ano@meta.data$new_cell_type <- Idents(scRNA_ano)
cell_type <- as.data.frame(table(scRNA_ano@meta.data$new_cell_type))
colnames(cell_type) <- c("cell_type","cell_num")
cell_type$origin_cell <- celltype

## remove undefine cell
scRNA_ano <- subset(scRNA_ano,new_cell_type != "undefine")

## cell type plot
if(T){
  DimPlot(scRNA_ano, reduction = "umap",label = TRUE,raster=FALSE) + 
    cluster_theme + 
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),axis.ticks.y = element_blank())
  function_plot(filename_prefix = str_glue("{output_file}/UMAP_annotation"), width = 200, height = 200)
  DimPlot(scRNA_ano, reduction = "umap",label = FALSE,raster=FALSE) + 
    cluster_theme + NoLegend() +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),axis.ticks.y = element_blank())
  function_plot(filename_prefix = str_glue("{output_file}/UMAP_annotation_nolabel"), width = 200, height = 200)
}

## data save
if(T){
  write_csv(DEGs_marker_again_filter,str_glue("{output_dir}/{celltype}/DEGs_marker_again_filter.csv"))
  write_csv(cell_type,str_glue("{output_dir}/{celltype}/cell_type.csv"))
  saveRDS(scRNA_ano,str_glue("{output_dir}/{celltype}/scRNA_annotation.RData"))
  dior::write_h5(scRNA_ano, file=str_glue("{output_dir}/{celltype}/scRNA_annotation.h5"),
                 assay.name = "RNA",object.type = 'seurat')

}
