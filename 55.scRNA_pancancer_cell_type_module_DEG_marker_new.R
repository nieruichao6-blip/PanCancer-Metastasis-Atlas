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
  theme(axis.text.x = element_text(angle = 90, hjust = 1 )) +
  theme(axis.ticks = element_line(linewidth = 1), axis.ticks.length.y = unit("2", "mm")) +
  theme(axis.text = element_text(color="black", size = 12, face = 'bold')) +
  theme(legend.text = element_text(color="black", size = 12, face = 'bold'),
        legend.title = element_text(size = 15, face = 'bold')) 

cluster_theme <- theme(plot.title = element_text(face = "bold",size = 12,color="black",hjust = 0.5),
                       axis.title = element_text(face = "bold",size = 12,color ="black"), 
                       axis.text = element_text(face = "bold",size= 12,color = "black"),
                       panel.grid.minor.y = element_blank(),
                       panel.grid.minor.x = element_blank(),
                       axis.text.x = element_text(angle = 90, hjust = 0.5 ),
                       panel.grid=element_blank(),
                       legend.text = element_text(face = "bold",size= 12),
                       legend.title= element_text(face = "bold",size= 12))

# data dir ----------------------------------------------------------------
input_dir <- " "
output_dir <- " "

read_list <- c("CM1","CM2","CM3","CM4","CM5","CM6","CM7","CM8","CM9")



# data combine ------------------------------------------------------------
## combine
scRNA_list <-list()
for(i in 1:length(read_list)){
  
  ## data read
  scRNA_list[[i]] <- read_h5(file = str_glue("{input_dir}/scRNA_{read_list[i]}_origin.h5"),
                             assay.name = 'RNA', 
                             target.object = 'seurat')
  Idents(scRNA_list[[i]]) <- scRNA_list[[i]]$final_cell_type
  scRNA_list[[i]] <- subset(scRNA_list[[i]], downsample = 5000)
  scRNA_list[[i]]@meta.data$type <- read_list[i]
  
}

## data merge
scRNA_list <- merge(x = scRNA_list[[1]],y=scRNA_list[-1])

## data normal
scRNA_list <- NormalizeData(scRNA_list, normalization.method = "LogNormalize", scale.factor = 10000)

# marker find -------------------------------------------------------------
## find
Idents(scRNA_list) <- scRNA_list$type
DEGs_marker <- FindAllMarkers(scRNA_list,,only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

## filter
DEGs_marker_filter <- DEGs_marker %>%
  rownames_to_column(var = "gene_name") %>% 
  filter(pct.1 >= 0.5) %>% 
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 100) %>%
  ungroup()


# marker plot -------------------------------------------------------------
## CM marker
CM_marker <- c("BANK1","RIPOR2","RUBCNL","CAMK1D","PRKCB",
               "COL4A1","COL4A2","SPARC","PTPRG","FN1",
               "CTSD","CTSB","SSR4","XBP1","LGALS1",
               "ISG15","MX1","IFI6","EPSTI1","ISG20",
               "HLA-DPB1","HLA-DQB1","HLA-DQA1","LTB","LIMD2",
               "ANKRD28","AOAH","FNDC3A","GLCCI1","ATXN1",
               "IGFBP5","GSN","MGP","SPARCL1","TIMP3",
               "CD247","SYTL3","STAT4","TNIK","SKAP1",
               "KLRB1","CCL4","GZMA","KLRD1","HCST")

## data plot
DotPlot(scRNA_list, features = CM_marker) + nrc_theme + ylab("") + xlab("") +
  scale_y_discrete(limits=read_list) +
  scale_size_continuous(range = c(0.1, 8),breaks = c(0,20,40,60,80)) +
  scale_color_gradient2(high = "#D7301F",mid = "#A1D99B",low = "#6BAED6")
function_plot(filename_prefix = str_glue("{output_dir}/CM_marker_dotplot"), width = 400, height = 100)

# data save ---------------------------------------------------------------
if(T){
  write_csv(DEGs_marker,str_glue("{output_dir}/CM_DEGs_marker_all.csv"))
  write_csv(DEGs_marker_filter,str_glue("{output_dir}/CM_DEGs_marker_filter.csv"))
}


















