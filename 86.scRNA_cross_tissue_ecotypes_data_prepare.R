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
cell_num_dir <- " "
output_dir <- " "


# data read ---------------------------------------------------------------
## read
all_cell_num_stat <- read_csv(str_glue("{cell_num_dir}/cell_type_tumor_tissue_num.csv"))
all_meta <- read_csv(str_glue("{input_dir}/all_cell_combine_meta.csv"))

## data transform
all_meta$final_cell_type <- gsub("^Ep_.*", "Epithelia_cell", all_meta$final_cell_type)
all_meta$final_cell_type <- gsub("Epithelia_tumor|Melanoma_tumor|Osteoblastic_tumor", "Tumor_cell", all_meta$final_cell_type)

## tissue per cell type filter
cross_tissue_cell_num_stat <- all_cell_num_stat[all_cell_num_stat$tumor_code_stat <= 5 | all_cell_num_stat$tissue_type_stat <= 5,]
cross_tissue_cell_type <- cross_tissue_cell_num_stat$cell_type


# cal proportion ----------------------------------------------------------
all_meta <- all_meta %>% 
  filter(final_cell_type != "SMC&Pericyte") %>% 
  filter(final_cell_type != "Bela_cell") %>%
  filter(final_cell_type != "Acinar_cell") %>%
  filter(final_cell_type != "Melanoma_cell") %>%
  filter(final_cell_type != "Osteoblastic_cell") %>%
  filter(final_cell_type != "Neuron") %>%
  filter(final_cell_type != "Glial_cell") %>%
  filter(final_cell_type != "Hepatocyte") %>%
  filter(final_cell_type %in% all_meta$final_cell_type[!all_meta$final_cell_type %in% cross_tissue_cell_type])

cell_type_sample_proportion <- all_meta %>%
  group_by(sample_ID, final_cell_type) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(clust_total = sum(count)) %>%
  dplyr::mutate(clust_prop = count / clust_total * 100) %>% 
  arrange(sample_ID, desc(clust_prop))

cell_type_sample_proportion_hub <- cell_type_sample_proportion %>% 
  filter(clust_total >= 500) %>% 
  filter(final_cell_type %in% cell_type_sample_proportion$final_cell_type[!cell_type_sample_proportion$final_cell_type %in% cross_tissue_cell_type])


# data save ---------------------------------------------------------------
if(T){
  write_csv(cell_type_sample_proportion,str_glue("{output_dir}/cell_type_sample_proportion.csv"))
  write_csv(cell_type_sample_proportion_hub,str_glue("{output_dir}/cell_type_sample_proportion_hub.csv"))
  
  
}





