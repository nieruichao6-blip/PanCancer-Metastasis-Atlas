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
library("Startrac")
library("corrplot")
library("Hmisc")
library("ComplexHeatmap")

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



# data dir ----------------------------------------------------------------
stRNA_dir <- " "
output_dir <- " "
RCTD_save <- data.frame(cell_type = as.character())

# CM  ---------------------------------------------------------------------
CM1 <- c("B_09_NR4A2","B_06_PDE4D","B_05_ANK3","B_12_RASSF6","B_14_CCSER1","B_03_CMSS1","T_CD4_03_ANK3","B_01_COL19A1","B_08_IGHM",
         "T_CD4_02_TSHZ2","Neu_06_CMTM2","Neu_02_CPPED1","Neu_01_S100A12","Neu_07_IFIT3")
CM2 <- c("EC_09_ADAMTSL1","EC_12_PGF","Fib_11_STMN1","Fib_06_POSTN","Fib_01_RUNX2","Fib_05_KIF26B","Mph_03_SLC16A10","Mph_10_HSPA1B",
         "EC_04_PDGFD","Pericyte_02_HTR1F","Mph_04_SELENOP_KCNMA1")
CM3 <- c("Plasma_01_IGHA2_IGHG3","Plasma_06_HSPA1B","Plasma_07_STMN1","T_STMN1","Fib_03_MMP11","T_CD4_Treg_02_TNFRSF18",
         "Mph_01_SELENOP_CCL18","T_CD8_08_HSPA1B","T_CD4_08_HSPA1B","Myeloid_STMN1","Mph_15_MT","Mph_07_SPP1")
CM4 <- c("Plasma_05_CD74","T_CD4_06_CXCL13","T_CD4_Treg_01_IKZF2","Mph_06_IL1B","Mph_09_ISG15","T_CD8_12_ISG15","B_16_ISG15",
         "T_CD4_10_ISG15","B_15_SOX5","T_CD8_01_CXCL13")
CM5 <- c("T_CD4_04_RTKN2","B_13_STMN1","B_10_RGS13","DC_02_CLEC9A","DC_01_CD1C","NK_05_SELL","T_CD8_09_CCR7","B_02_HSPA1B",
         "T_CD4_05_LTB","T_CD4_01_CCR7","B_04_TCL1A")
CM6 <- c("Plasma_03_CROCC","Plasma_02_SOX5","T_CD8_02_IL7R","Mono_02_IL1B","T_CD8_06_ANK3","T_CD8_05_ARIH1")
CM7 <- c("EC_01_ACKR1","EC_02_CXCL12","SMC","Fib_02_CFD","EC_03_CD36","EC_08_DNAJB1","Pericyte_03_CCL21","Fib_09_PTGDS",
         "Fib_04_PI16","EC_06_PROX1")
CM8 <- c("B_07_EGR1","T_CD4_09_ARIH1","Mast_cell","T_CD4_05_ANXA1","T_CD4_07_PPARG","T_CD8_07_P2RY8","NK_03_SPON2",
         "NK_04_FCGR3A_FGFBP2")
CM9 <- c("NK_01_CD160","T_CD8_10_SLC4A10","Mono_01_VCAN","Mono_03_FCGR3A","T_CD8_13_TRDV2","T_CD8_04_FGFBP2","NK_02_FCGR3A_CX3CR1")
all_CM <- c(CM1,CM2,CM3,CM4,CM5,CM6,CM7,CM8,CM9)
all_CM_list <- list("CM1"=CM1,"CM2"=CM2,"CM3"=CM3,"CM4"=CM4,"CM5"=CM5,"CM6"=CM6,"CM7"=CM7,"CM8"=CM8,"CM9"=CM9)


# data process ------------------------------------------------------------
tumor_list <- list.files(stRNA_dir)

for(tumor in tumor_list){
  # tumor <- "CRC"
  print(str_glue("############################### {tumor} ###############################"))
  
  ## output file create
  output_file <- str_glue("{output_dir}/{tumor}")
  unlink(output_file, recursive = T, force = T, expand = T)
  dir.create(output_file, recursive = T)
  
  ## sample list
  input_file <- str_glue("{stRNA_dir}/{tumor}")
  file_name <- list.files(input_file)
  
  
  ## RCTD process
  for(file in file_name){
    # file <- "GSM7058759.RData"
    print(str_glue("     ----------------- {file}"))
    ## stRNA read
    stRNA_data <- readRDS(str_glue("{stRNA_dir}/{tumor}/{file}"))
    
    ## CM combine
    for(CM_type in names(all_CM_list)){
      # CM_type <- "CM9"
      
      ### merge
      merge_list <- intersect(all_CM_list[[CM_type]], colnames(stRNA_data@meta.data))
      stRNA_data@meta.data[[CM_type]] <- rowSums(stRNA_data@meta.data[, merge_list, drop = FALSE])
      
      ### plot
      pt_size <- 4500/nrow(stRNA_data@meta.data)
      SpatialFeaturePlot(stRNA_data,features = CM_type,alpha = c(0.1, 1),pt.size.factor = pt_size)
      function_plot(filename_prefix = str_glue("{output_file}/{file}_{CM_type}"), width = 200, height = 220)
      
    }
    
    ## data save
    saveRDS(stRNA_data,str_glue("{output_file}/{file}"))

    
    
  }
  
}



