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
  # tumor <- "LUAD"
  print(str_glue("############################### {tumor} ###############################"))
  
  ## sample list
  input_file <- str_glue("{stRNA_dir}/{tumor}")
  file_name <- list.files(input_file)
  
 
  ## RCTD process
  for(file in file_name){
    # file <- "P07.RData"
    print(str_glue("     ----------------- {file}"))
    ## stRNA read
    stRNA_data <- readRDS(str_glue("{stRNA_dir}/{tumor}/{file}"))
    
    ## hub meta filter
    current_meta <- stRNA_data@meta.data %>% 
      dplyr::select(c("sample_ID","tumor_status",colnames(stRNA_data@meta.data)[colnames(stRNA_data@meta.data) %in% all_CM])) %>% 
      t(.) %>% 
      as.data.frame(.) %>% 
      rownames_to_column(var = "cell_type")
    
    ## all meta combine
    RCTD_save <- full_join(RCTD_save,current_meta,by = "cell_type")
  }
  
}


# CM correlation analysis -------------------------------------------------
for(CM_type in names(all_CM_list)){
  CM_type <- "CM9"
  
  ## data filter
  current_RCTD_result <- RCTD_save %>% 
    # filter()
    filter(cell_type %in% all_CM_list[[CM_type]])
    # dplyr::select(colnames(RCTD_save)[colnames(RCTD_save) %in% all_CM_list[[CM_type]]])
    # dplyr::select(all_CM_list[[CM_type]])
  
  ## cor analysis
  current_RCTD_result[is.na(current_RCTD_result)] <- 0
  current_RCTD_result <- current_RCTD_result %>% 
    column_to_rownames(var = "cell_type") %>% 
    t(.)
  cor_result <- rcorr(as.matrix(current_RCTD_result),type = "spearman")
  # cor_result <- rcorr(as.matrix(current_RCTD_result),type = "pearson")
  cor_matrix <- as.matrix(cor_result$r)
  cor_pvalue <- as.matrix(cor_result$P)
  diag(cor_pvalue) <- 1
  
  ## other para
  breaks <- c(seq(-0.6, 0, length = 50), seq(0.01,0.8, length = 50))
  col.heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))
  cols <- col.heat(100)
  col1 <- colorRamp3(breaks, cols)
  
  ## label plot
  # Cairo::CairoSVG(str_glue("{output_file}/{status}_netVisual_circle_weight.svg"))
  svg(str_glue("{output_dir}/cell_type_cor_{CM_type}_label_new.svg"),width = 7,height = 6)
  # pdf(str_glue("{output_dir}/cell_type_cor_{CM_type}_label.pdf"),width = 7,height = 6)
  Heatmap(
    cor_matrix,
    show_column_names = TRUE,
    show_row_names = TRUE,
    cluster_rows = T,
    cluster_columns = T,
    col = col1,
    color = cols,
    # breaks = seq(-0.5, 0.5, length.out = 100),
    # rect_gp = gpar(col = "black", lwd = .1),
    clustering_distance_rows = "euclidean",
    cluster_column_slices = "euclidean",
    clustering_method_columns = "complete",
    clustering_method_rows = "complete",
    cell_fun = function(j,i,x,y,w,h,fill){
      if(cor_pvalue[i,j] <= 0.01){
        grid.text("**",x,y,gp = gpar(fontsize = 14),vjust = 0.75,hjust = 0.5)
      }else if(cor_pvalue[i,j] <= 0.05){
        grid.text("*",x,y,gp = gpar(fontsize = 14),vjust = 0.75,hjust = 0.5)
      }
    }
    
  )
  dev.off()
}

