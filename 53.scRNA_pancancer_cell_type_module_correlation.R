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
library("Startrac")
library("corrplot")
library("Hmisc")
library("ComplexHeatmap")
(.packages())
# setwd(' ')
display.brewer.all()

function_plot <- function(filename_prefix, width, height, plot_obj){
  ggsave(filename = str_c(filename_prefix, '.pdf'), plot = last_plot(), device = 'pdf', width = width, height = height, units = 'mm', dpi = 500, plot = plot_obj)
  ggsave(filename = str_c(filename_prefix, '.svg'), plot = last_plot(), device = 'svg', width = width, height = height, units = 'mm', dpi = 500, plot = plot_obj)
  ggsave(filename = str_c(filename_prefix, '.png'), plot = last_plot(), device = 'png', width = width, height = height, units = 'mm', dpi = 2000, plot = plot_obj)
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

normalize_row <- function(x){
  (x - min(x))/(max(x) -min(x))
}

# data dir ----------------------------------------------------------------
sample_dir <- " "
input_dir <- " "
output_dir <- " "



# data read ---------------------------------------------------------------
sample_data <- read_csv(sample_dir)
hub_cell_sample_proportion <- read_csv(str_glue("{input_dir}/cell_in_cell_type_proportion_hub.csv"))
sample_cell_num <- read_csv(str_glue("{input_dir}/sample_proportion_hub.csv"))

## data combine
hub_cell_sample_proportion <- left_join(hub_cell_sample_proportion,sample_data,by = "sample_ID")

## hub sample filter
hub_sample <- unique(sample_cell_num[sample_cell_num$clust_total >= 500,]$sample_ID)

# cor analysis ------------------------------------------------------------
for(status in c("all","NAT","Tumor","Metastasis")){
  status <- "all"
  
  ## data prepare
  if(status == "all"){
    current_sample_proportion <- hub_cell_sample_proportion %>% 
      filter(sample_ID %in% hub_sample)
      # filter(tumor_status != "NAT") %>%
      # filter(tumor_status != "Healthy")
  }else{
    current_sample_proportion <- hub_cell_sample_proportion %>% 
      filter(sample_ID %in% hub_sample) %>% 
      filter(tumor_status == status)
  }
  
  current_sample_proportion_transform <- current_sample_proportion[,c(1,3,6)] %>%
    pivot_wider(names_from = new_cell_type,
                values_from = proportion,
                values_fill = 0) %>%
    column_to_rownames(var = "sample_ID")
  # current_sample_proportion_transform <- as.data.frame(t(apply(current_sample_proportion_transform,1,normalize_row)))
  
  ## high cell freq filter
  high_freq_cell <- melt(rownames_to_column(current_sample_proportion_transform,var = "sample_ID"))
  colnames(high_freq_cell) <- c("sample_ID","new_cell_type","prop")
  high_freq_cell <- high_freq_cell %>% 
    filter(prop > 0)
  high_freq_cell_stat <- as.data.frame(table(high_freq_cell$new_cell_type))
  colnames(high_freq_cell_stat) <- c("new_cell_type","num")
  
  ## cor cal
  cor_result <- rcorr(as.matrix(current_sample_proportion_transform),type = "pearson")
  cor_matrix <- as.matrix(cor_result$r)
  cor_pvalue <- as.matrix(cor_result$P)
  diag(cor_pvalue) <- 1
  
  ## other para
  breaks <- c(seq(min(cor_matrix), 0, length = 50), seq(0.01,max(cor_matrix)-0.2, length = 50))
  col.heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))
  cols <- col.heat(100)
  col1 <- colorRamp3(breaks, cols)
  
  ## label plot
  pdf(str_glue("{output_dir}/cell_type_cor_{status}_label.pdf"),width = 22,height = 20)
  Heatmap(
    cor_matrix,
    show_column_names = TRUE,
    show_row_names = TRUE,
    cluster_rows = T,
    cluster_columns = T,
    col = col1,
    color = cols,
    # breaks = seq(-0.4, 0.8, length.out = 100),
    # rect_gp = gpar(col = "black", lwd = .1),
    clustering_distance_rows = "euclidean",
    cluster_column_slices = "euclidean",
    # clustering_method_columns = "ward.D2",
    # clustering_method_rows = "ward.D2",
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
  
  ## nolabel plot
  pdf(str_glue("{output_dir}/cell_type_cor_{status}_nolabel.pdf"),width = 17,height = 15)
  Heatmap(
    cor_matrix,
    show_column_names = FALSE,
    show_row_names = FALSE,
    cluster_rows = T,
    cluster_columns = T,
    col = col1,
    color = cols,
    # breaks = seq(-0.4, 0.8, length.out = 100),
    # rect_gp = gpar(col = "black", lwd = .1),
    clustering_distance_rows = "euclidean",
    cluster_column_slices = "euclidean",
    # clustering_method_columns = "ward.D2",
    # clustering_method_rows = "ward.D2",
    clustering_method_columns = "complete",
    clustering_method_rows = "complete",
    cell_fun = function(j,i,x,y,w,h,fill){
      if(cor_pvalue[i,j] <= 0.01){
        grid.text("**",x,y,gp = gpar(fontsize = 13),vjust = 0.75,hjust = 0.5)
      }else if(cor_pvalue[i,j] <= 0.05){
        grid.text("*",x,y,gp = gpar(fontsize = 13),vjust = 0.75,hjust = 0.5)
      }
    }
    
  )
  dev.off()
  
  ## data save
  cor_pvalue_save <- rownames_to_column(as.data.frame(cor_pvalue),var = "cell_type")
  cor_matrix_save <- rownames_to_column(as.data.frame(cor_matrix),var = "cell_type")
  write_csv(cor_matrix_save,str_glue("{output_dir}/hub_cell_correlation.csv"))
  write_csv(cor_pvalue_save,str_glue("{output_dir}/hub_cell_pvalue.csv"))
}

