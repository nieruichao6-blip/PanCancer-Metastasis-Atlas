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
                       axis.text.x = element_text(angle = 0, hjust = 0.5 ),
                       panel.grid=element_blank(),
                       legend.text = element_text(face = "bold",size= 12),
                       legend.title= element_text(face = "bold",size= 12))


# data dir ----------------------------------------------------------------
input_dir <- " "
output_dir <- " "
CM_list <- c("CM1","CM2","CM3","CM4","CM5","CM6","CM7","CM8","CM9")
# CM_list <- c("CM4","CM5","CM6","CM7","CM8","CM9")

# data process ------------------------------------------------------------
for(CM_type in CM_list){
  print(str_glue("############################### {CM_type} process ###############################"))
  # CM_type <- "CM5"
  
  ## output file create
  output_file <- str_glue("{output_dir}/{CM_type}")
  unlink(output_file, recursive = T, force = T, expand = T)
  dir.create(output_file, recursive = T)
  
  ## data read
  scRNA_CM_data <- read_h5(file = str_glue("{input_dir}/scRNA_{CM_type}_origin.h5"),
                           assay.name = 'RNA', 
                           target.object = 'seurat')
  
  ## data nor
  scRNA_CM_data <- NormalizeData(scRNA_CM_data, normalization.method = "LogNormalize", scale.factor = 10000)
  
  for(status in c("All")){
  # for(status in c("All","NAT","Tumor","Metastasis")){
    status <- "NAT"
    print(str_glue("     --------------- {status} process"))
    
    all_commuication_stat <- data.frame(
      source = as.character(),
      target = as.character(),
      count = as.numeric()
    )
    
    ## data filter
    if(status == "All"){
      scRNA_CM_data_filter <- scRNA_CM_data
    }else{
      scRNA_CM_data_filter <- subset(scRNA_CM_data, tumor_status == status)
    }
    
    ## data process
    liana_result <- liana_wrap(scRNA_CM_data_filter, idents_col = 'final_cell_type',
                               method  = 'cellphonedb', resource = 'OmniPath',
                               min_cells = 5)
    
    ## liana result filter
    liana_result_filter <- liana_result %>% 
      filter(source != target) %>%
      # filter(pvalue <= 0.05) %>%
      filter(receptor.prop >= 0.3 & ligand.prop >= 0.3) %>% 
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
    breaks <- c(seq(0,60, length = 100))
    col.heat <- colorRampPalette(brewer.pal(n = 7, name = "RdBu")[4:7])
    cols <- col.heat(100)
    col1 <- colorRamp3(breaks, cols)
    
    svg(str_glue("{output_file}/{CM_type}_{status}_cell_type_communication.svg"),width = 7,height = 6)
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
    # function_plot(filename_prefix = str_glue("{output_file}/{CM_type}_cell_type_communication"), width = 120, height = 100)
    dev.off()
    
    ## ligand recepet stat
    if(status == "All"){
      liana_result_filter$communication_status <- str_c(liana_result_filter$ligand,liana_result_filter$receptor,sep = "_")
      lig_rec_stat <- as.data.frame(table(liana_result_filter$communication_status))
      colnames(lig_rec_stat) <- c("communication_status","num")
      lig_rec_stat <- arrange(lig_rec_stat,desc(num))
      lig_rec_stat <- head(lig_rec_stat,100)
      lig_rec_stat$num_correct <- lig_rec_stat$num / length(unique(liana_result_filter$source))
      
      ### plot
      ggplot(lig_rec_stat, aes(reorder(communication_status, - num_correct), num_correct)) + 
        geom_bar(stat = "identity", fill = "#2285FE") +
        theme_bw() +
        labs(x = "", y = "Number of ligand - receptor") +
        scale_y_continuous(expand = c(0.01,0)) +
        nrc_theme + NoLegend()
      function_plot(filename_prefix = str_glue("{output_file}/cell_type_hub_ligand_recepoter"), width = 450, height = 150)
        
    }
    
    
    ## data save
    write_csv(liana_result_filter,str_glue("{output_file}/{CM_type}_{status}_liana_result_filter.csv"))
    write_csv(all_commuication_stat,str_glue("{output_file}/{CM_type}_{status}_communication_num.csv"))
    
  }
  
  
}




