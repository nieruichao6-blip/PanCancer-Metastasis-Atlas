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
library("reticulate")
library("SCP")
(.packages())
display.brewer.all()


function_plot <- function(filename_prefix, width, height){
  ggsave(filename = str_c(filename_prefix, '.pdf'), plot = last_plot(), device = 'pdf', width = width, height = height, units = 'mm', dpi = 500)
  ggsave(filename = str_c(filename_prefix, '.png'), plot = last_plot(), device = 'png', width = width, height = height, units = 'mm', dpi = 2000)
  ggsave(filename = str_c(filename_prefix, '.svg'), plot = last_plot(), device = 'svg', width = width, height = height, units = 'mm', dpi = 500)
  try(ggsave(filename = str_c(filename_prefix, '.wmf'), plot = last_plot(), device = 'wmf', width = width, height = height, units = 'mm', dpi = 500),
      silent = T)
}


## plot theme
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



# data dir ----------------------------------------------------------------
input_dir <- " "
output_dir <- " "
group_data_dir <- " "
nFeature_RNA_min <- 500
nFeature_RNA_max <- 8000
nCount_RNA_min <- 1000
nCount_RNA_max <- 40000
percent_mt <- 15
percent_pHB <- 5

bar_color <- c("high_quality_cell" = "#4DAF4A", 
               "low_quality_cell" = "#E41A1C")


# read seurat object ------------------------------------------------------
cancer_file <- list.files(input_dir)

## read group data
group_data <- read.xlsx(str_glue("{group_data_dir} "),sheetName = "scRNA_seq")
colnames(group_data)[2] <- "sample_ID"

## stat collect
QC_data_stat <- data.frame(
  tumor_type = character(),
  origin_cell = numeric(),
  high_quality_cell = numeric()
)
all_sample_num <- data.frame(
  sample_ID = character(),
  num = numeric()
)

all_meta <- data.frame(
  sample_ID = character()
)
# QC_data_stat <- read_csv(str_glue("{output_dir}/QC_data_stat.csv"))

## run pre
for(file in 1:length(cancer_file)){
  # file <- 1
  print(str_glue("############################## {cancer_file[file]} Pre Processing ############################"))
  set.seed(1234)
  
  scRNA_list <- list()
  input_file <- str_glue("{input_dir}/{cancer_file[file]}")
  file_name <- list.files(input_file)
  
  ## dir create
  # output_file <- str_glue("{output_dir}/{cancer_file[file]}")
  # unlink(output_file, recursive = T, force = T, expand = T)
  # dir.create(output_file, recursive = T)
  
  ## seurat object read
  for(i in 1:length(file_name)){
    scRNA_list[[i]] <- readRDS(str_glue("{input_file}/{file_name[i]}"))
  }
  gc()
  
  ## data QC
  for(i in 1:length(scRNA_list)){
    sc <- scRNA_list[[i]]
    sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
    sc[["percent.pHB"]] <- PercentageFeatureSet(sc, pattern = "^HBA|HBB")
    sc[["percent.rb"]] <- PercentageFeatureSet(sc, pattern = "^RP[SL]")
    scRNA_list[[i]] <- sc
    
    current_meta <- sc@meta.data
    colnames(current_meta)[1] <- "sample_ID"
    all_meta <- rbind(all_meta,current_meta)
    
  }
  
  # ## data merge
  # scRNA_list <- merge(x = scRNA_list[[1]],y=scRNA_list[-1])
  # origin_cell <- nrow(scRNA_list@meta.data)

  
  ## data QC  
  # scRNA_list <- subset(scRNA_list, subset = nFeature_RNA > nFeature_RNA_min & 
  #                        nFeature_RNA < nFeature_RNA_max &
  #                        # percent.pHB < percent_pHB &
  #                        percent.mt < percent_mt &
  #                        nCount_RNA < nCount_RNA_max &
  #                        nCount_RNA > nCount_RNA_min
  #                      )
  # 
  # # join group data
  # scRNA_meta <- scRNA_list@meta.data
  # colnames(scRNA_meta)[1] <- "sample_ID"
  # scRNA_meta <- scRNA_meta %>%
  #   rownames_to_column(var = 'bind_id') %>%
  #   left_join(.,group_data,by = 'sample_ID') %>%
  #   column_to_rownames(var = 'bind_id')
  # scRNA_list@meta.data <- scRNA_meta
  # 
  # ## filter low cell sample
  # sample_num <- as.data.frame(table(scRNA_meta$sample_ID))
  # colnames(sample_num) <- c("sample_ID","num")
  # sample_num_filter <- sample_num %>% 
  #   filter(num >= 100)
  # scRNA_list <- subset(scRNA_list,sample_ID %in% sample_num_filter$sample_ID)
  # 
  # high_quality_cell <- nrow(scRNA_list@meta.data)
  # 
  # # ## seurat version transform
  # # scRNA_list[["RNA"]] <- as(scRNA_list[["RNA"]], "Assay")
  # 
  # ## data Statistics
  # new_data <- data.frame(
  #   tumor_type = str_glue("{cancer_file[file]}"),
  #   origin_cell = origin_cell,
  #   high_quality_cell = high_quality_cell
  # )
  # QC_data_stat <- rbind(QC_data_stat, new_data)
  # all_sample_num <- rbind(all_sample_num,sample_num_filter)
  # 
  # ## data save
  # saveRDS(scRNA_list,str_glue("{output_file}/scRNA_merge.RData"))
  # # saveRDS(scRNA_list@assays$RNA@counts,str_glue("{output_file}/scRNA_raw_count.RData"))
  # dior::write_h5(scRNA_list, file=str_glue("{output_file}/scRNA_merge.h5"),assay.name = "RNA",object.type = 'seurat')
  
  
  ## free space
  remove(scRNA_list)
  gc()
  
}


# QC plot -----------------------------------------------------------------
write_csv(all_meta,str_glue("{output_dir}/all_meta.csv"))

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
function_plot(filename_prefix = str_glue("{output_dir}/QC_before_plot"), width = 250, height = 150)



# sample filter -----------------------------------------------------------
group_data_filter <- group_data %>% 
  filter(sample_ID %in% all_sample_num$sample_ID)
write_csv(group_data_filter,str_glue("{group_data_dir}/scRNA_data_sample_filter.csv"))


# QC plot -----------------------------------------------------------------
QC_data_stat$low_quality_cell <- QC_data_stat$origin_cell - QC_data_stat$high_quality_cell
write_csv(QC_data_stat,str_glue("{output_dir}/QC_data_stat.csv"))

QC_data_stat <- QC_data_stat[,c(1,2,4,3)]
QC_plot_data <- melt(QC_data_stat[-2])
ggplot(QC_plot_data, aes(reorder(tumor_type, - value) ,value, fill = variable)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "", y = "Number of cells") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = bar_color) +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/QC_result"), width = 350, height = 200)
















