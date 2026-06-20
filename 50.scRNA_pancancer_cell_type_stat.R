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
# library('xlsx')
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
library("ggrastr")
library("ComplexHeatmap")
library("ggalluvial")
# library("scPDtools")
(.packages())
# setwd(' ')
display.brewer.all()

function_plot <- function(filename_prefix, width, height){
  # ggsave(filename = str_c(filename_prefix, '.pdf'), plot = last_plot(), device = 'pdf', width = width, height = height, units = 'mm', dpi = 500)
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


# color -------------------------------------------------------------------
new_status_color <- c("NAT"="#9BC53D","Primary_tumor"="#5BC0EB",
                      "Brain_metastasis"="#C69D3B","Liver_metastasis"="#EFC971","Lymph_node_metastasis"="#48A83E","Other_metastasis"="#7ABC68",
                      "Brain_Healthy"="#D86127","Liver_Healthy"="#F47D6E","Lymph_node_Healthy"="#F49897")

cell_type_color <- c("Tumor_cell"="#E41A1C",
                     "Epithelia"="#A6CEE3",
                     "Fibroblast"="#FDBF6F",
                     "Endothelia"="#B2DF8A",
                     "T_CD4_cell"="#984EA3",
                     "SMC&Pericyte"="#33A02C",
                     "T_CD8_cell"="#1F78B4",
                     "NK_cell"="#64abc0",
                     "B_cell"="#FF7F00",
                     "Plasma_cell"="#FB9A99",
                     # "Myeloid_Other"="#FB9A99",
                     "Myeloid_cell"="#CAB2D6",
                     "Mast_cell"="#6A3D9A",
                     "Neutrophils"="#FFFF99",
                     "Macrophage"="#9ec9e1",
                     "Monocyte"="#E66F00",
                     "Dendritic_cell"="#1d92c0"
)
tumor_type_color <- c("ACC"="#5A7B8F","ANS"="#EE4C97","BRCA"="#3D806F","CESC"="#A08634","CRC"="#F37C95","ESCC"="#608541","MA"="#FDC44D",
                      "HNSC"="#7D4E57","KIRC"="#BC3C29","LUAD"="#958056","LUSC"="#9FAFA3","NPC"="#6F99AD","PAAD"="#0072B5",
                      "PNET"="#CFC59A","SARC"="#E18727","SKCM"="#FFDC91","GC"="#718DAE","TGCT"="#7876B1","THCA"="#F9AC93",
                      "NSCLC"="#E3BC06","LSCC"="#2285FE","Liver_healthy"="#3E6086","Lymph_node_healthy"="#20854E","Brain_healthy"="#4A7985")

# data dir ----------------------------------------------------------------
sample_dir <- " "
input_dir <- " "
output_dir <- " "


# data read ---------------------------------------------------------------
all_cell_meta <- read_csv(str_glue("{input_dir}/all_cell_combine_meta.csv"))
sample_data <- read_csv(sample_dir)

# meta data transform -----------------------------------------------------
## cell type confirm
all_cell_meta <- all_cell_meta %>%
  mutate(most_cell_type = case_when(
    grepl("tumor", final_cell_type) ~ "Tumor_cell",
    grepl("^B_|^Plasma_", final_cell_type) ~ "B&Plasma_cell",
    grepl("^EC", final_cell_type) ~ "Endothelia",
    # grepl("^Plasma_", final_cell_type) ~ "Plasma_cell",
    grepl("^Fib_", final_cell_type) ~ "Fibroblast",
    # grepl("^Ep_", final_cell_type) ~ "Epithelia",
    grepl("^Pericyte|^SMC", final_cell_type) ~ "SMC&Pericyte",
    grepl("^T_CD8", final_cell_type) ~ "T_CD8_cell",
    grepl("^T_CD4", final_cell_type) ~ "T_CD4_cell",
    grepl("^NK_", final_cell_type) ~ "NK_cell",
    grepl("^Mono_", final_cell_type) ~ "Monocyte",
    grepl("^Mph_", final_cell_type) ~ "Macrophage",
    grepl("^DC_", final_cell_type) ~ "Dendritic_cell",
    # grepl("^Neu_", final_cell_type) ~ "Neutrophils",
    grepl("^Mast|^Myeloid|^Neu_", final_cell_type) ~ "Myeloid_Other",
    TRUE ~ "Other"
  ))

## cell type filter
# all_cell_meta_filter <- all_cell_meta %>%
#   filter(most_cell_type != "Other")

## prop cal
all_cell_meta_stat <- all_cell_meta %>% 
  dplyr::group_by(sample_ID,most_cell_type) %>% 
  dplyr::summarise(count = n(), .groups = "drop") %>% 
  dplyr::group_by(sample_ID) %>% 
  dplyr::mutate(percentage = count / sum(count)*100) %>% 
  ungroup

## cell type filter
all_cell_meta_stat <- all_cell_meta_stat %>%
  filter(most_cell_type != "Other")

## hub sample with zero
all_cell_meta_stat_transform <- all_cell_meta_stat[,c(1,2,4)] %>%
  pivot_wider(names_from = most_cell_type,
              values_from = percentage,
              values_fill = 0)
all_cell_meta_stat_final <- melt(all_cell_meta_stat_transform)
colnames(all_cell_meta_stat_final) <- c("sample_ID","cell_type","percentage")
# write_csv(all_cell_meta_stat_transform,str_glue("{output_dir}/all_cell_sample_percentage.csv"))


## combine sample info
all_cell_meta_stat_final <- left_join(all_cell_meta_stat_final,all_cell_meta[!duplicated(all_cell_meta$sample_ID),c("sample_ID","tumor_status","tissue_type","metastasis_tissue")],by = "sample_ID")
all_cell_meta_stat_final <- all_cell_meta_stat_final %>% 
  mutate(metastasis_tissue_transform = ifelse(
    tumor_status == "Metastasis" & metastasis_tissue == "Brain","Brain_metastasis",
    ifelse(tumor_status == "Metastasis" & metastasis_tissue == "Liver","Liver_metastasis",
           ifelse(tumor_status == "Metastasis" & metastasis_tissue == "Lymph_node","Lymph_node_metastasis",
                  ifelse(tumor_status == "Tumor","Primary_tumor",
                         ifelse(tumor_status == "NAT","NAT",
                                ifelse(tumor_status == "Healthy" & tissue_type == "Brain","Brain_Healthy",
                                       ifelse(tumor_status == "Healthy" & tissue_type == "Liver","Liver_Healthy",
                                              ifelse(tumor_status == "Healthy" & tissue_type == "Lymph_node","Lymph_node_Healthy",
                                                     "Other_metastasis")))))))
                                              
  ))


# tumor cell analysis -----------------------------------------------------
sample_data <- all_cell_meta[,c("sample_ID","tumor_code")]
sample_data <- sample_data[!duplicated(sample_data$sample_ID),]
tumor_cell_meta_stat <- all_cell_meta_stat_final %>% 
  filter(cell_type == "Tumor_cell") %>% 
  filter(tumor_status %in% c("Tumor","Metastasis")) %>% 
  left_join(.,sample_data,by = "sample_ID")
tumor_cell_meta_stat$tumor_code_status <- str_c(tumor_cell_meta_stat$tumor_code,tumor_cell_meta_stat$tumor_status,sep = "_")

ggplot(tumor_cell_meta_stat, aes(tumor_code_status,percentage,color = tumor_code)) +
  geom_boxplot(aes(fill = tumor_code),outlier.shape = NA,alpha = 0.2) +
  geom_jitter_rast(shape = 16, position = position_jitter(0.2),alpha = 0.5,size = 3) +
  scale_color_manual(values = tumor_type_color) +
  scale_fill_manual(values = tumor_type_color) +
  scale_y_continuous(breaks = seq(0,max(tumor_cell_meta_stat$percentage),by = max(tumor_cell_meta_stat$percentage)/3),
                     labels = scales::number_format(accuracy = 10)) +
  nrc_theme + NoLegend() + xlab("")
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_tumor_cell_tumor_code_percentage"), width = 600, height = 120)

# data plot ---------------------------------------------------------------
for(celltype in as.character(unique(all_cell_meta_stat_final$cell_type))){
  # celltype <- "B&Plasma_cell"
  
  ## data filter
  current_meta_stat <- all_cell_meta_stat_final %>% 
    filter(cell_type == celltype)
  
  ## data plot
  # ymax <- ceiling(max(current_meta_stat$percentage)/10)*10
  ggplot(current_meta_stat, aes(metastasis_tissue_transform,percentage,color = metastasis_tissue_transform)) +
    geom_boxplot(aes(fill = metastasis_tissue_transform),outlier.shape = NA,alpha = 0.2) +
    geom_jitter_rast(shape = 16, position = position_jitter(0.2),alpha = 0.5,size = 3) +
    scale_color_manual(values = new_status_color) +
    scale_fill_manual(values = new_status_color) +
    scale_x_discrete(limits = c("NAT","Primary_tumor","Brain_metastasis","Liver_metastasis","Lymph_node_metastasis","Other_metastasis",
                                "Brain_Healthy","Liver_Healthy","Lymph_node_Healthy")) +
    stat_compare_means(comparisons = list( c("Primary_tumor", "Brain_metastasis"), c("Primary_tumor", "Liver_metastasis"),
                                           c("Primary_tumor", "Lymph_node_metastasis"),c("Primary_tumor", "Other_metastasis"),
                                           c("NAT","Primary_tumor"),
                                           c("Brain_Healthy", "Brain_metastasis"), c("Liver_Healthy", "Liver_metastasis"),
                                           c("Lymph_node_Healthy", "Lymph_node_metastasis"),c("NAT", "Other_metastasis")
                                           ),
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),symbols = c("***", "**", "*", "ns")),
                       method="wilcox.test",label = "p.signif") +
    scale_y_continuous(breaks = seq(0,max(current_meta_stat$percentage),by = max(current_meta_stat$percentage)/3),
                       labels = scales::number_format(accuracy = 10)) +
    nrc_theme + NoLegend() + xlab("")
  function_plot(filename_prefix = str_glue("{output_dir}/boxplot_diff_status/scRNA_{celltype}_diff_tissue_percentage"), width = 120, height = 330)
}



# bar plot all cell -------------------------------------------------------
all_cell_meta_barplot <- all_cell_meta %>%
  mutate(most_cell_type = case_when(
    grepl("tumor", final_cell_type) ~ "Tumor_cell",
    grepl("^B_|^Plasma_", final_cell_type) ~ "B_cell",
    grepl("^EC", final_cell_type) ~ "Endothelia",
    # grepl("^Plasma_", final_cell_type) ~ "Plasma_cell",
    grepl("^Fib_", final_cell_type) ~ "Fibroblast",
    grepl("^Ep_", final_cell_type) ~ "Epithelia",
    grepl("^Pericyte|^SMC", final_cell_type) ~ "SMC&Pericyte",
    grepl("^T_CD8", final_cell_type) ~ "T_CD8_cell",
    grepl("^T_CD4", final_cell_type) ~ "T_CD4_cell",
    grepl("^NK_", final_cell_type) ~ "NK_cell",
    grepl("^Mono_", final_cell_type) ~ "Monocyte",
    grepl("^Mph_", final_cell_type) ~ "Macrophage",
    grepl("^DC_", final_cell_type) ~ "Dendritic_cell",
    grepl("^Neu_", final_cell_type) ~ "Neutrophils",
    grepl("^Mast_", final_cell_type) ~ "Mast_cell",
    TRUE ~ "Other"
  ))
all_cell_meta_barplot_filter <- all_cell_meta_barplot %>% 
  filter(most_cell_type != "Other")

## diff status stat
diff_status_counts <- all_cell_meta_barplot_filter %>%
  dplyr::group_by(tumor_status, tumor_code, most_cell_type) %>%
  dplyr::summarise(cell_count = n(), .groups = "drop") %>%
  dplyr::group_by(tumor_status, tumor_code) %>% 
  dplyr::mutate(total_cells_in_group = sum(cell_count),
         relative_percentage = (cell_count/total_cells_in_group) * 100) %>% 
  ungroup() %>% 
  arrange(tumor_status, tumor_code, desc(relative_percentage))
write_csv(diff_status_counts,str_glue("{output_dir}/diff_tumor_cell_type_percentage.csv"))


for(status in c("NAT","Tumor","Metastasis")){
  status <- "Metastasis"
  
  ## hub cell filter
  current_cell_meta <- diff_status_counts %>% 
    filter(tumor_status == status)
  
  ## cell plot
  ggplot(current_cell_meta, aes(x = tumor_code, stratum = most_cell_type, alluvium = most_cell_type, y = relative_percentage, fill = most_cell_type)) +
    geom_stratum(width = 0.5, col=NA) +
    geom_alluvium(width = 0.4, alpha = 0.4) +
    scale_fill_manual(values = cell_type_color) +
    labs(x = '', y = 'Proportion of Cell Type(%)')+
    theme_classic()+ nrc_theme + NoLegend()
  function_plot(filename_prefix = str_glue("{output_dir}/barplot_all_cell/scRNA_{status}_barplot_self"), width = 220, height = 80)
}

# diff metastasis stat ----------------------------------------------------
all_cell_meta_metastasis <- all_cell_meta_barplot[all_cell_meta_barplot$tumor_status == "Metastasis",]
all_cell_meta_metastasis$metastasis_status <- str_c(all_cell_meta_metastasis$metastasis_tissue,all_cell_meta_metastasis$tumor_code,sep = "_")
all_cell_meta_metastasis <- all_cell_meta_metastasis%>%
  filter(most_cell_type != "Other")
all_cell_meta_metastasis_stat <- all_cell_meta_metastasis %>% 
  dplyr::group_by(metastasis_status,most_cell_type) %>% 
  dplyr::summarise(count = n(), .groups = "drop") %>% 
  dplyr::group_by(metastasis_status) %>% 
  dplyr::mutate(percentage = count / sum(count)*100) %>% 
  ungroup
write_csv(all_cell_meta_metastasis_stat,str_glue("{output_dir}/diff_metastasis_tumor_cell_type_percentage.csv"))

## data plot
ggplot(all_cell_meta_metastasis_stat, aes(x = metastasis_status, stratum = most_cell_type, alluvium = most_cell_type, y = percentage, fill = most_cell_type)) +
  geom_stratum(width = 0.5, col=NA) +
  geom_alluvium(width = 0.4, alpha = 0.4) +
  scale_fill_manual(values = cell_type_color) +
  labs(x = '', y = 'Proportion of Cell Type(%)')+
  theme_classic()+ nrc_theme + NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/barplot_all_cell/scRNA_all_metastasis_tumor_code"), width = 350, height = 140)

# all cell proportion in diff status --------------------------------------
## data stat
all_cell_stat <- all_cell_meta %>% 
  dplyr::group_by(sample_ID,final_cell_type) %>% 
  dplyr::summarise(count = n(), .groups = "drop") %>% 
  dplyr::group_by(sample_ID) %>% 
  dplyr::mutate(percentage = count / sum(count)*100) %>% 
  ungroup
# colnames(all_cell_stat) <- c("sample_ID","cell_type","count","percentage")

## hub sample with zero
all_cell_stat <- all_cell_stat[,c(1,2,4)] %>%
  pivot_wider(names_from = final_cell_type,
              values_from = percentage,
              values_fill = 0)
all_cell_stat <- melt(all_cell_stat)
colnames(all_cell_stat) <- c("sample_ID","cell_type","percentage")

## combine sample info
all_cell_stat <- left_join(all_cell_stat,all_cell_meta[!duplicated(all_cell_meta$sample_ID),c("sample_ID","tumor_status","tissue_type","metastasis_tissue")],by = "sample_ID")
all_cell_stat <- all_cell_stat %>% 
  mutate(metastasis_tissue_transform = ifelse(
    tumor_status == "Metastasis" & metastasis_tissue == "Brain","Brain_metastasis",
    ifelse(tumor_status == "Metastasis" & metastasis_tissue == "Liver","Liver_metastasis",
           ifelse(tumor_status == "Metastasis" & metastasis_tissue == "Lymph_node","Lymph_node_metastasis",
                  ifelse(tumor_status == "Tumor","Primary_tumor",
                         ifelse(tumor_status == "NAT","NAT",
                                ifelse(tumor_status == "Healthy" & tissue_type == "Brain","Brain_Healthy",
                                       ifelse(tumor_status == "Healthy" & tissue_type == "Liver","Liver_Healthy",
                                              ifelse(tumor_status == "Healthy" & tissue_type == "Lymph_node","Lymph_node_Healthy",
                                                     "Other_metastasis")))))))
    
  ))
all_cell_stat_filter <- all_cell_stat %>% 
  filter(metastasis_tissue_transform %in% c("NAT","Primary_tumor","Brain_metastasis","Liver_metastasis","Lymph_node_metastasis","Other_metastasis"))


## data plot
for(celltype in as.character(unique(all_cell_stat_filter$cell_type))){
  # celltype <- "Mph_07_SPP1"
  
  ## data filter
  current_meta_stat <- all_cell_stat_filter %>% 
    filter(cell_type == celltype)
  
  ## data plot
  # ymax <- ceiling(max(current_meta_stat$percentage)/10)*10
  ggplot(current_meta_stat, aes(metastasis_tissue_transform,percentage,color = metastasis_tissue_transform)) +
    geom_boxplot(aes(fill = metastasis_tissue_transform),outlier.shape = NA,alpha = 0.2) +
    geom_jitter_rast(shape = 16, position = position_jitter(0.2),alpha = 0.5,size = 3) +
    scale_color_manual(values = new_status_color) +
    scale_fill_manual(values = new_status_color) +
    scale_x_discrete(limits = c("NAT","Primary_tumor","Brain_metastasis","Liver_metastasis","Lymph_node_metastasis","Other_metastasis")) +
    stat_compare_means(comparisons = list( c("Primary_tumor", "Brain_metastasis"), c("Primary_tumor", "Liver_metastasis"),
                                           c("Primary_tumor", "Lymph_node_metastasis"),c("Primary_tumor", "Other_metastasis"),
                                           c("NAT","Primary_tumor"),c("NAT", "Other_metastasis")
    ),
    symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),symbols = c("***", "**", "*", "ns")),
    method="wilcox.test",label = "p.signif") +
    scale_y_continuous(breaks = seq(0,max(current_meta_stat$percentage),by = max(current_meta_stat$percentage)/3),
                       labels = scales::number_format(accuracy = 10)) +
    nrc_theme + NoLegend() + xlab("")
  function_plot(filename_prefix = str_glue("{output_dir}/boxplot_all_cell_type/scRNA_{celltype}_diff_tissue_percentage"), width = 90, height = 330)
}



# new cell type in tumor status -------------------------------------------
all_cell_meta_transform <- all_cell_meta %>%
  mutate(hub_cell_type = case_when(
    # grepl("^B_", new_cell_type) ~ "B_cell",
    grepl("^B_|^Plasma_", final_cell_type) ~ "B&Plasma_cell",
    grepl("^EC", final_cell_type) ~ "Endothelia",
    # grepl("^Plasma_", new_cell_type) ~ "Plasma_cell",
    grepl("^Fib_", final_cell_type) ~ "Fibroblast",
    # grepl("^Ep_", new_cell_type) ~ "Epithelia",
    grepl("^Pericyte|^SMC", final_cell_type) ~ "SMC&Pericyte",
    # grepl("^T_", new_cell_type) ~ "T_cell",
    # grepl("^NK_", new_cell_type) ~ "NK_cell",
    grepl("^T_|^NK_", final_cell_type) ~ "T&NK_cell",
    grepl("^Mono_|^Mph_|^DC_|^Neu_|^Mast|^Myeloid", final_cell_type) ~ "Myeloid_cell",
    # grepl("^Mono_", new_cell_type) ~ "Monocyte",
    # grepl("^Mph_", new_cell_type) ~ "Macrophage",
    # grepl("^DC_", new_cell_type) ~ "Dendritic_cell",
    # grepl("^Neu_", new_cell_type) ~ "Neutrophils",
    # grepl("^Mast|^Myeloid", new_cell_type) ~ "Other_Myeloid",
    
    TRUE ~ "Other"
  ))

celltype_proportions <- all_cell_meta_transform %>%
  dplyr::group_by(sample_ID, hub_cell_type, final_cell_type) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  dplyr::group_by(sample_ID, hub_cell_type) %>%
  dplyr::mutate(total_per_sample_celltype = sum(count),
                proportion = count / total_per_sample_celltype * 100) %>%
  ungroup()

celltype_proportions <- left_join(celltype_proportions,sample_data,by = "sample_ID")

celltype_proportions <- celltype_proportions %>% 
  mutate(metastasis_tissue_transform = ifelse(
    tumor_status == "Metastasis" & metastasis_tissue == "Brain","Brain_metastasis",
    ifelse(tumor_status == "Metastasis" & metastasis_tissue == "Liver","Liver_metastasis",
           ifelse(tumor_status == "Metastasis" & metastasis_tissue == "Lymph_node","Lymph_node_metastasis",
                  ifelse(tumor_status == "Tumor","Primary_tumor",
                         ifelse(tumor_status == "NAT","NAT",
                                ifelse(tumor_status == "Healthy" & tissue_type == "Brain","Brain_Healthy",
                                       ifelse(tumor_status == "Healthy" & tissue_type == "Liver","Liver_Healthy",
                                              ifelse(tumor_status == "Healthy" & tissue_type == "Lymph_node","Lymph_node_Healthy",
                                                     "Other_metastasis")))))))
    
  ))
all_cell_stat_filter <- celltype_proportions %>% 
  filter(metastasis_tissue_transform %in% c("NAT","Primary_tumor","Brain_metastasis","Liver_metastasis","Lymph_node_metastasis","Other_metastasis"))



## data plot


# data save ---------------------------------------------------------------















