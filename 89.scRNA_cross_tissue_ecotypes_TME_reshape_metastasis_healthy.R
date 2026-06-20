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
library("ggalluvial")
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

# data color ----------------------------------------------------------------
metastasis_color <- c("Brain_Metastasis"="#33a02c","Vaginal_Metastasis"="#fb9a99","Liver_Metastasis"="#ff7f00","Lung_Metastasis"="#e31a1c",
                      "Lymph_node_Metastasis"="#a6761d","Ovary_Metastasis"="#b15928","Peritoneal_Metastasis"="#cab2d6","Bone_marrow_Metastasis"="#6a3d9a",
                      "Liver_Healthy"="#b2df8a","Lymph_node_Healthy"="#ffff99","Brain_Healthy"="#fdbf6f")
Ecotype_color <- c("ET_01"="#f97b72","ET_02"="#F6CF71","ET_03"="#3969AC","ET_04"="#80BA5A",
                   "ET_05"="#F2B701","ET_06"="#11A579","ET_07"="#CF1C90","ET_08"="#66C5CC")
tumor_type_color <- c("ACC"="#5A7B8F","ANS"="#EE4C97","BRCA"="#3D806F","CESC"="#A08634","CRC"="#F37C95","ESCC"="#608541",
                      "HNSC"="#7D4E57","KIRC"="#BC3C29","LUAD"="#958056","LUSC"="#9FAFA3","NPC"="#6F99AD","PAAD"="#0072B5",
                      "PNET"="#CFC59A","SARC"="#E18727","SKCM"="#FFDC91","GC"="#718DAE","TGCT"="#7876B1","THCA"="#F9AC93","MA"="#FDC44D",
                      "NSCLC"="#E3BC06","LSCC"="#2285FE","Liver_healthy"="#3E6086","Lymph_node_healthy"="#20854E","Brain_healthy"="#4A7985")
tissue_type_color <- c("Bone_marrow"="#CE1261","Kidney"="#5E4FA2","Lymph_node"="#8CA77B","Brain"="#00441B","Lung"="#ed3d30",
                       "Colon"="#DEDC00","Liver"="#B3DE69","Peritoneal"="#999999","Stomach"="#725e82","Ovary"="#8c4356",
                       "Head_and_neck"="#db5a6b","Pancreas"="#21a675","Vaginal"="#d9b611","Appendix"="#e29c45",
                       "Bone"="#4c8dae","Nasopharynx"="#003371","Salivary_gland"="#a1afc9","Testis"="#88ada6","Cervix"="#c93756",
                       "Skin"="#204B75","Breast"="#588257","Thyroid"="#B6DB7B","Esophagus"="#4285BF","Larynx"="#1CB232")
tumor_status_color <- c("Tumor" = "#5BC0EB","Metastasis"="#8DD3C7","NAT" = "#9BC53D","Healthy" = "#C3423F")
study_database_color <- c("GSE202813"="#55d151","GSE237425"="#0e09b7","GSE171306"="#f69a95","GSE152938"="#aa329a","GSE210038"="#1d96a9",
                          "GSE277742"="#EDE816","GSE143423"="#48dc88","GSE234832"="#b15b0a","GSE256490"="#ac66d2","GSE164789"="#dc1c6f",
                          "HRA006468"="#CB00FC","GSE178318"="#af3292","GSE183916"="#5da3ac","GSE231559"="#c6ab48","GSE261388"="#0c9986",
                          "GSE245552"="#2f52d1","GSE271913"="#aa152e","GSE201347"="#729e0c","GSE163558"="#be5793","GSE246662"="#3cc0c0",
                          "GSE181919"="#9ebe71","GSE173468"="#9a221c","GSE188737"="#FD9A97","GSE256136"="#0c9c54","GSE162708"="#C0AFE0",
                          "GSE156405"="#D283FA","GSE197177"="#A30097","GSE281288"="#FDC44D","GSE205013"="#737D58","GSE212966"="#9e2d71",
                          "GSE244031"="#5A2AB4","GSE270231"="#0D6094","GSE237070"="#9F4700","GSE249361"="#16A795","GSE280127"="#BF658B",
                          "GSE217084"="#59d185","GSE197778"="#FED8A3","GSE228501"="#7E9800","GSE208653"="#F82263","GSE174401"="#A35FFD",
                          "GSE169147"="#2285FE","GSE280584"="#F65FFE","GSE200218"="#75ADCD","GSE222446"="#959099","GSE225600"="#FB5A8A",
                          "GSE190870"="#BD856E","GSE184362"="#A61C79","GSE232237"="#FC7216","GSE199619"="#76AE86","GSE196756"="#F87AA3",
                          "GSE188900"="#E3DDEF","GSE206332"="#CEE769","GSE252490"="#E7ABFD","GSE243977"="#53475A","GSE185477"="#847A74",
                          "GSE243639"="#6FCDDC","GSE203067"="#1CB232")


# data dir ----------------------------------------------------------------
sample_dir <- " "
input_dir <- " "
output_dir <- " "

# data read ---------------------------------------------------------------
## read
all_meta <- read_csv(str_glue("{input_dir}/all_cell_combine_meta.csv"))
group_data <- read.xlsx(str_glue("{sample_dir} "),sheetName = "scRNA_seq")
colnames(group_data)[2] <- "sample_ID"

## data filter
all_meta_filter <- all_meta %>% 
  filter(tumor_status %in% c("Metastasis","Healthy"))
  # filter(metastasis_tissue %in% c("Lymph_node","Brain","Liver","None"))


# data stat ---------------------------------------------------------------
## remove no hub cell type
all_meta_stat_filter <- all_meta_filter %>% 
  filter(final_cell_type != "SMC&Pericyte") %>% 
  filter(final_cell_type != "Bela_cell") %>%
  filter(final_cell_type != "Acinar_cell") %>%
  # filter(final_cell_type != "Melanoma_cell") %>%
  filter(final_cell_type != "Osteoblastic_cell")
  # filter(final_cell_type != "Melanoma_tumor") %>%
  # filter(final_cell_type != "Osteoblastic_tumor") %>%
  # filter(final_cell_type != "Epithelia_tumor")
  # filter(final_cell_type != "Neuron") %>%
  # filter(final_cell_type != "Glial_cell") %>%
  # filter(final_cell_type != "Hepatocyte")

## cell type confirm
all_meta_stat_filter <- all_meta_stat_filter %>%
  mutate(most_cell_type = case_when(
    grepl("tumor", final_cell_type) ~ "Tumor_cell",
    grepl("^B_", final_cell_type) ~ "B_cell",
    grepl("^EC", final_cell_type) ~ "Endothelia",
    grepl("^Plasma_", final_cell_type) ~ "Plasma_cell",
    grepl("^Fib_", final_cell_type) ~ "Fibroblast",
    
    grepl("^Melanoma_", final_cell_type) ~ "Melanoma_cell",
    grepl("^Neuron", final_cell_type) ~ "Neuron",
    grepl("^Glial_", final_cell_type) ~ "Glial_cell",
    grepl("^Hepatocyte", final_cell_type) ~ "Hepatocyte",
    
    grepl("^Pericyte|^SMC", final_cell_type) ~ "SMC&Pericyte",
    grepl("^T_CD8", final_cell_type) ~ "T_CD8_cell",
    grepl("^T_CD4", final_cell_type) ~ "T_CD4_cell",
    grepl("^NK_", final_cell_type) ~ "NK_cell",
    grepl("^Mono_", final_cell_type) ~ "Monocyte",
    grepl("^Mph_", final_cell_type) ~ "Macrophage",
    grepl("^DC_", final_cell_type) ~ "Dendritic_cell",
    grepl("^Neu_", final_cell_type) ~ "Neutrophils",
    grepl("^Mast|^Myeloid", final_cell_type) ~ "Myeloid_Other",
    TRUE ~ "Other"
  ))

## prop cal
meta_stat_prop <- all_meta_stat_filter %>% 
  filter(most_cell_type != "Other") %>% 
  dplyr::group_by(sample_ID,most_cell_type) %>% 
  dplyr::summarise(count = n(), .groups = "drop") %>% 
  dplyr::group_by(sample_ID) %>% 
  dplyr::mutate(percentage = count / sum(count)*100) %>% 
  ungroup

## result transform
meta_stat_prop_transform <- meta_stat_prop[,c(1,2,4)] %>%
  pivot_wider(names_from = most_cell_type,
              values_from = percentage,
              values_fill = 0) %>% 
  column_to_rownames(var = "sample_ID")



# PCA analysis ------------------------------------------------------------
## data PCA
PCA_data <- as.matrix(meta_stat_prop_transform)
PCA_result <- prcomp(PCA_data,center = T,scale. = T)

## result
PCA_result_df <- data.frame(
  sample_ID = rownames(PCA_data),
  PC1 = PCA_result$x[,1],
  PC2 = PCA_result$x[,2]
)
PCA_result_df <- left_join(PCA_result_df,group_data,by = "sample_ID")
PCA_result_df$group <- str_c(PCA_result_df$tissue_type,PCA_result_df$tumor_status,sep = "_")



# result plot -------------------------------------------------------------
## tumor status
ggplot(PCA_result_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, aes(color = tumor_status)) + 
  scale_color_manual(values = tumor_status_color) +
  nrc_theme + NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_PCA_tumor_status"), width = 150, height = 150)

## metastasis tissue type
ggplot(PCA_result_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, aes(color = group)) + 
  scale_color_manual(values = metastasis_color) +
  nrc_theme + NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_PCA_metastasis_status"), width = 150, height = 150)

## tumor code
ggplot(PCA_result_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, aes(color = tumor_code)) + 
  scale_color_manual(values = tumor_type_color) +
  nrc_theme + NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_PCA_tumor_code"), width = 150, height = 150)

## tissue type
ggplot(PCA_result_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, aes(color = tissue_type)) + 
  scale_color_manual(values = tissue_type_color) +
  nrc_theme + NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_PCA_tissue_type"), width = 150, height = 150)








