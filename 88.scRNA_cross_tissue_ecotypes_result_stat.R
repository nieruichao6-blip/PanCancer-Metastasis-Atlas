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

# data color ----------------------------------------------------------------
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
ecotypes_result_dir <- " "
output_dir <- " "


# data read ---------------------------------------------------------------
ecotypes_result <- read_csv(str_glue("{ecotypes_result_dir}/cluster_ano.csv"))



# status stat -------------------------------------------------------------
ET_status_all <- data.frame(
  Ecotype = as.character(),
  NAT = as.numeric(),
  Tumor = as.numeric(),
  Metastasis = as.numeric(),
  Healthy = as.numeric()
)

## stat
for(ET_type in unique(ecotypes_result$Ecotype)){
  # ET_type <- "ET_01"
  
  ### data filter
  current_result <- ecotypes_result %>% 
    filter(Ecotype == ET_type)
  
  ### num
  current_NAT <- nrow(current_result[current_result$tumor_status == "NAT",])/111
  current_primary <- nrow(current_result[current_result$tumor_status == "Tumor",])/213
  current_metastasis <- nrow(current_result[current_result$tumor_status == "Metastasis",])/197
  current_Healthy <- nrow(current_result[current_result$tumor_status == "Healthy",])/21
  
  ### data collect
  current_collect <- data.frame(
    Ecotype = ET_type,
    NAT = current_NAT,
    Tumor = current_primary,
    Metastasis = current_metastasis,
    Healthy = current_Healthy
  )
  ET_status_all <- rbind(ET_status_all,current_collect)
  
}

## data transform
write_csv(ET_status_all,str_glue("{output_dir}/ET_status_all.csv"))
ET_status_melt <- melt(ET_status_all)
colnames(ET_status_melt) <- c("Ecotype","tumor_status","percentage")
ET_status_melt[ET_status_melt$percentage == 0,]$percentage <- 0.001

ggplot(ET_status_melt, aes(Ecotype ,percentage, fill = tumor_status)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  theme_bw() +
  labs(x = "", y = "Percentage of Ecotype") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = tumor_status_color) +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_status_Ecotype_all"), width = 280, height = 150)


# tumor code stat ---------------------------------------------------------
ET_tumor_code_all <- data.frame(
  Ecotype = as.character(),
  num = as.numeric(),
  all_num = as.numeric(),
  percentage = as.numeric(),
  tumor_code = as.character()
)

## data filter
ecotypes_result_filter <- ecotypes_result %>% 
  filter(tumor_status != "Healthy")
ecotypes_result_filter[ecotypes_result_filter$tumor_code == "LUAD" | ecotypes_result_filter$tumor_code == "LUSC",]$tumor_code <- "NSCLC"
all_ecotype_num <- as.data.frame(table(ecotypes_result_filter$Ecotype))
colnames(all_ecotype_num) <- c("Ecotype","all_num")

## stat
for(tumor in unique(ecotypes_result_filter$tumor_code)){
  # tumor <- "NSCLC"
  
  ### data filter
  current_result <- ecotypes_result_filter %>% 
    filter(tumor_code == tumor)
  
  ### num
  current_collect <- as.data.frame(table(current_result$Ecotype))
  colnames(current_collect) <- c("Ecotype","num")
  current_collect <- left_join(current_collect,all_ecotype_num,by = "Ecotype")
  current_collect$percentage <- current_collect$num/current_collect$all_num
  current_collect$tumor_code <- tumor
  
  ### data plot
  ET_tumor_code_all <- rbind(ET_tumor_code_all,current_collect)
  
}
ET_tumor_code_all_transform <- ET_tumor_code_all[,c(1,4,5)] %>%
  pivot_wider(names_from = tumor_code,
              values_from = percentage,
              values_fill = 0)
write_csv(ET_tumor_code_all_transform,str_glue("{output_dir}/ET_tumor_code_all_transform.csv"))

## result transform
ET_tumor_code_all_transform[,-1] <- apply(ET_tumor_code_all_transform[,-1], 2, function(x){
  x / sum(x)
})

ET_tumor_code_all_transform <- ET_tumor_code_all_transform %>%
  pivot_longer(
    cols = -Ecotype,
    names_to = "tumor_code",
    values_to = "percentage"
  )

ggplot(ET_tumor_code_all_transform, aes(x = tumor_code, stratum = Ecotype, alluvium = Ecotype, y = percentage, fill = Ecotype)) +
  geom_stratum(width = 0.5, col=NA) +
  geom_alluvium(width = 0.4, alpha = 0.4) +
  scale_fill_manual(values = Ecotype_color) +
  labs(x = '', y = 'Proportion of ET')+
  theme_classic()+ nrc_theme + NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_tumor_code_Ecotype_all"), width = 280, height = 150)


# tissue type stat --------------------------------------------------------
ET_tissue_all <- data.frame(
  Ecotype = as.character(),
  num = as.numeric(),
  all_num = as.numeric(),
  percentage = as.numeric(),
  tissue_type = as.character()
)

## data filter
all_ecotype_num <- as.data.frame(table(ecotypes_result$Ecotype))
colnames(all_ecotype_num) <- c("Ecotype","all_num")

## stat
for(tissue in unique(ecotypes_result$tissue_type)){
  # tumor <- "NSCLC"
  
  ### data filter
  current_result <- ecotypes_result %>% 
    filter(tissue_type == tissue)
  
  ### num
  current_collect <- as.data.frame(table(current_result$Ecotype))
  colnames(current_collect) <- c("Ecotype","num")
  current_collect <- left_join(current_collect,all_ecotype_num,by = "Ecotype")
  current_collect$percentage <- current_collect$num/current_collect$all_num
  current_collect$tissue_type <- tissue
  
  ### data plot
  ET_tissue_all <- rbind(ET_tissue_all,current_collect)
  
}

ET_tissue_all_transform <- ET_tissue_all[,c(1,4,5)] %>%
  pivot_wider(names_from = tissue_type,
              values_from = percentage,
              values_fill = 0)
write_csv(ET_tissue_all_transform,str_glue("{output_dir}/ET_tissue_all_transform.csv"))

## result transform
ET_tissue_all_transform[,-1] <- apply(ET_tissue_all_transform[,-1], 2, function(x){
  x / sum(x)
})

ET_tissue_all_transform <- ET_tissue_all_transform %>%
  pivot_longer(
    cols = -Ecotype,
    names_to = "tissue_type",
    values_to = "percentage"
  )

ggplot(ET_tissue_all_transform, aes(x = tissue_type, stratum = Ecotype, alluvium = Ecotype, y = percentage, fill = Ecotype)) +
  geom_stratum(width = 0.5, col=NA) +
  geom_alluvium(width = 0.4, alpha = 0.4) +
  scale_fill_manual(values = Ecotype_color) +
  labs(x = '', y = 'Proportion of ET')+
  theme_classic()+ nrc_theme + NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_tissue_type_Ecotype_all"), width = 350, height = 150)


# metastasis stat ---------------------------------------------------------
ET_metastasis_tissue_all <- data.frame(
  Ecotype = as.character(),
  num = as.numeric(),
  all_num = as.numeric(),
  percentage = as.numeric(),
  tissue_type = as.character()
)

ecotypes_result_metastasis <- ecotypes_result %>% 
  filter(tumor_status == "Metastasis")

## data stat
all_meta_ET <- as.data.frame(table(ecotypes_result_metastasis$Ecotype))
colnames(all_meta_ET) <- c("Ecotype","num_all")

for(tissue in unique(ecotypes_result_metastasis$tissue_type)){
  # tissue <- "Brain"
  
  ### data filter
  current_data <- ecotypes_result_metastasis %>% 
    filter(tissue_type ==  tissue)
  
  ### num
  current_collect <- as.data.frame(table(current_data$Ecotype))
  colnames(current_collect) <- c("Ecotype","num")
  current_collect <- left_join(current_collect,all_meta_ET,by = "Ecotype")
  current_collect$percentage <- current_collect$num/current_collect$num_all
  current_collect$tissue_type <- tissue
  
  ### data plot
  ET_metastasis_tissue_all <- rbind(ET_metastasis_tissue_all,current_collect)
  
}

## result transform
ET_metastasis_tissue_all_transform <- ET_metastasis_tissue_all[,c(1,4,5)] %>%
  pivot_wider(names_from = tissue_type,
              values_from = percentage,
              values_fill = 0)
write_csv(ET_metastasis_tissue_all_transform,str_glue("{output_dir}/ET_metastasis_tissue_all_transform.csv"))
ET_metastasis_tissue_all_transform_melt <- melt(ET_metastasis_tissue_all_transform)
colnames(ET_metastasis_tissue_all_transform_melt) <- c("Ecotype","tissue","percentage")
ET_metastasis_tissue_all_transform_melt[ET_metastasis_tissue_all_transform_melt$percentage == 0,]$percentage <- 0.005

ggplot(ET_metastasis_tissue_all_transform_melt, aes(tissue ,percentage, fill = Ecotype)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  theme_bw() +
  labs(x = "", y = "Percentage of Ecotype") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = Ecotype_color) +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_metastasis_Ecotype_tissue"), width = 320, height = 150)

# primary stat ------------------------------------------------------------
ET_primary_tissue_all <- data.frame(
  Ecotype = as.character(),
  num = as.numeric(),
  all_num = as.numeric(),
  percentage = as.numeric(),
  tumor_code = as.character()
)

ecotypes_result_primary <- ecotypes_result %>% 
  filter(tumor_status == "Tumor")

## data stat
all_meta_ET <- as.data.frame(table(ecotypes_result_primary$Ecotype))
colnames(all_meta_ET) <- c("Ecotype","num_all")

for(tumor in unique(ecotypes_result_primary$tumor_code)){
  # tissue <- "Brain"
  
  ### data filter
  current_data <- ecotypes_result_primary %>% 
    filter(tumor_code ==  tumor)
  
  ### num
  current_collect <- as.data.frame(table(current_data$Ecotype))
  colnames(current_collect) <- c("Ecotype","num")
  current_collect <- left_join(current_collect,all_meta_ET,by = "Ecotype")
  current_collect$percentage <- current_collect$num/current_collect$num_all
  current_collect$tumor_code <- tumor
  
  ### data plot
  ET_primary_tissue_all <- rbind(ET_primary_tissue_all,current_collect)
  
}

## result transform
ET_primary_tissue_all_transform <- ET_primary_tissue_all[,c(1,4,5)] %>%
  pivot_wider(names_from = tumor_code,
              values_from = percentage,
              values_fill = 0)
write_csv(ET_primary_tissue_all_transform,str_glue("{output_dir}/ET_primary_tissue_all_transform.csv"))

ET_primary_tissue_all_transform_melt <- melt(ET_primary_tissue_all_transform)
colnames(ET_primary_tissue_all_transform_melt) <- c("Ecotype","tumor_code","percentage")
ET_primary_tissue_all_transform_melt[ET_primary_tissue_all_transform_melt$percentage == 0,]$percentage <- 0.005

ggplot(ET_primary_tissue_all_transform_melt, aes(tumor_code ,percentage, fill = Ecotype)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  theme_bw() +
  labs(x = "", y = "Percentage of Ecotype") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = Ecotype_color) +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_priamry_Ecotype_tissue"), width = 550, height = 150)


# ET stat -----------------------------------------------------------------
ET_all <- data.frame(
  tissue = as.character(),
  num = as.numeric(),
  all_num = as.numeric(),
  percentage = as.numeric(),
  Ecotype = as.character()
)

ecotypes_result_metastasis <- ecotypes_result %>% 
  filter(tumor_status == "Metastasis")

## data stat
all_meta_ET <- as.data.frame(table(ecotypes_result_metastasis$tissue_type))
colnames(all_meta_ET) <- c("tissue","num_all")

for(tumor in unique(ecotypes_result_metastasis$Ecotype)){
  # tumor <- "ET_01"
  
  ### data filter
  current_data <- ecotypes_result_metastasis %>% 
    filter(Ecotype ==  tumor)
  
  ### num
  current_collect <- as.data.frame(table(current_data$tissue_type))
  colnames(current_collect) <- c("tissue","num")
  current_collect <- left_join(current_collect,all_meta_ET,by = "tissue")
  current_collect$percentage <- current_collect$num/current_collect$num_all
  current_collect$Ecotype <- tumor
  
  ### data plot
  ET_all <- rbind(ET_all,current_collect)
  
}

## result transform
ET_all_transform <- ET_all[,c(1,4,5)] %>%
  pivot_wider(names_from = Ecotype,
              values_from = percentage,
              values_fill = 0)
# write_csv(ET_primary_tissue_all_transform,str_glue("{output_dir}/ET_primary_tissue_all_transform.csv"))

ET_all_transform_melt <- melt(ET_all_transform)
colnames(ET_all_transform_melt) <- c("tissue","Ecotype","percentage")
ET_all_transform_melt[ET_all_transform_melt$percentage == 0,]$percentage <- 0.005

ggplot(ET_all_transform_melt, aes(Ecotype ,percentage, fill = tissue)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  theme_bw() +
  labs(x = "", y = "Percentage of Ecotype") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = tissue_type_color) +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_tissue_metastasis_Ecotype_tissue"), width = 350, height = 150)


UNC5B_data <- readRDS(" ")
cell_marker <- c("CD3D","CD3E","CD2","NKG7","GNLY","KLRF1","CD79A","CD79B","MZB1","IGHA1","IGHG1","LYZ","AIF1","C1QC","CSF1R",
                 "IL1R2","CXCL8","TPSAB1","CPA3","DCN","COL1A1","RGS5","ACTA2","PECAM1","VWF","KRT8","KRT18","KRT19","UNC5B")
direct <- c("T_cell","NK_cell","B_cell","Plasma_cell","Myeloid_cell","Mast_cell","Fibroblast","SMC&Pericyte","Endothelia",
            "Epithelia","Epithelia_tumor")
UNC5B_data <- NormalizeData(UNC5B_data,normalization.method = "LogNormalize",scale.factor = 10000)
DotPlot(UNC5B_data, feature = cell_marker,group.by = "new_cell_type") + nrc_theme + ylab("") + xlab("") +
  scale_y_discrete(limits = direct) +
  # scale_size_continuous(range = c(0.1,8), breaks = c(0,20,40,60,80,100)) +
  scale_color_gradient2(high = "#D7301F", mid = "#A1D99B", low = "#6BAED6")

UNC5B_data_primary <- subset(UNC5B_data,tumor_status %in% c("Tumor"))
UNC5B_data_primary <- NormalizeData(UNC5B_data_primary,normalization.method = "LogNormalize",scale.factor = 10000)
DotPlot(UNC5B_data_primary, feature = cell_marker,group.by = "new_cell_type") + nrc_theme + ylab("") + xlab("") +
  scale_y_discrete(limits = direct) +
  scale_size_continuous(range = c(0.1,8), breaks = c(0,20,40,60,80,100)) +
  scale_color_gradient2(high = "#D7301F", mid = "#A1D99B", low = "#6BAED6")

