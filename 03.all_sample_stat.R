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
# .libPaths(c(" "))
# .libPaths(c(" "))

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
library("ggalluvial")
library("networkD3")
library("ggsankey")
library("ComplexHeatmap")
(.packages())
setwd(' ')
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
input_dir <- " "
output_dir <- " "

# data color --------------------------------------------------------------
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
bulk_tumor_code_color <- c("UVM"="#55d151","UCEC"="#0e09b7","UCS"="#f69a95","THCA"="#F9AC93","THYM"="#aa329a","TGCT"="#7876B1",
                           "STAD"="#718DAE","SKCM"="#FFDC91","SARC"="#E18727","READ"="#1d96a9","PRAD"="#EDE816","PCPG"="#48dc88",
                           "PAAD"="#b15b0a","OV"="#ac66d2","MESO"="#dc1c6f","LUSC"="#9FAFA3","LUAD"="#958056","LIHC"="#CB00FC",
                           "KIRP"="#af3292","KIRC"="#BC3C29","KICH"="#af3292","HNSC"="#7D4E57","GBM"="#5da3ac","ESCA"="#c6ab48",
                           "DLBC"="#0c9986","COAD"="#2f52d1","CHOL"="#aa152e","CESC"="#A08634","BRCA"="#3D806F","LGG"="#729e0c",
                           "BLCA"="#be5793","ACC"="#5A7B8F","LAML"="#3cc0c0")


# data read ---------------------------------------------------------------
scRNA_sample <- read_csv(str_glue("{input_dir}/scRNA_data_sample_filter.csv"))
stRNA_sample <- read.xlsx(str_glue("{input_dir} "),sheetName = "stRNA_seq")
bulkRNA_sample <- read_csv(" ")

# scRNA sample stat -------------------------------------------------------
## stat
tumor_code_stat <- as.data.frame(table(scRNA_sample$tumor_code))
colnames(tumor_code_stat) <- c("tumor_code","num")
tissue_type_stat <- as.data.frame(table(scRNA_sample$tissue_type))
colnames(tissue_type_stat) <- c("tissue_type","num")
tumor_status_stat <- as.data.frame(table(scRNA_sample$tumor_status))
colnames(tumor_status_stat) <- c("tumor_status","num")
study_ID_stat <- as.data.frame(table(scRNA_sample$study_ID))
colnames(study_ID_stat) <- c("study_ID","num")
scRNA_sample_metastasis <- scRNA_sample[scRNA_sample$tumor_status == "Metastasis",]
tissue_metastasis_stat <- as.data.frame(table(scRNA_sample_metastasis$tissue_type))
colnames(tissue_metastasis_stat) <- c("tissue_type","num")

## data plot
ggplot(tumor_code_stat, aes(reorder(tumor_code, - num), num, fill = tumor_code)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = num), vjust = -0.5, size = 5) +
  theme_bw() +
  labs(x = "", y = "Number of sample") +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = tumor_type_color) +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_tumor_type_stat"), width = 400, height = 200)

ggplot(tissue_type_stat, aes(reorder(tissue_type, - num), num, fill = tissue_type)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = num), vjust = -0.5, size = 5) +
  theme_bw() +
  labs(x = "", y = "Number of sample") +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = tissue_type_color) +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_tissue_type_stat"), width = 400, height = 200)

ggplot(study_ID_stat, aes(reorder(study_ID, - num), num, fill = study_ID)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = num), vjust = -0.5, size = 5) +
  theme_bw() +
  labs(x = "", y = "Number of sample") +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = study_database_color) +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_study_type_stat"), width = 400, height = 200)

## primary to metastasis plot
ggsankey_plot <- scRNA_sample %>%
  filter(tumor_status == "Metastasis") %>% 
  make_long(primary_tissue, metastasis_tissue)
ggplot(ggsankey_plot, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_manual(values = tissue_type_color) +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",plot.title = element_text(hjust = .5))
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_ggsankey_plot"), width = 200, height = 200)

## pheatmap plot
scRNA_pheatmap <- scRNA_sample %>% 
  arrange(study_ID,tumor_code,tissue_type,tumor_status)
  # arrange(tumor_status,tissue_type,tumor_code,study_ID)

ha <- HeatmapAnnotation(
  study_ID = scRNA_pheatmap$study_ID,
  tumor_code = scRNA_pheatmap$tumor_code,
  tissue_type = scRNA_pheatmap$tissue_type,
  
  tumor_status = scRNA_pheatmap$tumor_status,
  
  col = list(
    study_ID = study_database_color,tumor_status = tumor_status_color,tissue_type = tissue_type_color,tumor_code = tumor_type_color
  )
)
mat <- matrix(
  0,nrow = 1,ncol = nrow(scRNA_pheatmap)
)
colnames(mat) <- rownames(scRNA_pheatmap)

pdf(str_glue("{output_dir}/sample_heatmap.pdf"),width = 13,height = 9)
Heatmap(
  mat,
  name = "1",
  top_annotation = ha,
  show_row_names = F,
  show_column_names = F,
  cluster_rows = F,
  cluster_columns = F
  
  
)
dev.off()

## data save
if(T){
  write_csv(tumor_code_stat,str_glue("{output_dir}/scRNA_tumor_type_stat.csv"))
  write_csv(tissue_type_stat,str_glue("{output_dir}/scRNA_tissue_type_stat.csv"))
  write_csv(tumor_status_stat,str_glue("{output_dir}/scRNA_tumor_status_stat.csv"))
  write_csv(study_ID_stat,str_glue("{output_dir}/scRNA_study_ID_stat.csv"))
  write_csv(tissue_metastasis_stat,str_glue("{output_dir}/scRNA_tissue_metastasis_stat.csv"))
}


# stRNA sample stat -------------------------------------------------------
## stat
sttumor_code_stat <- as.data.frame(table(stRNA_sample$tumor_code))
colnames(sttumor_code_stat) <- c("tumor_code","num")
sttumor_status_stat <- as.data.frame(table(stRNA_sample$tumor_status))
colnames(sttumor_status_stat) <- c("tumor_status","num")

## data plot
ggplot(sttumor_code_stat, aes(reorder(tumor_code, - num), num, fill = tumor_code)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = num), vjust = -0.5, size = 5) +
  theme_bw() +
  labs(x = "", y = "Number of sample") +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = tumor_type_color) +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/stRNA_tumor_type_stat"), width = 200, height = 200)

## data save
write_csv(sttumor_code_stat,str_glue("{output_dir}/stRNA_tumor_type_stat.csv"))
write_csv(sttumor_status_stat,str_glue("{output_dir}/stRNA_tumor_status_stat.csv"))

# bulkRNA sample stat -----------------------------------------------------
bulkRNA_sample <- bulkRNA_sample[!duplicated(bulkRNA_sample$sample),]
bulkRNA_sample <- bulkRNA_sample[bulkRNA_sample$tumor_code != "other",]

bulkRNA_tumor_code_stat <- as.data.frame(table(bulkRNA_sample$tumor_code))
colnames(bulkRNA_tumor_code_stat) <- c("tumor_code","num")
bulkRNA_tissue_type_stat <- as.data.frame(table(bulkRNA_sample$tissue_type))
colnames(bulkRNA_tissue_type_stat) <- c("tissue_type","num")

## data plot
ggplot(bulkRNA_tumor_code_stat, aes(reorder(tumor_code, - num), num, fill = tumor_code)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = num), vjust = -0.5, size = 5) +
  theme_bw() +
  labs(x = "", y = "Number of sample") +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = bulk_tumor_code_color) +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/bulkRNA_tumor_type_stat"), width = 450, height = 200)

## data save
write_csv(bulkRNA_tumor_code_stat,str_glue("{output_dir}/bulkRNA_tumor_type_stat.csv"))

