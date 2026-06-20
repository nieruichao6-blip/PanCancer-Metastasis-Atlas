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
library("SCP")
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
CIN70_geneset_dir <- " "
NMF_MP_dir <- " "
input_dir <- " "
output_dir <- " "
tumor_code_list <- c("ACC","BRCA","CESC","CRC","ESCC","HNSC","KIRC","NPC","PAAD","PNET","GC",
                     "THCA","LUAD","LSCC")

# data color --------------------------------------------------------------
tumor_type_color <- c("ACC"="#5A7B8F","ANS"="#EE4C97","BRCA"="#3D806F","CESC"="#A08634","CRC"="#F37C95","ESCC"="#608541",
                      "HNSC"="#7D4E57","KIRC"="#BC3C29","LUAD"="#958056","LUSC"="#9FAFA3","NPC"="#6F99AD","PDAC"="#0072B5",
                      "PNET"="#CFC59A","SARC"="#E18727","SKCM"="#FFDC91","GC"="#718DAE","TGCT"="#7876B1","THCA"="#F9AC93",
                      "NSCLC"="#E3BC06","LSCC"="#2285FE","Liver_healthy"="#3E6086","Lymph_node_healthy"="#20854E","Brain_healthy"="#4A7985")
tissue_type_color <- c("Bone_marrow"="#CE1261","Kidney"="#5E4FA2","Lymph_node"="#8CA77B","Brain"="#00441B","Lung"="#ed3d30",
                       "Colon"="#DEDC00","Liver"="#B3DE69","Peritoneal"="#999999","Stomach"="#725e82","Ovary"="#8c4356",
                       "Head_and_neck"="#db5a6b","Pancreas"="#21a675","Vaginal"="#d9b611","Appendix"="#e29c45",
                       "Bone"="#4c8dae","Nasopharynx"="#003371","Salivary_gland"="#a1afc9","Testis"="#88ada6","Cervix"="#c93756",
                       "Skin"="#204B75","Breast"="#588257","Thyroid"="#B6DB7B","Esophagus"="#4285BF","larynx"="#1CB232")
tumor_status_color <- c("Tumor" = "#5BC0EB","Metastasis"="#8DD3C7","NAT" = "#9BC53D","Healthy" = "#C3423F")



# data read ---------------------------------------------------------------
NMF_MP_result <- read_csv(NMF_MP_dir)
CIN70_geneset <- read.xlsx(CIN70_geneset_dir,sheetName = "CIN70")
CIN70_geneset <- list(c(CIN70_geneset$CIN70))
names(CIN70_geneset) <- "CIN70"

scRNA_tumor <- read_h5(file = str_glue("{input_dir}/scRNA_bbknn.h5"),
                       assay.name = 'RNA', 
                       target.object = 'seurat')

# addmodulescore ---------------------------------------------------------
scRNA_tumor <- AddModuleScore(scRNA_tumor,
                       features = CIN70_geneset,
                       ctrl = 5,
                       name = "CIN_addmodulescore")
scaled <- scale(scRNA_tumor$CIN_addmodulescore1)
scRNA_tumor$score_normalized <- (scaled - min(scaled)) / (max(scaled) - min(scaled))

# AUCell ------------------------------------------------------------------
cells_rankings <- AUCell_buildRankings(scRNA_tumor@assays$RNA@data,splitByBlocks=TRUE)
cells_AUC <- AUCell_calcAUC(CIN70_geneset, cells_rankings, 
                            aucMaxRank=nrow(cells_rankings)*0.1)
auc_score <- as.data.frame(getAUC(cells_AUC)["CIN70",])
colnames(auc_score) <- "CIN_aucell"
scRNA_tumor <- AddMetaData(scRNA_tumor,auc_score)

# addmodulescore result plot -------------------------------------------------------------
addmodulescore_meta <- scRNA_tumor@meta.data
# addmodulescore_meta_all <- addmodulescore_meta %>%
#   dplyr::group_by(sample_ID,tumor_status) %>%
#   dplyr::summarise(mean_sample_score = mean(score_normalized,na.rm = T),.groups = 'drop') %>% 
#   ungroup()
## diff tissye
# addmodulescore_meta <- read_csv(str_glue("{output_dir}/CIN70_addmodulescore_meta.csv"))
ggviolin(addmodulescore_meta, 
         x="tumor_status", y="score_normalized", 
         width = 0.8,color = "black",
         fill="tumor_status",
         xlab = F,
         add = "boxplot",
         add.params = list(outlier.shape = NA),
         bxp.errorbar=T, 
         bxp.errorbar.width=0.05,
         size=0.5,
         palette = tumor_status_color, 
         legend = NULL) + nrc_theme +
    stat_compare_means(comparisons = list(c("Metastasis","Tumor"))) + NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/addmodulescore_diff_tissue"), width = 100, height = 180)

## diff tumor code
addmodulescore_meta_hub <- addmodulescore_meta %>% 
  filter(tumor_code %in% tumor_code_list)
ggviolin(addmodulescore_meta_hub,
         x="tumor_code", y="score_normalized", 
         width = 0.8,color = "black",
         fill="tumor_status",
         xlab = F,
         add = "boxplot",
         add.params = list(outlier.shape = NA),
         bxp.errorbar=T, 
         bxp.errorbar.width=0.05,
         size=0.5,
         palette = tumor_status_color, 
         legend = NULL) +
    stat_compare_means(aes(group=tumor_status), method="wilcox.test",
        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
        symbols = c("***", "**", "*", " ")), label = "p.signif") + nrc_theme +
    NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/addmodulescore_tumor_code"), width = 300, height = 180)

## diff metastasis tissue
addmodulescore_meta_metastasis <- addmodulescore_meta %>% 
  filter(tumor_status == "Metastasis")
addmodulescore_meta_metastasis <- addmodulescore_meta_metastasis %>% 
  mutate(metastasis_tissue = fct_reorder(metastasis_tissue,score_normalized,.desc = T,.fun = median))
ggviolin(addmodulescore_meta_metastasis,
         x="metastasis_tissue", y="score_normalized", 
         width = 0.8,color = "black",
         fill="metastasis_tissue",
         xlab = F,
         add = "boxplot",
         add.params = list(outlier.shape = NA),
         bxp.errorbar=T,
         bxp.errorbar.width=0.05,
         size=0.5,
         palette = tissue_type_color, 
         legend = NULL) +
    # stat_compare_means(comparisons = list( c("Lymph_node", "Brain"), c("Lymph_node", "Liver"), c("Brain", "Liver"))) + 
    nrc_theme +
    NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/addmodulescore_metastasis_tissue"), width = 200, height = 140)

## diff primary tissue
addmodulescore_meta_primary <- addmodulescore_meta %>% 
  filter(tumor_status == "Tumor")
addmodulescore_meta_primary <- addmodulescore_meta_primary %>% 
  mutate(tumor_code = fct_reorder(tumor_code,score_normalized,.desc = T,.fun = median))
ggviolin(addmodulescore_meta_primary,
         x="tumor_code", y="score_normalized", 
         width = 0.8,color = "black",
         fill="tumor_code",
         xlab = F,
         add = "boxplot",
         add.params = list(outlier.shape = NA),
         bxp.errorbar=T, 
         bxp.errorbar.width=0.05,
         size=0.5,
         palette = tumor_type_color, 
         legend = NULL) +
  # stat_compare_means(comparisons = list( c("Lymph_node", "Brain"), c("Lymph_node", "Liver"), c("Brain", "Liver"))) + 
  nrc_theme +
  NoLegend()

# AUCell result plot -------------------------------------------------------------
aucell_meta <- scRNA_tumor@meta.data
## diff tissye
ggviolin(aucell_meta, 
         x="tumor_status", y="CIN_aucell", 
         width = 0.8,color = "black",
         fill="tumor_status",
         xlab = F,
         add = "boxplot",
         add.params = list(outlier.shape = NA),
         bxp.errorbar=T, 
         bxp.errorbar.width=0.05,
         size=0.5,
         palette = tumor_status_color, 
         legend = NULL) +
    stat_compare_means(comparisons = list(c("Metastasis","Tumor"))) + NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/aucell_diff_tissue"), width = 80, height = 180)

## diff tumor code
aucell_meta <- aucell_meta %>% 
  filter(tumor_code %in% tumor_code_list)
ggviolin(aucell_meta,
         x="tumor_code", y="CIN_aucell", 
         width = 0.8,color = "black",
         fill="tumor_status",
         xlab = F,
         add = "boxplot",
         add.params = list(outlier.shape = NA),
         bxp.errorbar=T, 
         bxp.errorbar.width=0.05,
         size=0.5,
         palette = tumor_status_color, 
         legend = NULL) +
    stat_compare_means(aes(group=tumor_status), method="wilcox.test",
        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
        symbols = c("***", "**", "*", " ")), label = "p.signif") + nrc_theme +
    NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/aucell_tumor_code"), width = 300, height = 180)



# diff metastasis tissue plot ---------------------------------------------
for(tissue in c("Brain","Liver","Lymph_node","Other")){
  # tissue <- "Other"
  
  if(tissue == "Other"){
    current_meta_metastasis <- addmodulescore_meta_hub %>% 
      filter(metastasis_tissue %in% c("Lung","Peritoneal","Ovary","Bone_marrow","Vaginal"))
  }else{
    current_meta_metastasis <- addmodulescore_meta_hub %>% 
      filter(metastasis_tissue == tissue)
  }
  hub_tumor <- as.character(unique(current_meta_metastasis$tumor_code))
  current_meta_primary <- addmodulescore_meta_hub %>% 
    filter(tumor_code %in% hub_tumor) %>% 
    filter(tumor_status == "Tumor")
  current_meta <- rbind(current_meta_primary,current_meta_metastasis)
  
  ggviolin(current_meta,
           x="tumor_code", y="score_normalized", 
           width = 0.8,color = "black",
           fill="tumor_status",
           xlab = F,
           add = "boxplot",
           add.params = list(outlier.shape = NA),
           bxp.errorbar=T, 
           bxp.errorbar.width=0.05,
           size=0.5,
           palette = tumor_status_color, 
           legend = NULL) +
    stat_compare_means(aes(group=tumor_status), method="wilcox.test",
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", " ")), label = "p.signif") + nrc_theme +
    NoLegend()
  function_plot(filename_prefix = str_glue("{output_dir}/addmodulescore_diff_metastasis_tissue_{tissue}"), width = length(hub_tumor) * 30, height = 150)
  
}



# CIN70 NMF cor -----------------------------------------------------------
CIN_NMF_cor_result <- data.frame(
  final_MP = as.character(),
  cor = as.numeric(),
  pvalue = as.numeric()
)

## data preapre
MP_filter <- NMF_MP_result %>% 
  dplyr::select(c(18,24:51))
MP_name <- colnames(MP_filter)[-1]
CIN_socre <- addmodulescore_meta %>% 
  dplyr::select(cell_id,score_normalized)
CIN_MP_data <- left_join(CIN_socre,MP_filter,by = "cell_id")

## cor process
for(MP in MP_name){
  # MP <- "pearson"
  cor_result <- cor.test(CIN_MP_data[["score_normalized"]],CIN_MP_data[[MP]],method = "pearson")
  
  ## data combine
  current_cor <- data.frame(
    final_MP = MP,
    cor = cor_result$estimate,
    pvalue = cor_result$p.value
  )
  CIN_NMF_cor_result <- rbind(CIN_NMF_cor_result,current_cor)
}

# data save ---------------------------------------------------------------
write_csv(addmodulescore_meta,str_glue("{output_dir}/CIN70_addmodulescore_meta.csv"))





