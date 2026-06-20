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
hub_metastasis_geneset_dir <- " "
NMF_MP_dir <- " "
input_dir <- " "
output_dir <- " "
tumor_code_list <- c("BRCA","CESC","CRC","ESCC","HNSC","KIRC","NPC","PAAD","PNET","GC",
                     "THCA","NSCLC","LSCC")

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
new_status_color <- c("NAT"="#9BC53D","Primary_tumor"="#5BC0EB",
                      "Brain_metastasis"="#C69D3B","Liver_metastasis"="#EFC971","Lymph_node_metastasis"="#48A83E","Other_metastasis"="#7ABC68",
                      "Brain_Healthy"="#D86127","Liver_Healthy"="#F47D6E","Lymph_node_Healthy"="#F49897")
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

# data read ---------------------------------------------------------------
## scRNA read
scRNA_tumor <- read_h5(file = str_glue("{input_dir}/scRNA_cell_type_cluster_only.h5"),
                       assay.name = 'RNA', 
                       target.object = 'seurat')

## metastasis hub gene read
metastasis_gene <- read_csv(str_glue("{hub_metastasis_geneset_dir}/high_per_DEG_all_metastasis.csv"))
colnames(metastasis_gene) <- "gene_name"
metastasis_gene_list <- list(metastasis_hub_gene = metastasis_gene$gene_name)


# HLA gene exp ------------------------------------------------------------
## sample desc
sample_data <- scRNA_tumor@meta.data[,c("sample_ID","tumor_code","tissue_type","metastasis_tissue","tumor_type","tumor_status","study_ID")]
sample_data <- sample_data[!duplicated(sample_data$sample_ID),]
sample_data <- sample_data %>% 
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
    
  )) %>% 
  arrange(metastasis_tissue_transform)

## gene expression
Idents(scRNA_tumor) <- scRNA_tumor$sample_ID
all_gene_expression <- AverageExpression(scRNA_tumor,group.by = "sample_ID")
HLA_gene <- c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F",
              "HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB")
# HLA_gene_expression <- all_gene_expression$RNA[HLA_gene,sample_data$sample_ID]
sample_data$mean_expr <- colMeans(HLA_gene_expression)
sample_data <- sample_data[
  order(
    sample_data$metastasis_tissue_transform,
    -sample_data$mean_expr
  ),
]
HLA_gene_expression <- all_gene_expression$RNA[HLA_gene,sample_data$sample_ID]
cluster_ano <- sample_data[,c("tumor_status","metastasis_tissue_transform","tumor_code","tissue_type","study_ID")]
breaks <- c(seq(-2, 0, length = 50), seq(0.01,2, length = 50))
col.heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))
cols <- col.heat(100)
col1 <- colorRamp3(breaks, cols)
pdf(str_glue("{output_dir}/HLA_gene_expression.pdf"),width = 13,height = 4)
ComplexHeatmap::pheatmap(
  as.matrix(HLA_gene_expression),
  annotation_col = cluster_ano,
  show_colnames = FALSE,
  show_rownames = TRUE,
  treeheight_row = 20,
  treeheight_col = 50,
  clustering_method = 'ward.D2',
  # scale = "row",
  cluster_rows = F,
  cluster_cols = F,
  cutree_cols = 20,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  gaps_col = c(50,95,136,159),
  gaps_row = c(5),
  scale = "row",
  annotation_colors = list(tumor_code = c(tumor_type_color),metastasis_tissue_transform = c(new_status_color),
                           tissue_type = c(tissue_type_color),tumor_status = c(tumor_status_color),study_ID = c(study_database_color)),
  legend = TRUE,
  color = cols,
  breaks = seq(-1.5, 2, length.out = 101)
)
dev.off()

# addmodulescore ---------------------------------------------------------
scRNA_tumor <- AddModuleScore(scRNA_tumor,
                       features = metastasis_gene_list,
                       ctrl = 5,
                       name = "metastasis_hub_gene")
scaled <- scale(scRNA_tumor$metastasis_hub_gene1)
scRNA_tumor$score_normalized <- (scaled - min(scaled)) / (max(scaled) - min(scaled))

# addmodulescore result plot -------------------------------------------------------------
addmodulescore_meta <- scRNA_tumor@meta.data
addmodulescore_meta_all <- addmodulescore_meta %>%
  dplyr::group_by(sample_ID,tumor_status) %>%
  dplyr::summarise(mean_sample_score = mean(score_normalized,na.rm = T),.groups = 'drop') %>% 
  ungroup()
hub_meta <- addmodulescore_meta[,c("sample_ID","tumor_code")]
hub_meta <- hub_meta[!duplicated(hub_meta$sample_ID),]
addmodulescore_meta_all <- left_join(addmodulescore_meta_all,hub_meta,by = "sample_ID")

## diff tissye
ggviolin(addmodulescore_meta_all, 
         x="tumor_status", y="mean_sample_score", 
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
addmodulescore_meta[addmodulescore_meta$tumor_code %in% c("NSCLC","LUAD","LUSC"),]$tumor_code <- "NSCLC"
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
function_plot(filename_prefix = str_glue("{output_dir}/addmodulescore_primary_tissue"), width = 200, height = 140)

# diff metastasis tissue plot ---------------------------------------------
for(tissue in c("Brain","Liver","Lymph_node","Other")){
  tissue <- "Lymph_node"
  
  if(tissue == "Other"){
    current_meta_metastasis <- addmodulescore_meta %>% 
      filter(metastasis_tissue %in% c("Lung","Peritoneal","Ovary","Bone_marrow","Vaginal"))
  }else{
    current_meta_metastasis <- addmodulescore_meta %>% 
      filter(metastasis_tissue == tissue)
  }
  hub_tumor <- as.character(unique(current_meta_metastasis$tumor_code))
  current_meta_primary <- addmodulescore_meta %>% 
    filter(tumor_code %in% hub_tumor) %>% 
    filter(tumor_status == "Tumor")
  current_meta <- rbind(current_meta_primary,current_meta_metastasis)
  current_meta[current_meta$tumor_code %in% c("NSCLC","LUAD","LUSC"),]$tumor_code <- "NSCLC"
  current_meta <- current_meta[current_meta$tumor_code != "TGCT",]
  
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


