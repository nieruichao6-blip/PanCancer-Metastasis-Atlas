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
library("ggrastr")
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
# MP1 = '#F6CF71', MP2 = '#3969AC', MP3 = '#80BA5A', MP4 = '#F2B701', MP5 = '#11A579', MP6 = '#CF1C90', MP7 = '#66C5CC', MP8 = '#f97b72', MP9 = '#ed5887'
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
                     "Myeloid_Other"="#CAB2D6",
                     "Mast_cell"="#6A3D9A",
                     "Neutrophils"="#FFFF99",
                     "Macrophage"="#9ec9e1",
                     "Monocyte"="#E66F00",
                     "Dendritic_cell"="#1d92c0"
)


# data dir ----------------------------------------------------------------
sample_dir <- " "
input_dir <- " "
output_dir <- " "


# data read ---------------------------------------------------------------
scRNA_cross_tissue_proportion <- read_csv(str_glue("{input_dir}/cell_type_sample_proportion_hub.csv"))
group_data <- read.xlsx(str_glue("{sample_dir} "),sheetName = "scRNA_seq")

# data transform ----------------------------------------------------------
scRNA_cross_tissue_proportion_filter <- scRNA_cross_tissue_proportion %>% 
  left_join(.,group_data,by = "sample_ID") %>% 
  filter(tumor_status %in% c("NAT","Tumor","Metastasis","Healthy"))
  # filter(tumor_status %in% c("Metastasis"))

scRNA_cross_tissue_proportion_transform <- scRNA_cross_tissue_proportion_filter[,c(1,2,5)] %>%
  pivot_wider(names_from = sample_ID,
              values_from = clust_prop,
              values_fill = 0)
  # arrange(desc(final_cell_type)) %>%
  # column_to_rownames(var = "final_cell_type")


## reorder data
order_prefix <- c(
  "Tumor_cell","Epithelia_","T_CD8","T_CD4","T_ST","NK_","B_","Plasma_","Fib_","EC_","SMC","Pericyte_",
  "Mph_","Mono_","DC_","Myeloid","Neu_","Mast"
  # "Glial_",
  # "Hepatocyte"
)

scRNA_cross_tissue_proportion_transform <- scRNA_cross_tissue_proportion_transform %>%
  mutate(main_type = case_when(
    str_starts(final_cell_type, "Tumor_cell") ~ "Tumor_cell",
    str_starts(final_cell_type, "Epithelia_") ~ "Epithelia_",
    str_starts(final_cell_type, "T_CD8") ~ "T_CD8",
    str_starts(final_cell_type, "T_CD4") ~ "T_CD4",
    str_starts(final_cell_type, "T_ST") ~ "T_ST",
    str_starts(final_cell_type, "NK_") ~ "NK_",
    str_starts(final_cell_type, "B_") ~ "B_",
    str_starts(final_cell_type, "Plasma_") ~ "Plasma_",
    str_starts(final_cell_type, "Fib_") ~ "Fib_",
    str_starts(final_cell_type, "EC_") ~ "EC_",
    str_starts(final_cell_type, "SMC") ~ "SMC",
    str_starts(final_cell_type, "Pericyte_") ~ "Pericyte_",
    str_starts(final_cell_type, "Mph_") ~ "Mph_",
    str_starts(final_cell_type, "Mono_") ~ "Mono_",
    str_starts(final_cell_type, "DC_") ~ "DC_",
    str_starts(final_cell_type, "Myeloid") ~ "Myeloid",
    str_starts(final_cell_type, "Neu_") ~ "Neu_",
    str_starts(final_cell_type, "Mast") ~ "Mast",
    str_starts(final_cell_type, "Glial_") ~ "Glial_",
    str_starts(final_cell_type, "Hepatocyte") ~ "Hepatocyte",
    TRUE ~ "Other"
  )) %>% 
  column_to_rownames(var = "final_cell_type")
scRNA_cross_tissue_proportion_transform$main_type <- factor(scRNA_cross_tissue_proportion_transform$main_type,levels = order_prefix)
scRNA_cross_tissue_proportion_transform <- scRNA_cross_tissue_proportion_transform[order(scRNA_cross_tissue_proportion_transform$main_type),]
scRNA_cross_tissue_proportion_transform <- scRNA_cross_tissue_proportion_transform[,-ncol(scRNA_cross_tissue_proportion_transform)]


# proportion boxplot ------------------------------------------------------
scRNA_cross_tissue_proportion_transform_boxplot <- scRNA_cross_tissue_proportion_transform %>% 
  t(.) %>% 
  as.data.frame(.) %>% 
  rownames_to_column(var = "sample_ID")
scRNA_cross_tissue_proportion_transform_boxplot <- melt(scRNA_cross_tissue_proportion_transform_boxplot)
colnames(scRNA_cross_tissue_proportion_transform_boxplot) <- c("sample_ID","cell_type","proportion")

ggplot(scRNA_cross_tissue_proportion_transform_boxplot, aes(cell_type,proportion),color = "#F47D6E") +
  geom_boxplot(fill = "#F47D6E",outlier.shape = NA,alpha = 0.2) +
  geom_jitter_rast(shape = 16, position = position_jitter(0.2),alpha = 0.5,size = 0.5) +
  scale_x_discrete(limits = rownames(scRNA_cross_tissue_proportion_transform)) +
  nrc_theme + NoLegend() + xlab("")
function_plot(filename_prefix = str_glue("{output_dir}/boxplot_diff_status/scRNA_{celltype}_diff_tissue_percentage"), width = 120, height = 300)

# hub cluster confirm -----------------------------------------------------
## scale data
scRNA_cross_tissue_proportion_transform <- scale(scRNA_cross_tissue_proportion_transform)

## cluster confirm
dist_mat <- dist(t(scRNA_cross_tissue_proportion_transform))
hc <- hclust(dist_mat, method = "ward.D2")
clusters <- cutree(hc, k = 20)
cluster_label <- paste0(clusters)
cluster_df <- data.frame(
  sample_ID = names(clusters),
  EC_class = cluster_label,
  row.names = NULL
)
cluster_df$EC_class <- as.numeric(cluster_df$EC_class)
cluster_df <- as.data.frame(cluster_df) %>% 
  arrange(EC_class)

## final confirm
cluster_df <- cluster_df %>%
  mutate(Ecotype = case_when(
    EC_class %in% c(1) ~ "ET_01",
    EC_class %in% c(5) ~ "ET_02",
    EC_class %in% c(12) ~ "ET_03",
    EC_class %in% c(3) ~ "ET_04",
    EC_class %in% c(11) ~ "ET_05",
    EC_class %in% c(17,4,8,10,9,13) ~ "ET_06",
    EC_class %in% c(18,2,7,19) ~ "ET_07",
    EC_class %in% c(6,15,14,20,16) ~ "ET_08"
    # MetaProgramName %in% c() ~ "EC_09",
    # MetaProgramName %in% c() ~ "EC_10"
  ))


# result plot -------------------------------------------------------------
## ano data
cluster_ano <- group_data[,c("sample_ID","study_ID","tumor_code","tissue_type","tumor_status")] %>% 
  filter(sample_ID %in% unique(scRNA_cross_tissue_proportion_filter$sample_ID)) %>% 
  left_join(cluster_df[,c(1,3)],by = "sample_ID") %>% 
  arrange(Ecotype) %>% 
  column_to_rownames(var = "sample_ID")


svg(str_glue("{output_dir}/cross_tissue_ecotypes.svg"),width = 35,height = 16)
ComplexHeatmap::pheatmap(
  as.matrix(scRNA_cross_tissue_proportion_transform[,rownames(cluster_ano)]),
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
  gaps_col = c(142,186,203,234,274,384,459),
  annotation_colors = list(tumor_code = c(tumor_type_color),Ecotype = c(Ecotype_color),
                           tissue_type = c(tissue_type_color),tumor_status = c(tumor_status_color),study_ID = c(study_database_color)),
  legend = TRUE,
  breaks = seq(-0.5, 2, length.out = 333)
)
dev.off()

pdf(str_glue("{output_dir}/cross_tissue_ecotypes.pdf"),width = 35,height = 16)
ComplexHeatmap::pheatmap(
  as.matrix(scRNA_cross_tissue_proportion_transform[,rownames(cluster_ano)]),
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
  gaps_col = c(142,186,203,234,274,382,457),
  annotation_colors = list(tumor_code = c(tumor_type_color),Ecotype = c(Ecotype_color),
                           tissue_type = c(tissue_type_color),tumor_status = c(tumor_status_color),study_ID = c(study_database_color)),
  legend = TRUE,
  breaks = seq(-0.7, 1.5, length.out = 333)
)
dev.off()



# proportion bar plot -----------------------------------------------------
## cell type confirm
barplot_meta_data <- read_csv(str_glue(" "))
barplot_meta_data_filter <- barplot_meta_data %>%
  mutate(barplot_cell_type = case_when(
    grepl("tumor", final_cell_type) ~ "Tumor_cell",
    grepl("^B_", final_cell_type) ~ "B_cell",
    grepl("^EC", final_cell_type) ~ "Endothelia",
    grepl("^Plasma_", final_cell_type) ~ "Plasma_cell",
    grepl("^Fib_", final_cell_type) ~ "Fibroblast",
    grepl("^Ep_", final_cell_type) ~ "Epithelia",
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

## proportion stat
barplot_meta_data_filter <- barplot_meta_data_filter %>%
  filter(final_cell_type != "Other") %>% 
  group_by(sample_ID, barplot_cell_type) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(clust_total = sum(count)) %>%
  dplyr::mutate(clust_prop = count / clust_total * 100) %>% 
  arrange(sample_ID, desc(clust_prop))


## data plot
ggplot(barplot_meta_data_filter, aes(x = sample_ID, y = clust_prop, fill = barplot_cell_type)) + 
  geom_bar(stat = "identity",position = "fill") +
  theme_bw() +
  labs(x = "", y = "Proportion of Cell Type(%)") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = cell_type_color) +
  scale_x_discrete(limits = rownames(cluster_ano)) +
  nrc_theme + NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/proportion_barplot"), width = 400, height = 130)

# data save ---------------------------------------------------------------
cluster_save <- cluster_ano %>% 
  rownames_to_column(var = "sample_ID")
write_csv(cluster_save,str_glue("{output_dir}/cluster_ano.csv"))

