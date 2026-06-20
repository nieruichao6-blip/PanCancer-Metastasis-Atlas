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
library("msigdbr")
library("fgsea")
library("GeneNMF")
library("UCell")
library("SCP")
library("Startrac")
library("corrplot")
library("Hmisc")
library("beyondcell")
library("ComplexHeatmap")
library("ggrastr")
library("ggraph")
(.packages())
# setwd(' ')
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
display.brewer.all()

function_plot <- function(filename_prefix, width, height){
  ggsave(filename = str_c(filename_prefix, '.pdf'), plot = last_plot(), device = 'pdf', width = width, height = height, units = 'mm', dpi = 500)
  ggsave(filename = str_c(filename_prefix, '.svg'), plot = last_plot(), device = 'svg', width = width, height = height, units = 'mm', dpi = 500)
  ggsave(filename = str_c(filename_prefix, '.png'), plot = last_plot(), device = 'png', width = width, height = height, units = 'mm', dpi = 1000)
}

jaccard <- function (a, b) {
  intersection = length ( intersect (a,b))
  union = length (a) + length (b) - intersection
  return (intersection/union)
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

all_drug_MP_cor <- data.frame(
  drug = as.character(),
  MP = as.character(),
  cor = as.numeric(),
  pvalue = as.numeric()
)

# data color ----------------------------------------------------------------
final_MP <- c("Cell.Cycle.G2M"="#e0bc58","Cell.Cycle.G1S"="#64abc0","Cell.Cycle.HMG"="#fab37f","Chromatin"="#e98741",
              "Stress"="#8fc0dc","Hypoxia"="#967568","Stress.in.vitro"="#f2d3ca","Protein.maturation"="#eebd85",
              "EMT.I"="#82c785","MHC"="#edeaa4","Epithelial.Senescence"="#cdaa9f","MYC"="#794976",
              "Respriration"="#bcacd3","Secreted.I"="#889b5d","Cilia.I"="#4e9592","PDAC.classcial"="#dbad5f",
              "Alveolar"="#64ae79","PDAC.related"="#ac5092","Lymphocyte.Activation"="#dc8e97","Cell.Signaling"="#e3d1db",
              "Coagulation"="#74a893","Neuron"="#ac9141",
              "Cell.Cycle.single.nucleus"="#5ac6e9","EMT.II"="#ebce8e",
              "EMT.VI"="#e5c06e","Secreted.II"="#7587b1","Cilia.II"="#c7deef","Colon.related"="#e97371","Other"="#1CB232"
              # ""="#e1a4c6",""="#916ba6",
              # ""="#cb8f82",""="#7db3af",""="#d2e0ac"
)
new_status_color <- c("NAT"="#9BC53D","Primary_tumor"="#5BC0EB",
                      "Brain_metastasis"="#C69D3B","Liver_metastasis"="#EFC971","Lymph_node_metastasis"="#48A83E","Other_metastasis"="#7ABC68",
                      "Brain_Healthy"="#D86127","Liver_Healthy"="#F47D6E","Lymph_node_Healthy"="#F49897")

# data dir ----------------------------------------------------------------
input_dir <- " "
scRNA_tumor_dir <- " "
output_dir <- " "


# data read ---------------------------------------------------------------
## origin data read
scRNA_ano <- readRDS(str_glue("{scRNA_tumor_dir}/scRNA_ano_MP.RData"))
scRNA_tumor_cell <- read_h5(file = str_glue("{input_dir}/scRNA_gene_filter.h5"),
                            assay.name = 'RNA',
                            target.object = 'seurat')

origin_meta <- scRNA_ano@meta.data
origin_meta <- rownames_to_column(origin_meta,var = "cell_name")

scRNA_tumor_cell@meta.data <- scRNA_ano@meta.data
scRNA_tumor_cell@reductions <- scRNA_ano@reductions
remove(scRNA_ano)
gc()

## high gene filter
# scRNA_tumor_cell <- subset(scRNA_tumor_cell, subset = nFeature_RNA >= 1500)
# gc()


# data pre process --------------------------------------------------------
scRNA_tumor_cell <- NormalizeData(scRNA_tumor_cell,normalization.method = 'LogNormalize',scale.factor = 10000)
# scRNA_tumor_cell <- FindVariableFeatures(scRNA_tumor_cell, selection = "vst", nfeatures = 3000)
# scRNA_tumor_cell <- ScaleData(scRNA_tumor_cell,features = head(VariableFeatures(scRNA_tumor_cell),5000),dim = 1:30)
# scRNA_tumor_cell <- subset(scRNA_tumor_cell, tumor_code %in% c("NSCLC","LUAD"))
gc()

# data run ----------------------------------------------------------------
## prepare
DefaultAssay(scRNA_tumor_cell) <- "RNA"
gs <- GetCollection(SSc)

## run
bc <- bcScore(scRNA_tumor_cell,gs,expr.thres = 0.1)
bc <- bcUMAP(bc, pc = 30, k.neighbors = 3, res = 0.2)

bcSignatures(bc, UMAP = "beyondcell", signatures = list(values = "sig-21161"))


# result transform --------------------------------------------------------
## data filter and combine
MP_meta <- bc@meta.data[,c(24:51)] %>% 
  rownames_to_column(var = "cell_name")
drug_score <- bc@scaled %>% 
  t(.) %>% 
  as.data.frame(.) %>% 
  rownames_to_column(var = "cell_name")
MP_drug_score <- left_join(MP_meta,drug_score,by = "cell_name")

## hub name
MP_name_list <- colnames(MP_meta)[-1]
drug_name_list <- colnames(drug_score)[-1]

## drug info save
drug_info <- as.data.frame(FindDrugs(bc,colnames(drug_score)[-1]))

# cor analysis ------------------------------------------------------------
for(drug_name in drug_name_list){
  # drug_name <- "sig-21407"
  for(MP_name in MP_name_list){
    # MP_name <- "Cell.Cycle.G2M"
    # cor_result <- cor.test(MP_drug_score[[drug_name]],MP_drug_score[[MP_name]],method = "spearman")
    cor_result <- cor.test(MP_drug_score[[drug_name]],MP_drug_score[[MP_name]],method = "pearson")
    
    ## result combine
    current_result <- data.frame(
      drug = drug_name,
      MP = MP_name,
      cor = cor_result$estimate,
      pvalue = cor_result$p.value
    )
    all_drug_MP_cor <- rbind(all_drug_MP_cor,current_result)
    
  }
}


# diff metastasis tissue drug ---------------------------------------------
sample_meta <- bc@meta.data[,c(1:18)] %>% 
  rownames_to_column(var = "cell_name")
sample_drug_score <- left_join(sample_meta,drug_score,by = "cell_name") %>% 
  dplyr::select(c("sample_ID",drug_name_list))
sample_drug_score_averge <-sample_drug_score %>% 
  dplyr::group_by(sample_ID) %>% 
  dplyr::summarise(across(where(is.numeric),mean,na.rm = T))
sample_info <- sample_meta[!duplicated(sample_meta$sample_ID),c("sample_ID","tumor_code","tumor_status","metastasis_tissue")] %>% 
  mutate(metastasis_tissue_transform = ifelse(
    tumor_status == "Metastasis" & metastasis_tissue == "Brain","Brain_metastasis",
    ifelse(tumor_status == "Metastasis" & metastasis_tissue == "Liver","Liver_metastasis",
           ifelse(tumor_status == "Metastasis" & metastasis_tissue == "Lymph_node","Lymph_node_metastasis",
                  ifelse(tumor_status == "Tumor","Primary_tumor","Other_metastasis")))
    
  ))

## data plot
sample_drug_score_averge_info <- left_join(sample_info,sample_drug_score_averge,by = "sample_ID")
for(drug_name in drug_name_list){
  # drug_name <- "sig-21327"
  
  ## data filter
  current_sample_drug_score_averge_info <- sample_drug_score_averge_info %>% 
    dplyr::select(c("sample_ID",drug_name,"metastasis_tissue_transform"))
  colnames(current_sample_drug_score_averge_info)[2] <- "score"
  pathway_name <- gsub("[ ,/\\\\\\-]","_",drug_name)
  
  ## data plot
  ggplot(current_sample_drug_score_averge_info, aes(metastasis_tissue_transform,score,color = metastasis_tissue_transform)) +
    geom_boxplot(aes(fill = metastasis_tissue_transform),outlier.shape = NA,alpha = 0.2) +
    geom_jitter_rast(shape = 16, position = position_jitter(0.2),alpha = 0.5,size = 3) +
    scale_color_manual(values = new_status_color) +
    scale_fill_manual(values = new_status_color) +
    scale_x_discrete(limits = c("Primary_tumor","Brain_metastasis","Liver_metastasis","Lymph_node_metastasis","Other_metastasis")) +
    stat_compare_means(comparisons = list( c("Primary_tumor", "Brain_metastasis"), c("Primary_tumor", "Liver_metastasis"),
                                           c("Primary_tumor", "Lymph_node_metastasis"),c("Primary_tumor", "Other_metastasis")),
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),symbols = c("***", "**", "*", "ns")),
                       method="wilcox.test",label = "p.signif") +
    nrc_theme + NoLegend() + xlab("")
  function_plot(filename_prefix = str_glue("{output_dir}/drug_boxplot/drug_{pathway_name}"), width = 100, height = 150)
}

# data plot ---------------------------------------------------------------
all_drug_MP_cor_transform <- all_drug_MP_cor[,c(1,2,3)]
all_drug_MP_cor_transform <- all_drug_MP_cor_transform[!grepl("CancerSEA",all_drug_MP_cor_transform$drug),]
all_drug_MP_cor_transform <- all_drug_MP_cor_transform %>% 
  pivot_wider(names_from = MP,values_from = cor) %>% 
  column_to_rownames(var = "drug") %>% 
  as.matrix(.) %>% 
  t(.)

if(T){
  dominant_clusters <- apply(all_drug_MP_cor_transform, 2, function(x) rownames(all_drug_MP_cor_transform)[which.max(x)])
  ordered_samples <- character(0)
  for (cluster in rownames(all_drug_MP_cor_transform)) {
    candidate_samples <- names(dominant_clusters)[dominant_clusters == cluster]
    candidate_samples <- setdiff(candidate_samples, ordered_samples)
    samples_sorted <- names(sort(all_drug_MP_cor_transform[cluster, candidate_samples], decreasing = TRUE))
    
    ordered_samples <- c(ordered_samples, samples_sorted)
  }
  sorted_mat <- all_drug_MP_cor_transform[, ordered_samples]
  
}

all_cor_filter_3 <- all_drug_MP_cor
all_cor_filter_3 <- all_cor_filter_3[!grepl("CancerSEA",all_cor_filter_3$drug),]
all_cor_filter_3 <- all_cor_filter_3[!grepl("WHITFIELD",all_cor_filter_3$drug),]
all_cor_filter_3 <- all_cor_filter_3[!grepl("EMT",all_cor_filter_3$drug),]
all_cor_filter_3 <- all_cor_filter_3[!grepl("senescence",all_cor_filter_3$drug),]
all_cor_filter_3 <- all_cor_filter_3 %>% 
  dplyr::group_by(MP) %>%
  dplyr::arrange(desc(cor),.by_group = T) %>% 
  slice_head(n = 2)
ano_col <- data.frame(row.names = unique(all_drug_MP_cor$MP),
                                  final_MP=unique(all_drug_MP_cor$MP))

pos <- which(rownames(t(sorted_mat)) %in% unique(all_cor_filter_3$drug))
label <- rowAnnotation(Zscore = anno_mark(at = pos,labels = unique(all_cor_filter_3$drug),labels_gp = gpar(fontSize = 8)))
pdf(str_glue("{output_dir}/all_drug_result.pdf"),width = 15,height = 15)
pheatmap(t(sorted_mat),
         cluster_row = F, cluster_cols = F, show_rownames = F,
         # annotation_row = ano_row,
         # annotation_col = ano_col,
         color = colorRampPalette(c("#3B4CC0","white","#B40426"))(100),
         breaks = seq(-1,1,length.out = 101),
         show_colnames = T, legend = T, annotation_names_row = F, annotation_names_col = F,
         annotation_col = ano_col,
         right_annotation = label,
         cellwidth = 20,cellheight = 1.5,
         fontsize = 5,
         annotation_colors = list(final_MP = c(final_MP)),
         use_raster = F
)
dev.off()

# data save ---------------------------------------------------------------
all_drug_MP_save <- all_drug_MP_cor %>%
  group_by(MP) %>%
  arrange(desc(cor)) %>%
  slice_head(n = 50) %>%
  ungroup()

# bc <- readRDS(str_glue("{output_dir}/bc.RData"))
# all_drug_MP_cor <- read_csv(str_glue("{output_dir}/all_drug_MP_cor_origin.csv"))
saveRDS(bc,str_glue("{output_dir}/bc.RData"))
all_drug_MP_cor <- read_csv(str_glue("{output_dir}/all_drug_MP_cor_origin.csv"))
write_csv(all_drug_MP_cor,str_glue("{output_dir}/all_drug_MP_cor_origin.csv"))
write_csv(all_drug_MP_save,str_glue("{output_dir}/all_drug_MP_cor_filter.csv"))
write_csv(drug_info,str_glue("{output_dir}/drug_info.csv"))







