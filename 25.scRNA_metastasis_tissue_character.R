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
library("progeny")
library("gghalves")
library("gridExtra")
library("scMetabolism")
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
input_dir <- " "
output_dir <- " "
hall_dir <- " "
tumor_code_list <- c("ACC","BRCA","CESC","CRC","ESCC","HNSC","KIRC","NPC","PAAD","PNET","GC",
                     "THCA","NSCLC","LSCC")


# data color --------------------------------------------------------------
tumor_type_color <- c("ACC"="#5A7B8F","ANS"="#EE4C97","BRCA"="#3D806F","CESC"="#A08634","CRC"="#F37C95","ESCC"="#608541",
                      "HNSC"="#7D4E57","KIRC"="#BC3C29","LUAD"="#958056","LUSC"="#9FAFA3","NPC"="#6F99AD","PAAD"="#0072B5",
                      "PNET"="#CFC59A","SARC"="#E18727","SKCM"="#FFDC91","GC"="#718DAE","TGCT"="#7876B1","THCA"="#F9AC93",
                      "NSCLC"="#E3BC06","LSCC"="#2285FE","Liver_healthy"="#3E6086","Lymph_node_healthy"="#20854E","Brain_healthy"="#4A7985")
tissue_type_color <- c("Bone_marrow"="#CE1261","Kidney"="#5E4FA2","Lymph_node"="#8CA77B","Brain"="#00441B","Lung"="#ed3d30",
                       "Colon"="#DEDC00","Liver"="#B3DE69","Peritoneal"="#999999","Stomach"="#725e82","Ovary"="#8c4356",
                       "Head_and_neck"="#db5a6b","Pancreas"="#21a675","Vaginal"="#d9b611","Appendix"="#e29c45",
                       "Bone"="#4c8dae","Nasopharynx"="#003371","Salivary_gland"="#a1afc9","Testis"="#88ada6","Cervix"="#c93756",
                       "Skin"="#204B75","Breast"="#588257","Thyroid"="#B6DB7B","Esophagus"="#4285BF","larynx"="#1CB232")


# data read ---------------------------------------------------------------
genesets <- getGmt(str_glue("{hall_dir}/h.all.v2023.1.Hs.symbols.gmt"))
scRNA_tumor_cell <- read_h5(file = str_glue("{input_dir}/scRNA_bbknn.h5"),
                            assay.name = 'RNA', 
                            target.object = 'seurat')
scRNA_tumor_cell@meta.data$tumor_tissue_type <- str_c(scRNA_tumor_cell@meta.data$tumor_code,
                                                      scRNA_tumor_cell@meta.data$tissue_type,
                                                      scRNA_tumor_cell@meta.data$tumor_status,sep = "_")

# data process ------------------------------------------------------------
Idents(scRNA_tumor_cell) <- scRNA_tumor_cell@meta.data$tumor_tissue_type

## avg exp
ave_exp <- as.matrix(AverageExpression(scRNA_tumor_cell)[[1]])
ave_exp <- ave_exp[rowSums(ave_exp) > 1,]

## GSVA
# gsvaPar <- gsvaParam(ave_exp,genesets,kcdf = "Gaussian",maxDiff = T,sparse = F)
hall_GSVA <- gsva(ave_exp,genesets,kcdf = "Gaussian",mx.diff = T)

## result combine
hall_GSVA <- as.data.frame(hall_GSVA)

## col name filter
col_sorted <- c("BRCA_Breast_Tumor","GC_Stomach_Tumor","CRC_Colon_Tumor","ACC_Salivary_gland_Tumor","CESC_Cervix_Tumor",
                "LUAD_Lung_Tumor","PAAD_Pancreas_Tumor","THCA_Thyroid_Tumor","NPC_Nasopharynx_Tumor","KIRC_Kidney_Tumor",
                "ESCC_Esophagus_Tumor","HNSC_Head_and_neck_Tumor","LSCC_Larynx_Tumor","PNET_Pancreas_Tumor",
                "LSCC_Lymph_node_Metastasis","TGCT_Lymph_node_Metastasis","BRCA_Lymph_node_Metastasis","ESCC_Lymph_node_Metastasis",
                "HNSC_Lymph_node_Metastasis","LUAD_Lymph_node_Metastasis","LUSC_Lymph_node_Metastasis","THCA_Lymph_node_Metastasis",
                "BRCA_Brain_Metastasis","CRC_Brain_Metastasis","NSCLC_Brain_Metastasis","ESCC_Brain_Metastasis","LUAD_Brain_Metastasis",
                "THCA_Brain_Metastasis",
                "CESC_Liver_Metastasis","CRC_Liver_Metastasis","NPC_Liver_Metastasis","PAAD_Liver_Metastasis","BRCA_Liver_Metastasis",
                "PNET_Liver_Metastasis","GC_Liver_Metastasis",
                "GC_Ovary_Metastasis","CRC_Peritoneal_Metastasis","GC_Peritoneal_Metastasis","PAAD_Lung_Metastasis","ACC_Lung_Metastasis",
                "PAAD_Vaginal_Metastasis","ANS_Peritoneal_Metastasis","KIRC_Bone_marrow_Metastasis")

hall_GSVA_sorted <- hall_GSVA %>% 
  dplyr::select(col_sorted)

## data plot
pdf(str_glue("{output_dir}/GSVA_pathway_hallmarker.pdf"),width = 10,height = 7)
pheatmap(as.matrix(hall_GSVA_sorted),scale = "row",
         show_colnames = T,fontsize_row = 8, fontsize_col =8,cluster_cols= F,cluster_rows = T,
         color = colorRampPalette(c("#2285FE", "white", "#e53f31"))(50),border_color = "black",
         gaps_col = c(14,22,28,35))
dev.off()


# scMetabolism process ----------------------------------------------------
kegg_meta_geneset <- getGmt(str_glue("{hall_dir}/KEGG_metabolism_nc.gmt"))
reactome_meta_geneset <- getGmt(str_glue("{hall_dir}/REACTOME_metabolism.gmt"))

## KEGG meta
# gsvaPar_kegg <- gsvaParam(ave_exp,kegg_meta_geneset,kcdf = "Gaussian",maxDiff = T,sparse = F)
hall_GSVA_kegg <- gsva(ave_exp,kegg_meta_geneset,kcdf = "Gaussian",mx.diff = T)
hall_GSVA_kegg <- as.data.frame(hall_GSVA_kegg)
hall_GSVA_kegg_sorted <- hall_GSVA_kegg %>% 
  dplyr::select(col_sorted)
pdf(str_glue("{output_dir}/GSVA_pathway_kegg_meta.pdf"),width = 11,height = 12)
pheatmap(as.matrix(hall_GSVA_kegg_sorted),scale = "row",
         show_colnames = T,fontsize_row = 8, fontsize_col =8,cluster_cols= F,cluster_rows = T,
         color = colorRampPalette(c("#2285FE", "white", "#e53f31"))(50),border_color = "black",
         gaps_col = c(14,22,28,35))
dev.off()

## REACTOME meta
hall_GSVA_reactome <- gsva(ave_exp,reactome_meta_geneset,kcdf = "Gaussian",mx.diff = T)
hall_GSVA_reactome <- as.data.frame(hall_GSVA_reactome)
hall_GSVA_reactome_sorted <- hall_GSVA_reactome %>% 
  dplyr::select(col_sorted)
pdf(str_glue("{output_dir}/GSVA_pathway_reactome_meta.pdf"),width = 11,height = 12)
pheatmap(as.matrix(hall_GSVA_reactome_sorted),scale = "row",cellwidth = 9,
         show_colnames = T,fontsize_row = 8, fontsize_col =8,cluster_cols= F,cluster_rows = T,
         color = colorRampPalette(c("#2285FE", "white", "#e53f31"))(50),border_color = "black",
         gaps_col = c(14,22,28,35))
dev.off()


# scRNA_tumor_cell_meta <- subset(scRNA_tumor_cell,tumor_status == "Metastasis")
# scRNA_tumor_Metabolism <- sc.metabolism.Seurat(obj = scRNA_tumor_cell_meta,
#                                                method = "VISION",
#                                                imputation = F,
#                                                ncores = 5,
#                                                metabolism.type = "KEGG")
# DotPlot.metabolism(scRNA_tumor_Metabolism)


# PROGENy process ---------------------------------------------------------------
## run
scRNA_tumor_cell <- progeny(scRNA_tumor_cell, scale=FALSE, organism="Human", top=500, perm=1, return_assay = TRUE)
scRNA_tumor_cell <- Seurat::ScaleData(scRNA_tumor_cell, assay = "progeny")

## result transform
progeny_scores <- as.data.frame(t(GetAssayData(scRNA_tumor_cell, slot = "scale.data", assay = "progeny"))) %>%
    rownames_to_column("Cell") %>%
    gather(Pathway, Activity, -Cell) 

## add info
CellsClusters <- data.frame(Cell = names(Idents(scRNA_tumor_cell)), 
    CellType = as.character(Idents(scRNA_tumor_cell)),
    stringsAsFactors = FALSE)
progeny_scores <- inner_join(progeny_scores, CellsClusters)

## stat
summarized_progeny_scores <- progeny_scores %>% 
    dplyr::group_by(Pathway, CellType) %>%
    dplyr::summarise(avg = mean(Activity), std = sd(Activity))
summarized_progeny_scores_df <- summarized_progeny_scores %>%
    dplyr::select(-std) %>%   
    spread(Pathway, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
summarized_progeny_scores_df <- summarized_progeny_scores_df %>% 
    t(.) %>% 
    as.data.frame(.) %>% 
    dplyr::select(col_sorted) %>% 
    rownames_to_column(var = "pathway") %>% 
    # filter(pathway %in% c("TNFa","NFkB","TGFb","p53","JAK-STAT","Androgen","WNT","Trail","PI3K","Estrogen")) %>% 
    column_to_rownames(var = "pathway")

## data plot
paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(as.matrix(summarized_progeny_scores_df),fontsize=14, 
                        fontsize_row = 10, 
                        color=colorRampPalette(c("#2285FE", "white", "#e53f31"))(100), 
                        # breaks = progenyBreaks, 
                        main = "PROGENy (500)", angle_col = 45,
                        treeheight_col = 0,  border_color = NA)
pheatmap(as.matrix(summarized_progeny_scores_df),
         show_colnames = T,fontsize_row = 8, fontsize_col =8,cluster_cols= F,cluster_rows = T,
         color = colorRampPalette(c("#2285FE", "white", "#e53f31"))(100),border_color = "black",
         breaks = progenyBreaks, 
         gaps_col = c(14,22,28,35))



# data save ---------------------------------------------------------------
hall_GSVA_save <- rownames_to_column(hall_GSVA,var = "pathway")

if(T){
  write_csv(hall_GSVA_save,str_glue("{output_dir}/hall_GSVA_result.csv"))
  
}











