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
library("Startrac")
library("ggrastr")
library("fgsea")
library("GeneNMF")
# library("scPDtools")
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

# color -------------------------------------------------------------------
new_status_color <- c("NAT"="#9BC53D","Primary_tumor"="#5BC0EB",
                      "Brain_metastasis"="#C69D3B","Liver_metastasis"="#EFC971","Lymph_node_metastasis"="#48A83E","Other_metastasis"="#7ABC68",
                      "Brain_Healthy"="#D86127","Liver_Healthy"="#F47D6E","Lymph_node_Healthy"="#F49897")
tumor_type_color <- c("ACC"="#5A7B8F","ANS"="#EE4C97","BRCA"="#3D806F","CESC"="#A08634","CRC"="#F37C95","ESCC"="#608541",
                      "HNSC"="#7D4E57","KIRC"="#BC3C29","LUAD"="#958056","LUSC"="#9FAFA3","NPC"="#6F99AD","PAAD"="#0072B5",
                      "PNET"="#CFC59A","SARC"="#E18727","SKCM"="#FFDC91","GC"="#718DAE","TGCT"="#7876B1","THCA"="#F9AC93",
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
CM_color <- c("CM1"="","CM2"="","CM3"="","CM4"="","CM5"="","CM6"="","CM7"="","CM8"="","CM9"="")


# data dir ----------------------------------------------------------------
input_dir <- " "
all_gene_dir <- " "
output_dir <- " "

# CM  ---------------------------------------------------------------------
CM1 <- c("B_09_NR4A2","B_06_PDE4D","B_05_ANK3","B_12_RASSF6","B_14_CCSER1","B_03_CMSS1","T_CD4_03_ANK3","B_01_COL19A1","B_08_IGHM",
         "T_CD4_02_TSHZ2","Neu_06_CMTM2","Neu_02_CPPED1","Neu_01_S100A12","Neu_07_IFIT3")
CM2 <- c("EC_09_ADAMTSL1","EC_12_PGF","Fib_11_STMN1","Fib_06_POSTN","Fib_01_RUNX2","Fib_05_KIF26B","Mph_03_SLC16A10","Mph_10_HSPA1B",
         "EC_04_PDGFD","Pericyte_02_HTR1F","Mph_04_SELENOP_KCNMA1")
CM3 <- c("Plasma_01_IGHA2_IGHG3","Plasma_06_HSPA1B","Plasma_07_STMN1","T_STMN1","Fib_03_MMP11","T_CD4_Treg_02_TNFRSF18",
         "Mph_01_SELENOP_CCL18","T_CD8_08_HSPA1B","T_CD4_08_HSPA1B","Myeloid_STMN1","Mph_15_MT","Mph_07_SPP1")
CM4 <- c("Plasma_05_CD74","T_CD4_06_CXCL13","T_CD4_Treg_01_IKZF2","Mph_06_IL1B","Mph_09_ISG15","T_CD8_12_ISG15","B_16_ISG15",
         "T_CD4_10_ISG15","B_15_SOX5","T_CD8_01_CXCL13")
CM5 <- c("T_CD4_04_RTKN2","B_13_STMN1","B_10_RGS13","DC_02_CLEC9A","DC_01_CD1C","NK_05_SELL","T_CD8_09_CCR7","B_02_HSPA1B",
         "T_CD4_05_LTB","T_CD4_01_CCR7","B_04_TCL1A")
CM6 <- c("Plasma_03_CROCC","Plasma_02_SOX5","T_CD8_02_IL7R","Mono_02_IL1B","T_CD8_06_ANK3","T_CD8_05_ARIH1")
CM7 <- c("EC_01_ACKR1","EC_02_CXCL12","SMC","Fib_02_CFD","EC_03_CD36","EC_08_DNAJB1","Pericyte_03_CCL21","Fib_09_PTGDS",
         "Fib_04_PI16","EC_06_PROX1")
CM8 <- c("B_07_EGR1","T_CD4_09_ARIH1","Mast_cell","T_CD4_05_ANXA1","T_CD4_07_PPARG","T_CD8_07_P2RY8","NK_03_SPON2",
         "NK_04_FCGR3A_FGFBP2")
CM9 <- c("NK_01_CD160","T_CD8_10_SLC4A10","Mono_01_VCAN","Mono_03_FCGR3A","T_CD8_13_TRDV2","T_CD8_04_FGFBP2","NK_02_FCGR3A_CX3CR1")
all_CM <- list("CM1"=CM1,"CM2"=CM2,"CM3"=CM3,"CM4"=CM4,"CM5"=CM5,"CM6"=CM6,"CM7"=CM7,"CM8"=CM8,"CM9"=CM9)
all_CM_DEG <- list()

# data read ---------------------------------------------------------------
## DEG read
all_cell_type_DEG <- read_csv(str_glue("{input_dir}/CM_DEGs_marker_all.csv"))

## all gene read
all_gene_data <- read_h5(file = all_gene_dir,
                         assay.name = 'RNA', 
                         target.object = 'seurat')
all_gene <- rownames(all_gene_data)

# hub DEG filter -------------------------------------------
for(CM_type in c("CM1","CM2","CM3","CM4","CM5","CM6","CM7","CM8","CM9")){
  # CM_type <- "CM1"
  
  ## hub DEG filter
  current_DEG <- all_cell_type_DEG %>% 
    filter(cluster %in% CM_type) %>% 
    filter(avg_log2FC >= log2(1.5) & pct.1 >= 0.5) %>% 
    arrange(desc(avg_log2FC)) %>% 
    head(100)
  
  ## hub gene
  hub_DEG <- unique(current_DEG$gene)
  
  ## collect all CM DEG
  all_CM_DEG[[CM_type]] <- hub_DEG
  
  
}

# hub DEG enrichment ------------------------------------------------------
enrich_GO_BP <- lapply(all_CM_DEG, function(program){runGSEA(program, universe=all_gene,
                                                             category = "C5", subcategory = "GO:BP")})
enrich_GO_CC <- lapply(all_CM_DEG, function(program){runGSEA(program, universe=all_gene,
                                                             category = "C5", subcategory = "GO:CC")})
enrich_GO_MF <- lapply(all_CM_DEG, function(program){runGSEA(program, universe=all_gene,
                                                             category = "C5", subcategory = "GO:MF")})
enrich_Hallmarker <- lapply(all_CM_DEG, function(program){runGSEA(program, universe=all_gene,
                                                                  category = "H")})

## enrich save
all_data <- data.frame(
  pathway = as.character(),
  pval = as.numeric(),
  padj = as.numeric(),
  foldEnrichment = as.numeric(),
  overlap = as.numeric(),
  size = as.numeric(),
  overlapGenes = as.character(),
  MP = as.character(),
  type = as.character()
)

hall_data <- data.frame(
  pathway = as.character(),
  pval = as.numeric(),
  padj = as.numeric(),
  foldEnrichment = as.numeric(),
  overlap = as.numeric(),
  size = as.numeric(),
  overlapGenes = as.character(),
  MP = as.character(),
  type = as.character()
)
go_data <- data.frame(
  pathway = as.character(),
  pval = as.numeric(),
  padj = as.numeric(),
  foldEnrichment = as.numeric(),
  overlap = as.numeric(),
  size = as.numeric(),
  overlapGenes = as.character(),
  MP = as.character(),
  type = as.character()
)
for(MP_name in names(all_CM_DEG)){
  # MP_name <- "Alveolar"
  # hall_current <- head(enrich_Hallmarker[[MP_name]],10)
  hall_current <- enrich_Hallmarker[[MP_name]]
  hall_current$MP <- MP_name
  hall_current$type <- "HALL"
  hall_data <- rbind(hall_data,hall_current)
  
  # go_current_BP <- head(enrich_GO_BP[[MP_name]],15)
  go_current_BP <- enrich_GO_BP[[MP_name]]
  # go_current_CC <- head(enrich_GO_CC[[MP_name]],15)
  go_current_CC <- enrich_GO_CC[[MP_name]]
  # go_current_MF <- head(enrich_GO_MF[[MP_name]],15)
  go_current_MF <- enrich_GO_MF[[MP_name]]
  go_current_BP$MP <- MP_name
  go_current_BP$type <- "GO_BP"
  go_current_CC$MP <- MP_name
  go_current_CC$type <- "GO_CC"
  go_current_MF$MP <- MP_name
  go_current_MF$type <- "GO_MF"
  go_data <- rbind(go_data,go_current_BP)
  go_data <- rbind(go_data,go_current_CC)
  go_data <- rbind(go_data,go_current_MF)
  
  ## all combine
  all_data <- rbind(all_data,hall_current)
  all_data <- rbind(all_data,go_current_BP)
  all_data <- rbind(all_data,go_current_CC)
  all_data <- rbind(all_data,go_current_MF)
}

## result plot
# all_data <- read_csv(str_glue("{output_dir}/enrichment_all_pathway.csv"))
hub_pathway <- read_csv(str_glue("{output_dir}/enrichment_all_pathway_self_filter.csv"))
all_data_filter <- all_data %>% 
  filter(pathway %in% unique(hub_pathway$pathway))
all_data_filter$log10p <- -log10(all_data_filter$pval)
ggplot(all_data_filter, aes(x = pathway, y = MP)) +
  geom_point(aes(color = type,size = log10p)) +
  scale_size_continuous(range = c(2,8),breaks = c(0,8,16,24,32,40)) +
  scale_y_discrete(limits = unique(all_data_filter$MP)) +
  scale_x_discrete(limits = unique(hub_pathway$pathway)) +
  nrc_theme + scale_color_manual(values = c("GO_BP"="#FDBF6F","GO_CC"="#B2DF8A","GO_MF"="#A6CEE3","HALL"="#33A02C")) +
  xlab("") + ylab("") 
function_plot(filename_prefix = str_glue("{output_dir}/NMF_MP_enrich_dotplot"), width = 350, height = 180)
ggplot(all_data_filter, aes(x = pathway, y = MP)) +
  geom_point(aes(color = type,size = log10p)) +
  scale_size_continuous(range = c(2,8),breaks = c(0,8,16,24,32,40)) +
  scale_y_discrete(limits = unique(all_data_filter$MP)) +
  scale_x_discrete(limits = unique(hub_pathway$pathway)) +
  nrc_theme + scale_color_manual(values = c("GO_BP"="#FDBF6F","GO_CC"="#B2DF8A","GO_MF"="#A6CEE3","HALL"="#33A02C")) +
  xlab("") + ylab("") + coord_flip()
function_plot(filename_prefix = str_glue("{output_dir}/NMF_MP_enrich_dotplot_transform"), width = 280, height = 250)

# data save ---------------------------------------------------------------
if(T){
  saveRDS(all_CM_DEG,str_glue("{output_dir}/all_CM_DEG.RData"))
  write_csv(all_data,str_glue("{output_dir}/enrichment_all_pathway.csv"))
  write_csv(all_data,str_glue("{output_dir}/enrichment_all_pathway_self_filter.csv"))
}


















