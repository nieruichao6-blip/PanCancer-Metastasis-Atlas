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
library("Startrac")
library("ggrastr")
library("ggraph")
library("tidygraph")
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
cell_type_color <- c("Epithelia_tumor"="#E41A1C",
                     "Melanoma_tumor"="#E41A1C",
                     "Osteoblastic_tumor"="#E41A1C",
                     "Epithelia"="#A6CEE3",
                     "Fibroblast"="#FDBF6F",
                     "Endothelia"="#B2DF8A",
                     "Acinar_cell"="#984EA3",
                     "SMC&Pericyte"="#33A02C",
                     "T_cell"="#1F78B4",
                     "NK_cell"="#64abc0",
                     "B_cell"="#FF7F00",
                     "Plasma_cell"="#FB9A99",
                     "Myeloid_cell"="#CAB2D6",
                     "Mast_cell"="#6A3D9A",
                     "Neutrophils"="#FFFF99",
                     "Hepatocyte"="#9ec9e1",
                     "Melanoma_cell"="#E66F00",
                     "Neuron"="#1d92c0",
                     "Glial_cell"="#B15928",
                     "Osteoblastic_cell"="#fcc5c1",
                     "Alpha_cell"="#42aa5e",
                     "Bela_cell"="#c22b86"
)
CM_color <- c("CM1"="#ea5c6f","CM2"="#f7905a","CM3"="#3187cb","CM4"="#fb948d","CM5"="#e2b159",
              "CM6"="#ebed6f","CM7"="#b2db87","CM8"="#7ee7bb","CM9"="#64cccf")
new_status_color <- c("NAT"="#9BC53D","Primary_tumor"="#5BC0EB",
                      "Brain_metastasis"="#C69D3B","Liver_metastasis"="#EFC971","Lymph_node_metastasis"="#48A83E","Other_metastasis"="#7ABC68",
                      "Brain_Healthy"="#D86127","Liver_Healthy"="#F47D6E","Lymph_node_Healthy"="#F49897")
tumor_type_color <- c("ACC"="#5A7B8F","ANS"="#EE4C97","BRCA"="#3D806F","CESC"="#A08634","CRC"="#F37C95","ESCC"="#608541","MA"="#FDC44D",
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
correlation_dir <- " "
input_dir <- " "
output_dir <- " "


# CM  ---------------------------------------------------------------------
CM1 <- c("B_09_NR4A2","B_06_PDE4D","B_05_ANK3","B_12_RASSF6","B_14_CCSER1","B_03_CMSS1","T_CD4_03_ANK3","B_01_COL19A1","B_08_IGHM",
         "T_CD4_02_TSHZ2","Neu_06_CMTM2","Neu_02_CPPED1","Neu_01_S100A12","Neu_07_IFIT3")
# CM1 <- c("B_09_NR4A2","B_06_PDE4D","B_05_ANK3","B_12_RASSF6","B_14_CCSER1","B_03_CMSS1","T_CD4_03_ANK3","B_01_COL19A1","B_08_IGHM",
#          "T_CD4_02_TSHZ2")
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
all_CM <- c(CM1,CM2,CM3,CM4,CM5,CM6,CM7,CM8,CM9)
all_CM_list <- list("CM1"=CM1,"CM2"=CM2,"CM3"=CM3,"CM4"=CM4,"CM5"=CM5,"CM6"=CM6,"CM7"=CM7,"CM8"=CM8,"CM9"=CM9)

## CM dataframe save
CM_save <- bind_rows(lapply(all_CM_list,function(x) as.data.frame(t(x))))
CM_save <- as.data.frame(t(CM_save))
colnames(CM_save) <- names(all_CM_list)
write_csv(CM_save,str_glue("{output_dir}/CM_save.csv"))

# data read ---------------------------------------------------------------
cell_type_correlation <- read_csv(str_glue("{correlation_dir}/hub_cell_correlation.csv"))
cell_type_pvalue <- read_csv(str_glue("{correlation_dir}/hub_cell_pvalue.csv"))
scRNA_all_cell_meta <- read_csv(str_glue("{input_dir}/all_cell_combine_meta.csv"))

## add CM meta
scRNA_all_cell_meta_CM <- scRNA_all_cell_meta %>% 
  mutate(CM = ifelse(final_cell_type %in% CM1,"CM1",
                     ifelse(final_cell_type %in% CM2,"CM2",
                            ifelse(final_cell_type %in% CM3,"CM3",
                                   ifelse(final_cell_type %in% CM4,"CM4",
                                          ifelse(final_cell_type %in% CM5,"CM5",
                                                 ifelse(final_cell_type %in% CM6,"CM6",
                                                        ifelse(final_cell_type %in% CM7,"CM7",
                                                               ifelse(final_cell_type %in% CM8,"CM8",
                                                                      ifelse(final_cell_type %in% CM9,"CM9","Other"))))))))))
scRNA_all_cell_meta_CM_filter <- scRNA_all_cell_meta_CM %>%
  filter(CM != "Other")


# CM OR -------------------------------------------------------------------
CM_tissue_OR <- calTissueDist(scRNA_all_cell_meta_CM,byPatient = F,colname.cluster = "CM",colname.patient = "sample_ID",
                              colname.tissue = "tissue_type",method = "chisq",min.rowSum = 0)
CM_tissue_OR <- CM_tissue_OR[c("CM1","CM2","CM3","CM4","CM5","CM6","CM7","CM8","CM9"),]
pdf(str_glue("{output_dir}/tissue_OR_pheatmap.pdf"), width = ncol(CM_tissue_OR)*0.4+2, height = nrow(CM_tissue_OR)*0.54+1)
pheatmap(CM_tissue_OR,
         cluster_rows = T,
         color = c("#FEE8C8", "#FDBB84", "#FC8D59", "#EF6548"),
         breaks = c(0, 1, 1.5, 3, 4),
         cluster_cols = T,
         # angle_col = 45,
         fontsize = 18, border_color = "white",
         display_numbers = matrix(ifelse(CM_tissue_OR > 3, "+++", ifelse(CM_tissue_OR > 1.5, "++", ifelse(CM_tissue_OR > 1, "+", ""))), nrow(CM_tissue_OR)),
         number_color = "black"
)
dev.off()

CM_tissue_OR_save <- as.data.frame(CM_tissue_OR)
colnames(CM_tissue_OR_save) <- c("CM_type","tissue_type","OR")
write_csv(CM_tissue_OR_save,str_glue("{output_dir}/CM_tissue_OR.csv"))

# CM correlation network plot ---------------------------------------------
for(CM_type in names(all_CM_list)){
  CM_type <- "CM1"
  
  ## data filter
  current_correlation <- cell_type_correlation %>% 
    column_to_rownames(var = "cell_type")
  current_correlation <- current_correlation[all_CM_list[[CM_type]],all_CM_list[[CM_type]]]
  current_pvalue <- cell_type_pvalue %>% 
    column_to_rownames(var = "cell_type")
  current_pvalue <- current_pvalue[all_CM_list[[CM_type]],all_CM_list[[CM_type]]]
  
  ## data transform
  current_correlation_transform <- rownames_to_column(current_correlation,var = "cell_type")
  current_correlation_transform <- melt(current_correlation_transform)
  colnames(current_correlation_transform) <- c("cell_type1","cell_type2","correlation")
  current_correlation_transform$bind_id <- str_c(current_correlation_transform$cell_type1,current_correlation_transform$cell_type2,sep = "_")
  current_pvalue_transform <- rownames_to_column(current_pvalue,var = "cell_type")
  current_pvalue_transform <- melt(current_pvalue_transform)
  colnames(current_pvalue_transform) <- c("cell_type1","cell_type2","pvalue")
  current_pvalue_transform$bind_id <- str_c(current_pvalue_transform$cell_type1,current_pvalue_transform$cell_type2,sep = "_")
  
  ## cor and pavlue cpmbine
  current_network_data <- left_join(current_correlation_transform,current_pvalue_transform[,c("bind_id","pvalue")],by = "bind_id")
  current_network_data[current_network_data$pvalue > 0.5,]$correlation <- NA
  
  ## plot prepare
  nodeDF <- as.data.frame(unique(current_network_data$cell_type1))
  nodeDF$id <- c(1:nrow(nodeDF))
  colnames(nodeDF)[1] <- "cell_type"
  nodeDF$label <- nodeDF$cell_type
  nodeDF <- nodeDF %>%
    mutate(hub_cell_type = case_when(
      grepl("^B_|^Plasma_", cell_type) ~ "B_cell",
      grepl("^^Plasma_", cell_type) ~ "Plasma_cell",
      grepl("^EC", cell_type) ~ "Endothelia",
      grepl("^Fib_", cell_type) ~ "Fibroblast",
      grepl("^Pericyte|^SMC", cell_type) ~ "SMC&Pericyte",
      grepl("^T_", cell_type) ~ "T_cell",
      grepl("^^NK_", cell_type) ~ "NK_cell",
      grepl("^Mono_|^Mph_|^DC_|^Neu_|^Mast|^Myeloid", cell_type) ~ "Myeloid_cell",
      grepl("^Mast", cell_type) ~ "Mast_cell",
      grepl("^Neu_", cell_type) ~ "Neutrophils",
      TRUE ~ "Other"
    ))
  linDF <- current_network_data
  linDF <- linDF[!is.na(linDF$correlation),]
  # linDF[linDF$correlation >= 0.5,]$correlation <- 0.5
  linDF <- linDF[linDF$correlation >= 0.1,]
  colnames(linDF)[1:2] <- c("from","to")
  
  graph_data <- graph_from_data_frame(linDF, vertices = nodeDF, directed = FALSE)
  
  ## result plot
  ggraph(graph_data,layout = "linear", circular = TRUE) +
    geom_edge_link(aes(edge_colour = correlation),edge_width = 1) +
    geom_node_point(aes(fill = hub_cell_type,color = hub_cell_type),size = 24,shape = 21) +
    geom_node_text(aes(label = label),size = 2,repel = FALSE) +
    scale_edge_color_gradient2(low = "#6D9EC1",mid = 'white', high = "#E46726") + 
    scale_fill_manual(values = cell_type_color) +
    scale_color_manual(values = cell_type_color) +
    theme_void() +
    NoLegend()
    # scale_edge_color_continuous(range = c(0,0.5))
  function_plot(filename_prefix = str_glue("{output_dir}/CM_network_plot/{CM_type}_network"), width = 100, height = 100)
  
  
}



# CM in status ------------------------------------------------------------
## sample percentage
CM_status_stat <- scRNA_all_cell_meta_CM %>% 
  dplyr::group_by(sample_ID,CM) %>% 
  dplyr::summarise(count = n(), .groups = "drop") %>% 
  dplyr::group_by(sample_ID) %>% 
  dplyr::mutate(percentage = count / sum(count)*100) %>% 
  ungroup

CM_status_stat <- CM_status_stat %>%
  filter(CM != "Other")

## hub sample with zero
CM_status_stat_transform <- CM_status_stat[,c(1,2,4)] %>%
  pivot_wider(names_from = CM,
              values_from = percentage,
              values_fill = 0)
CM_status_stat_final <- melt(CM_status_stat_transform)
colnames(CM_status_stat_final) <- c("sample_ID","CM","percentage")

write_csv(CM_status_stat_transform,str_glue("{output_dir}/CM_status_stat_transform.csv"))

## combine sample info
CM_status_stat_final <- left_join(CM_status_stat_final,scRNA_all_cell_meta[!duplicated(scRNA_all_cell_meta$sample_ID),c("sample_ID","tumor_code","tumor_status","tissue_type","metastasis_tissue")],by = "sample_ID")
CM_status_stat_final <- CM_status_stat_final %>% 
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

## result plot
for(CM_type in as.character(unique(CM_status_stat_final$CM))){
  # CM_type <- "CM1"
  
  ## data filter
  current_meta_stat <- CM_status_stat_final %>% 
    filter(CM == CM_type)
  
  ## data plot
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
  function_plot(filename_prefix = str_glue("{output_dir}/CM_diff_status/scRNA_{CM_type}_diff_tissue_percentage"), width = 120, height = 330)
}


# percentage in sample ----------------------------------------------------
CM_sample_data <- CM_status_stat[,c(1,2,4)] %>% 
  pivot_wider(names_from = sample_ID,
              values_from = percentage,
              values_fill = list(percentage = 0)) %>% 
  column_to_rownames(var = "CM")
# CM_sample_data_filter <- CM_sample_data
CM_sample_data_filter <- CM_sample_data[,colSums(CM_sample_data >= 5) > 0]

## new sorted
if(T){
  dominant_clusters <- apply(CM_sample_data_filter, 2, function(x) rownames(CM_sample_data_filter)[which.max(x)])
  ordered_samples <- character(0)
  for (cluster in rownames(CM_sample_data_filter)) {
    candidate_samples <- names(dominant_clusters)[dominant_clusters == cluster]
    candidate_samples <- setdiff(candidate_samples, ordered_samples)
    samples_sorted <- names(sort(CM_sample_data_filter[cluster, candidate_samples], decreasing = TRUE))
    
    ordered_samples <- c(ordered_samples, samples_sorted)
  }
  sorted_mat <- CM_sample_data_filter[, ordered_samples]
  
}

## data plot
ano_col <- data.frame(sample_ID = colnames(sorted_mat))
sample_info <- scRNA_all_cell_meta_CM[,c("sample_ID","study_ID","tumor_code","tissue_type","tumor_status")]
sample_info <- sample_info[!duplicated(sample_info$sample_ID),]
ano_col <- left_join(ano_col,sample_info)
ano_col <- column_to_rownames(sample_info,var = "sample_ID")
ano_row <- data.frame(row.names = rownames(sorted_mat),CM = rownames(sorted_mat))
pdf(str_glue("{output_dir}/CM_percentage_in_sample.pdf"), width = 15, height = 6)
pheatmap(sorted_mat,
         cluster_row = F, cluster_cols = F, show_rownames = T, 
         annotation_row = ano_row,
         annotation_col = ano_col,
         show_colnames = F, legend = T, annotation_names_row = F, annotation_names_col = F,
         annotation_colors = list(tumor_code = c(tumor_type_color),CM = CM_color,
                                  tissue_type = c(tissue_type_color),tumor_status = c(tumor_status_color),study_ID = c(study_database_color)),
         cellwidth = 1,cellheight = 30,
         fontsize = 5,
         breaks = seq(0,50,length.out = 101),
         use_raster = F
)
dev.off()


# diff tumor status -------------------------------------------------------
for(CM_type in as.character(unique(CM_status_stat_final$CM))){
  # CM_type <- "CM3"
  
  ## data filter
  current_meta_stat <- scRNA_all_cell_meta_CM %>% 
    filter(CM == CM_type)
  
  current_meta_stat <- current_meta_stat %>% 
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
  
  ## data plot
  CM_status_OR <- calTissueDist(current_meta_stat,byPatient = F,colname.cluster = "final_cell_type",colname.patient = "sample_ID",
                                       colname.tissue = "metastasis_tissue_transform",method = "chisq",min.rowSum = 0)
  CM_status_OR <- CM_status_OR[,c("NAT","Primary_tumor","Brain_metastasis","Liver_metastasis","Lymph_node_metastasis","Other_metastasis",
                                  "Brain_Healthy","Liver_Healthy","Lymph_node_Healthy")]
  pdf(str_glue("{output_dir}/CM_diff_cell_type/{CM_type}_pheatmap.pdf"), width = ncol(CM_status_OR)*0.6+3, height = nrow(CM_status_OR)*0.6+1)
  pheatmap(CM_status_OR,
           cluster_rows = TRUE,
           color = c("#FEE8C8", "#FDBB84", "#FC8D59", "#EF6548"),
           breaks = c(0, 1, 1.5, 3, max(CM_status_OR)),
           cluster_cols = FALSE,
           angle_col = 45,
           fontsize = 18, border_color = "white",
           display_numbers = matrix(ifelse(CM_status_OR > 3, "+++", ifelse(CM_status_OR > 1.5, "++", ifelse(CM_status_OR > 1, "+", "+/-"))), nrow(CM_status_OR)),
           number_color = "black"
  )
  dev.off()

}



# CM percentage in diff tumor code ----------------------------------------
CM_tumor_code <- scRNA_all_cell_meta_CM_filter %>% 
  dplyr::group_by(tumor_code,CM) %>% 
  dplyr::summarise(cell_n = n(), .groups = "drop") %>% 
  dplyr::group_by(tumor_code) %>% 
  dplyr::mutate(total_n = sum(cell_n),
                proportion = cell_n/total_n*100) %>% 
  ungroup()

                                                  
metadata_summary_sorted <- CM_tumor_code %>%
  filter(tumor_code %in% unique(CM_tumor_code$tumor_code)[!unique(CM_tumor_code$tumor_code) %in% c("Brain_healthy","Liver_healthy","Lymph_node_healthy")]) %>% 
  mutate(label = ifelse(proportion >= 10, ">10%",
                        ifelse(proportion >= 5, "5%-10%",
                               ifelse(proportion >= 1, "1%-5%","<1%"))))
metadata_summary_sorted$tumor_code <- as.character(metadata_summary_sorted$tumor_code)
metadata_summary_sorted$CM <- factor(metadata_summary_sorted$CM)

mycols <- brewer.pal(9, "GnBu")[c(1, 3, 5, 7)]
names(mycols) <- c("<1%", "1%-5%", "5%-10%", ">10%")
ggplot(data = metadata_summary_sorted, mapping = aes_string(x = "tumor_code", y = "CM")) +
  geom_point_rast(mapping = aes_string(size = "proportion", color = "label")) +
  scale_size_continuous(breaks = c(0,1,5,10,15,20)) +
  scale_x_discrete(limits = c("ACC","ANS","BRCA","CESC","CRC","ESCC","GC","HNSC","KIRC","LSCC","NSCLC","LUAD","LUSC","MA","NPC","PAAD","PNET")) +
  theme_bw() +
  theme(
    strip.text.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(colour = "grey", linetype = "dashed", size = 0.2)
  ) +
  scale_color_manual(values = mycols)
function_plot(filename_prefix = str_glue("{output_dir}/cell_type_proportion_dotplot"), width = 150, height = 70)



# CM in tissue type -------------------------------------------------------
for(CM_type in as.character(unique(CM_status_stat_final$CM))){
  # CM_type <- "CM1"
  
  ## data filter
  current_meta_stat <- CM_status_stat_final %>% 
    filter(CM == CM_type)
  
  ## data plot
  ggplot(current_meta_stat, aes(tissue_type,percentage,color = tissue_type)) +
    geom_boxplot(aes(fill = tissue_type),outlier.shape = NA,alpha = 0.2) +
    geom_jitter_rast(shape = 16, position = position_jitter(0.2),alpha = 0.5,size = 2) +
    scale_color_manual(values = tissue_type_color) +
    scale_fill_manual(values = tissue_type_color) +
    nrc_theme + NoLegend() + xlab("")
  function_plot(filename_prefix = str_glue("{output_dir}/CM_diff_tissue/scRNA_{CM_type}_diff_tissue_percentage"), width = 300, height = 120)
}


# CM in tumor code --------------------------------------------------------
for(CM_type in as.character(unique(CM_status_stat_final$CM))){
  # CM_type <- "CM1"
  
  ## data filter
  current_meta_stat <- CM_status_stat_final %>% 
    filter(CM == CM_type) %>% 
    filter(tumor_status != "Healthy")
  
  ## data plot
  ggplot(current_meta_stat, aes(tumor_code,percentage,color = tumor_code)) +
    geom_boxplot(aes(fill = tumor_code),outlier.shape = NA,alpha = 0.2) +
    geom_jitter_rast(shape = 16, position = position_jitter(0.2),alpha = 0.5,size = 2) +
    scale_x_discrete(limits = c("ANS","BRCA","CESC","CRC","ESCC","GC","HNSC","KIRC","LSCC","NSCLC","LUAD","LUSC","MA","NPC","PAAD","PNET",
                                "SARC","TGCT","THCA","ACC")) +
    scale_color_manual(values = tumor_type_color) +
    scale_fill_manual(values = tumor_type_color) +
    nrc_theme + NoLegend() + xlab("")
  function_plot(filename_prefix = str_glue("{output_dir}/CM_diff_tumor/scRNA_{CM_type}_diff_tissue_percentage"), width = 260, height = 120)
}





