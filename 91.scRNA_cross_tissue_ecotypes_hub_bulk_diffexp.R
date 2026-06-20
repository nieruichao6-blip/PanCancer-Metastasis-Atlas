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
scRNA_data_dir <- " "
ET_ano_dir <- " "
output_dir <- " "
all_sample_exp <- data.frame(
  gene_id = as.character()
)
all_DEG_result <- data.frame(
  gene_id = as.character(),
  FDR = as.numeric(),
  Pvalue = as.numeric(),
  logFC = as.numeric(),
  Ecotype = as.character()
)

# data read ---------------------------------------------------------------
ET_ano <- read_csv(str_glue("{ET_ano_dir}/cluster_ano.csv"))


# scRNA transform bulkRNA -------------------------------------------------
scRNA_lsit <- list.files(scRNA_data_dir)
for(tumor_type in scRNA_lsit){
  # tumor_type <- "ANS"
  print(str_glue("################################## {tumor_type} ######################################"))
  
  ## data read
  current_data <- readRDS(str_glue("{scRNA_data_dir}/{tumor_type}/scRNA_annotation_new.RData"))
  Idents(current_data) <- current_data$sample_ID
  
  ## data normal
  current_data <- NormalizeData(current_data,normalization.method = 'LogNormalize',scale.factor = 10000)
  
  ## data transform
  ave_exp <- as.matrix(AverageExpression(current_data)[[1]])
  ave_exp <- as.data.frame(ave_exp) %>% 
    rownames_to_column(var = "gene_id")
  
  ## all sample data collect
  all_sample_exp <- full_join(all_sample_exp,ave_exp,by = "gene_id")
  
}


# all exp data transform --------------------------------------------------
## add NA data
all_sample_exp_transform <- all_sample_exp
all_sample_exp_transform[is.na(all_sample_exp_transform)] <- 0

## keep hub gene
all_sample_exp_transform <- all_sample_exp_transform[!grepl("^ENSG",all_sample_exp_transform$gene_id),]
all_sample_exp_transform <- all_sample_exp_transform[!grepl("^LINC",all_sample_exp_transform$gene_id),]

## hub sample filter
all_sample_exp_hub <- all_sample_exp_transform[,c("gene_id",ET_ano$sample_ID)]

## remove low exp gene
expr_count <- rowSums(all_sample_exp_hub[,-1] > 0,na.rm = T)
keep <- expr_count >= (542/5)
all_sample_exp_hub <- all_sample_exp_hub[keep,]

## exp save
saveRDS(all_sample_exp,str_glue("{output_dir}/sample_origin_exp.RData"))


# diffexp analysis --------------------------------------------------------
for(ET_type in unique(ET_ano$Ecotype)){
  # ET_type <- "ET_01"
  print(str_glue("################################## {ET_type} ######################################"))
  
  ## group confirm
  treating_label <- "hub_ET"
  treating_sample_list <- ET_ano[ET_ano$Ecotype == ET_type,]$sample_ID
  ctl_label <- "other_ET"
  ctl_sample_list <- ET_ano[ET_ano$Ecotype != ET_type,]$sample_ID
  group <- factor(c(rep(ctl_label,time = length(ctl_sample_list)),
                    rep(treating_label,time = length(treating_sample_list))),
                  levels = c(ctl_label,treating_label))
  
  ## exp data transform
  all_sample_exp_limma <- all_sample_exp_hub %>% 
    dplyr::group_by(gene_id) %>% 
    dplyr::summarise(across(where(is.numeric),mean,na.rm = T))
    
  all_sample_exp_limma <- all_sample_exp_limma %>% 
    column_to_rownames(var = "gene_id") %>% 
    dplyr::select(ctl_sample_list,treating_sample_list)
  
  ## limma run
  design <- model.matrix(~0+factor(group))
  colnames(design) <- levels(factor(group))
  rownames(design) <- colnames(all_sample_exp_limma)
  
  contrast.matrix <- makeContrasts(hub_ET-other_ET,levels = design)
  
  fit <- lmFit(all_sample_exp_limma,design)
  fit2 <- contrasts.fit(fit,contrast.matrix)
  fit2 <- eBayes(fit2)
  DEG <- topTable(fit2,coef = 1,n = Inf,sort.by = "logFC")
  DEG <- na.omit(DEG)
  DEG <- subset(DEG,select = c("adj.P.Val","P.Value","logFC"))
  colnames(DEG) <- c("FDR","Pvalue","logFC")
  
  ## result transform
  DEG <- rownames_to_column(DEG,var = "gene_id")
  # DEG_filter <- DEG[DEG$FDR <= 0.0001 & DEG$logFC >= log2(1.2),]
  DEG_filter <- DEG
  DEG_filter$Ecotype <- ET_type
  
  ## all DEG collect
  all_DEG_result <- rbind(all_DEG_result,DEG_filter)
  
}



# all DEG filter ----------------------------------------------------------
all_DEG_result_filter <- all_DEG_result %>% 
  filter(FDR <= 0.001 & logFC > log2(1)) %>% 
  group_by(Ecotype) %>%
  arrange(desc(logFC)) %>%
  slice_head(n = 50) %>%
  ungroup()

## data save
write_csv(all_DEG_result,str_glue("{output_dir}/all_DEG_result.csv"))
write_csv(all_DEG_result_filter,str_glue("{output_dir}/all_DEG_result_filter.csv"))












