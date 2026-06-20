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
library("survival")
library("survminer")
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
ano_colors <- list(
  final_MP = Ecotype_color,
  tumor_code = c("ACC"="#5A7B8F","PRAD"="#EE4C97","BRCA"="#3D806F","CESC"="#A08634","CRC"="#F37C95","ESCA"="#608541",
                 "HNSC"="#7D4E57","KIRC"="#BC3C29","LUAD"="#958056","LUSC"="#9FAFA3","PCPG"="#6F99AD","PAAD"="#0072B5",
                 "READ"="#CFC59A","SARC"="#E18727","SKCM"="#FFDC91","GC"="#718DAE","TGCT"="#7876B1","THCA"="#F9AC93",
                 "NSCLC"="#E3BC06","THYM"="#2285FE","UVM"="#3E6086","UCEC"="#20854E","UCS"="#4A7985","GBM"="#9F4700",
                 "OV"="#55d151","MESO"="#f69a95","LIHC"="#aa329a","KIRP"="#EDE816","KICH"="#48dc88","DLBC"="#48dc88",
                 "CHOL"="#9ebe71","LGG"="#FDC44D","LAML"="#16A795","STAD"="#6FCDDC","COAD"="#E7ABFD","BLCA"="#CEE769")
  
)

all_cox_result <- data.frame(
  MP_name = as.character(),
  tumor_code = as.character(),
  HR = as.numeric(),
  L95CI = as.numeric(),
  H95CI = as.numeric(),
  pvalue = as.numeric()
)

# data dir ----------------------------------------------------------------
bulkRNA_dir <- " "
ET_hub_gene_dir <- " "
output_dir <- " "
ET_hub_gene_list <- list()

# data read ---------------------------------------------------------------
bulkRNA_fpkm <- read_csv(str_glue("{bulkRNA_dir}/bulkRNA_fpkm_filter.csv"))
bulkRNA_sample <- read_csv(str_glue("{bulkRNA_dir}/bulkRNA_smaple_filter.csv"))
bulkRNA_survival <- read_csv(str_glue("{bulkRNA_dir}/bulkRNA_survival_filter.csv"))
ET_hub_gene <- read_csv(str_glue("{ET_hub_gene_dir}/all_DEG_result.csv"))


# data transform ----------------------------------------------------------
## bulk data
bulkRNA_sample_filter <- bulkRNA_sample %>% 
  filter(study == "TCGA" & tissue_type == "tumor")
bulkRNA_fpkm_filter <- bulkRNA_fpkm %>% 
  dplyr::select(c("gene_name",bulkRNA_sample_filter$sample))
bulkRNA_fpkm_filter <- bulkRNA_fpkm_filter[!duplicated(bulkRNA_fpkm_filter$gene_name),]
bulkRNA_fpkm_filter <- column_to_rownames(bulkRNA_fpkm_filter,var = "gene_name")

## ET hub gene
ET_hub_gene <- ET_hub_gene[!grepl("^RP",ET_hub_gene$gene_id),]
ET_hub_gene <- ET_hub_gene[!grepl("^MT-",ET_hub_gene$gene_id),]
ET_hub_gene_filter <- ET_hub_gene %>% 
  filter(FDR <= 0.001 & logFC > log2(1)) %>% 
  group_by(Ecotype) %>%
  arrange(desc(logFC)) %>%
  slice_head(n = 50) %>%
  ungroup()
for(ET_type in unique(ET_hub_gene_filter$Ecotype)){
  # ET_type <- "ET_01"
  
  current_ET_gene <- ET_hub_gene_filter[ET_hub_gene_filter$Ecotype == ET_type,]$gene_id
  ET_hub_gene_list[[ET_type]] <- current_ET_gene
  
}


# GSVA analysis -----------------------------------------------------------
## run gsva
ssGSEA_matrix <- gsva(expr = as.matrix(bulkRNA_fpkm_filter), 
                      gset.idx.list = ET_hub_gene_list,
                      method = 'ssgsea',kcdf = 'Gaussian',abs.ranking = TRUE)

## data normal
scale_gsva_matrix <- t(scale(t(ssGSEA_matrix)))
scale_gsva_matrix[scale_gsva_matrix< -2] <- -2
scale_gsva_matrix[scale_gsva_matrix>2] <- 2
normalization <- function(x){return((x-min(x))/(max(x)-min(x)))} 
nor_gsva_matrix <- normalization(scale_gsva_matrix)

## new sorted
if(T){
  dominant_clusters <- apply(nor_gsva_matrix, 2, function(x) rownames(nor_gsva_matrix)[which.max(x)])
  ordered_samples <- character(0)
  for (cluster in rownames(nor_gsva_matrix)) {
    candidate_samples <- names(dominant_clusters)[dominant_clusters == cluster]
    candidate_samples <- setdiff(candidate_samples, ordered_samples)
    samples_sorted <- names(sort(nor_gsva_matrix[cluster, candidate_samples], decreasing = TRUE))
    
    ordered_samples <- c(ordered_samples, samples_sorted)
  }
  sorted_mat <- nor_gsva_matrix[, ordered_samples]
  
}

## result plot
ano_row <-  as.data.frame(rownames(sorted_mat))
rownames(ano_row) <- rownames(sorted_mat)
colnames(ano_row) <- "final_MP"
ano_col <- as.data.frame(colnames(sorted_mat))
colnames(ano_col) <- "sample"
ano_col <- left_join(ano_col,bulkRNA_sample_filter[,c(1,8)],by = "sample")
rownames(ano_col) <- bulkRNA_sample_filter$sample
ano_col$tumor_code <- factor(ano_col$tumor_code)
ano_col <- ano_col[-1]

pdf(str_glue("{output_dir}/bulkRNA_GSVA_ET.pdf"), width = 10, height = 6)
pheatmap(sorted_mat,
         cluster_row = F, cluster_cols = F, show_rownames = T, 
         annotation_row = ano_row,
         annotation_col = ano_col,
         show_colnames = F, legend = T, annotation_names_row = F, annotation_names_col = F,
         annotation_colors = ano_colors,
         cellwidth = 0.05,cellheight = 14,
         fontsize = 5,
         use_raster = F
)
dev.off()

bulkRNA_gsva_sorted <- as.data.frame(sorted_mat)
bulkRNA_gsva_sorted <- rownames_to_column(bulkRNA_gsva_sorted,var = "Ecotype")

# GSVA analysis -----------------------------------------------------------
GSVA_data_transform <- bulkRNA_gsva_sorted %>% 
  column_to_rownames(var = "Ecotype") %>% 
  t(.) %>% 
  as.data.frame(.) %>% 
  rownames_to_column(var = "sample") %>% 
  left_join(.,bulkRNA_survival,by = "sample") %>% 
  left_join(.,bulkRNA_sample[,c(1,8)],by = "sample")

## survival
status_cutoff <- 0.1
for(module in unique(bulkRNA_gsva_sorted$Ecotype)){
  # tumor <- "BRCA"
  
  for(tumor in unique(GSVA_data_transform$tumor_code)){
    # module <- "Cell.Cycle.G2M"
    if(tumor == "LAML"){
      next
    }
    
    # if(tumor == "DLBC" & module == "ET_05"){
    #   next
    # }
    # if(tumor == "CHOL" & module == "Coagulation"){
    #   next
    # }
    
    current_GSVA <- GSVA_data_transform %>% 
      filter(tumor_code == tumor)
    
    
    ## best cutoff find
    best_cutoff <- surv_cutpoint(current_GSVA,time = "survival_time", event = "survival_status", variables = module, minprop = status_cutoff)
    current_GSVA <- surv_categorize(best_cutoff)
    colnames(current_GSVA)[3] <- "status"
    
    # ## data transform
    # high_cutoff <- quantile(current_GSVA[[module]], status_cutoff)
    # low_cutoff <- quantile(current_GSVA[[module]], 1-status_cutoff)
    # # high_cutoff <- median(current_GSVA[[module]])
    # # low_cutoff <- median(current_GSVA[[module]])
    # current_GSVA <- current_GSVA %>%
    #   mutate(status = ifelse(get(module) >= high_cutoff,"high",
    #                          ifelse(get(module) < low_cutoff,"low","None"))) %>%
    #   filter(status != "None") %>%
    #   arrange(status)
    if(length(unique(current_GSVA$status)) < 2){
      next
    }
    
    ### survival
    surv_formula <- as.formula('Surv(survival_time, survival_status)~status')
    # surv_formula <- as.formula('Surv(PFI.time, PFI)~status')
    current_cox <- coxph(surv_formula,data = current_GSVA)
    current_cox_result <- summary(current_cox)
    
    ### result combine
    current_result <- data.frame(
      MP_name = module,
      tumor_code = tumor,
      HR = current_cox_result$coefficients[,2],
      L95CI = current_cox_result$conf.int[,3],
      H95CI = current_cox_result$conf.int[,4],
      pvalue = current_cox_result$coefficients[,5]
    )
    all_cox_result <- rbind(all_cox_result,current_result)
    
  }
}

## result plot
all_cox_result_transform <- all_cox_result %>% 
  mutate(survival_status = ifelse(HR > 1 & pvalue <= 0.05, "good", 
                                  ifelse(HR < 1 & pvalue <= 0.05, "poor", "other")))
all_cox_result_transform <- all_cox_result_transform[!is.na(all_cox_result_transform$pvalue),]

## pheatmap plot
ggplot(all_cox_result_transform, aes(x = tumor_code, y = MP_name, fill = survival_status)) +
  geom_tile(color = "white") +
  theme_minimal()+
  theme(plot.title = element_text(face = "bold",size = 12,color="black",hjust = 0.5),
        axis.title = element_text(face = "bold",size = 12,color ="black"), 
        axis.text = element_text(face = "bold",size= 16,color = "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1 ),
        panel.grid=element_blank(),
        legend.text = element_text(face = "bold",size= 12),
        legend.title= element_text(face = "bold",size= 12)
  ) + labs(fill = "pearson") + xlab("") + ylab("")
function_plot(filename_prefix = str_glue("{output_dir}/ET_cor"), width = 420, height = 200)

## survival poor stat
survival_poor_data <- all_cox_result_transform %>% 
  filter(survival_status == "poor")
survival_poor_stat <- as.data.frame(table(survival_poor_data$MP_name))
colnames(survival_poor_stat) <- c("MP_name","num")
ggplot(survival_poor_stat, aes(MP_name ,num, fill = MP_name)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "", y = "Number of survival num") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = Ecotype_color) +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/ET_survival_poor_stat"), width = 420, height = 250)

## survival good stat
survival_poor_data <- all_cox_result_transform %>% 
  filter(survival_status == "good")
survival_poor_stat <- as.data.frame(table(survival_poor_data$MP_name))
colnames(survival_poor_stat) <- c("MP_name","num")
ggplot(survival_poor_stat, aes(MP_name ,num, fill = MP_name)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "", y = "Number of survival num") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = Ecotype_color) +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/ET_survival_good_stat"), width = 420, height = 250)
