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
library("survival")
library("survminer")
# library("scPDtools")
library("ComplexHeatmap")

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
bulkRNA_dir <- " "
scRNA_dir <- " "
output_dir <- " "


# data color --------------------------------------------------------------
ano_colors <- list(
  CM = c("CM1"="#ea5c6f","CM2"="#f7905a","CM3"="#3187cb","CM4"="#fb948d","CM5"="#e2b159",
         "CM6"="#ebed6f","CM7"="#b2db87","CM8"="#7ee7bb","CM9"="#64cccf"
  ),
  tumor_code = c("ACC"="#5A7B8F","PRAD"="#EE4C97","BRCA"="#3D806F","CESC"="#A08634","CRC"="#F37C95","ESCA"="#608541",
                 "HNSC"="#7D4E57","KIRC"="#BC3C29","LUAD"="#958056","LUSC"="#9FAFA3","PCPG"="#6F99AD","PAAD"="#0072B5",
                 "READ"="#CFC59A","SARC"="#E18727","SKCM"="#FFDC91","GC"="#718DAE","TGCT"="#7876B1","THCA"="#F9AC93",
                 "NSCLC"="#E3BC06","THYM"="#2285FE","UVM"="#3E6086","UCEC"="#20854E","UCS"="#4A7985","GBM"="#9F4700",
                 "OV"="#55d151","MESO"="#f69a95","LIHC"="#aa329a","KIRP"="#EDE816","KICH"="#48dc88","DLBC"="#48dc88",
                 "CHOL"="#9ebe71","LGG"="#FDC44D","LAML"="#16A795","STAD"="#6FCDDC","COAD"="#E7ABFD","BLCA"="#CEE769"),
  tissue_type = c("tumor" = "#5BC0EB","normal" = "#9BC53D")
  
)
tumor_code = c("ACC"="#5A7B8F","PRAD"="#EE4C97","BRCA"="#3D806F","CESC"="#A08634","CRC"="#F37C95","ESCA"="#608541",
               "HNSC"="#7D4E57","KIRC"="#BC3C29","LUAD"="#958056","LUSC"="#9FAFA3","PCPG"="#6F99AD","PAAD"="#0072B5",
               "READ"="#CFC59A","SARC"="#E18727","SKCM"="#FFDC91","GC"="#718DAE","TGCT"="#7876B1","THCA"="#F9AC93",
               "NSCLC"="#E3BC06","THYM"="#2285FE","UVM"="#3E6086","UCEC"="#20854E","UCS"="#4A7985","GBM"="#9F4700",
               "OV"="#55d151","MESO"="#f69a95","LIHC"="#aa329a","KIRP"="#EDE816","KICH"="#48dc88","DLBC"="#48dc88",
               "CHOL"="#9ebe71","LGG"="#FDC44D","LAML"="#16A795","STAD"="#6FCDDC","COAD"="#E7ABFD","BLCA"="#CEE769")

all_cox_result <- data.frame(
  MP_name = as.character(),
  tumor_code = as.character(),
  HR = as.numeric(),
  L95CI = as.numeric(),
  H95CI = as.numeric(),
  pvalue = as.numeric()
)

# data read ---------------------------------------------------------------
## DEG read
module_gene_list <- readRDS(str_glue("{scRNA_dir}/all_CM_DEG.RData"))

## bulkRNA fpkm read
bulkRNA_fpkm <- read_csv(str_glue("{bulkRNA_dir}/bulkRNA_fpkm_filter.csv"))

## bulkRNA sample read
bulkRNA_sample <- read_csv(str_glue("{bulkRNA_dir}/bulkRNA_smaple_filter.csv"))

## bulkRNA survival read
bulkRNA_survival <- read_csv(str_glue("{bulkRNA_dir}/bulkRNA_survival_filter.csv"))


# data filter -------------------------------------------------------------
bulkRNA_sample_filter <- bulkRNA_sample
bulkRNA_sample_filter <- bulkRNA_sample_filter[!duplicated(bulkRNA_sample_filter$sample),]
  # filter(study == "TCGA" & tissue_type == "tumor")
bulkRNA_sample_tumor <- bulkRNA_sample_filter %>% 
  filter(study == "TCGA" & tissue_type == "tumor")
bulkRNA_fpkm_filter <- bulkRNA_fpkm %>% 
  dplyr::select(c("gene_name",bulkRNA_sample_filter$sample))


# bulkRNA remove duplicate gene -------------------------------------------
bulkRNA_fpkm_filter <- bulkRNA_fpkm_filter[!duplicated(bulkRNA_fpkm_filter$gene_name),]
bulkRNA_fpkm_filter <- column_to_rownames(bulkRNA_fpkm_filter,var = "gene_name")

# GSVA --------------------------------------------------------------------
## run gsva
ssGSEA_matrix <- gsva(expr = as.matrix(bulkRNA_fpkm_filter), 
                      gset.idx.list = module_gene_list,
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
colnames(ano_row) <- "CM"
ano_col <- as.data.frame(colnames(sorted_mat))
colnames(ano_col) <- "sample"
ano_col <- left_join(ano_col,bulkRNA_sample_filter[,c(1,8,9)],by = "sample")
rownames(ano_col) <- bulkRNA_sample_filter$sample
ano_col$tumor_code <- factor(ano_col$tumor_code)
ano_col <- ano_col[-1]

pdf(str_glue("{output_dir}/04.bulkRNA_CM_GSVA/bulkRNA_GSVA.pdf"), width = 19, height = 6)
pheatmap(sorted_mat,
         cluster_row = F, cluster_cols = F, show_rownames = T, 
         annotation_row = ano_row,
         annotation_col = ano_col,
         show_colnames = F, legend = T, annotation_names_row = F, annotation_names_col = F,
         annotation_colors = ano_colors,
         cellwidth = 0.05,cellheight = 20,
         fontsize = 5,
         use_raster = F
)
dev.off()



# survival analysis -------------------------------------------------------
## data combine
GSVA_data <- as.data.frame(sorted_mat) %>% 
  dplyr::select(bulkRNA_sample_tumor$sample)
GSVA_data_transform <- GSVA_data %>% 
  t(.) %>% 
  as.data.frame(.) %>% 
  rownames_to_column(var = "sample") %>% 
  left_join(.,bulkRNA_survival,by = "sample") %>% 
  left_join(.,bulkRNA_sample[,c(1,8)],by = "sample")

## survival
status_cutoff <- 0.1
for(module in c("CM1","CM2","CM3","CM4","CM5","CM6","CM7","CM8","CM9")){
  # tumor <- "BRCA"
  
  for(tumor in unique(GSVA_data_transform$tumor_code)){
    # module <- "Cell.Cycle.G2M"
    if(tumor == "LAML"){
      next
    }
    
    if(tumor == "PCPG" & module == "Cell.Signaling"){
      next
    }
    if(tumor == "CHOL" & module == "Coagulation"){
      next
    }
    
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

## result stat
all_cox_result_transform <- all_cox_result %>% 
  mutate(survival_status = ifelse(HR > 1 & pvalue <= 0.05, "good", 
                                  ifelse(HR < 1 & pvalue <= 0.05, "poor", "other")))
all_cox_result_transform <- all_cox_result_transform[!is.na(all_cox_result_transform$pvalue),]

## pheatmap plot
ggplot(all_cox_result_transform, aes(x = tumor_code, y = MP_name, fill = survival_status)) +
  geom_tile(color = "white") +
  # scale_fill_gradient2(low = "#6D9EC1",mid = '#D9D9D9', high = "#E46726") +
  theme_minimal()+
  # scale_y_discrete(limits = colnames(MP_3CA_hub)) +
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
function_plot(filename_prefix = str_glue("{output_dir}/05.bulkRNA_CM_survival/MP_cor"), width = 420, height = 240)

# hub MP per --------------------------------------------------------------
bulkRNA_gsva_sorted <- as.data.frame(sorted_mat)
bulkRNA_gsva_sorted <- rownames_to_column(bulkRNA_gsva_sorted,var = "CM")
# bulkRNA_gsva_sorted <- read_csv(str_glue("{output_dir}/bulkRNA_gsva_sorted.csv"))

## data transform
TCGA_plot_data <- bulkRNA_gsva_sorted %>% 
  column_to_rownames(var = "final_MP") %>% 
  t(.) %>% 
  as.data.frame(.) %>% 
  rownames_to_column(var = "sample") %>% 
  left_join(bulkRNA_sample_filter[,c("sample","tumor_code")],by = "sample")

## data plot
for(MP in names(module_gene_list)){
  # MP <- "Cell.Cycle.G2M"
  
  current_data <- TCGA_plot_data %>% 
    dplyr::select(c("sample",MP,"tumor_code"))
  colnames(current_data)[2] <- "score"
  current_data <- current_data %>% 
    mutate(tumor_code = fct_reorder(tumor_code,score,.desc = T,.fun = median))
  
  ggviolin(current_data, 
           x="tumor_code", y="score", 
           width = 1.5,color = "black",
           fill="tumor_code",
           xlab = F,
           add = "boxplot",
           add.params = list(outlier.shape = NA),
           bxp.errorbar=T, 
           bxp.errorbar.width=0.02,
           size=0.5,
           palette = tumor_code, 
           legend = NULL) + nrc_theme + NoLegend()
  function_plot(filename_prefix = str_glue("{output_dir}/GSVA_{MP}_ggviolin"), width = 450, height = 150)
}


# data save ---------------------------------------------------------------


write_csv(bulkRNA_gsva_sorted,str_glue("{output_dir}/bulkRNA_gsva_sorted.csv"))
