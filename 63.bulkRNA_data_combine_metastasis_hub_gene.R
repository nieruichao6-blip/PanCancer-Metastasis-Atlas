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
# library("scPDtools")
library("ComplexHeatmap")
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

all_cox_result <- data.frame(
  tumor_code = as.character(),
  HR = as.numeric(),
  L95CI = as.numeric(),
  H95CI = as.numeric(),
  pvalue = as.numeric()
)

tumor_code_color <- c("ACC"="#5A7B8F","PRAD"="#EE4C97","BRCA"="#3D806F","CESC"="#A08634","CRC"="#F37C95","ESCA"="#608541",
               "HNSC"="#7D4E57","KIRC"="#BC3C29","LUAD"="#958056","LUSC"="#9FAFA3","PCPG"="#6F99AD","PAAD"="#0072B5",
               "READ"="#CFC59A","SARC"="#E18727","SKCM"="#FFDC91","GC"="#718DAE","TGCT"="#7876B1","THCA"="#F9AC93",
               "NSCLC"="#E3BC06","THYM"="#2285FE","UVM"="#3E6086","UCEC"="#20854E","UCS"="#4A7985","GBM"="#9F4700",
               "OV"="#55d151","MESO"="#f69a95","LIHC"="#aa329a","KIRP"="#EDE816","KICH"="#48dc88","DLBC"="#48dc88",
               "CHOL"="#9ebe71","LGG"="#FDC44D","LAML"="#16A795","STAD"="#6FCDDC","COAD"="#E7ABFD","BLCA"="#CEE769")

# data dir ----------------------------------------------------------------
bulkRNA_dir <- " "
scRNA_dir <- " "
output_dir <- " "


# data read ---------------------------------------------------------------
## bulkRNA fpkm read
bulkRNA_fpkm <- read_csv(str_glue("{bulkRNA_dir}/bulkRNA_fpkm_filter.csv"))

## bulkRNA sample read
bulkRNA_sample <- read_csv(str_glue("{bulkRNA_dir}/bulkRNA_smaple_filter.csv"))

## bulkRNA survival read
bulkRNA_survival <- read_csv(str_glue("{bulkRNA_dir}/bulkRNA_survival_filter.csv"))

## metastasis hub gene read
metastasis_gene <- read_csv(str_glue("{scRNA_dir}/high_per_DEG_all_metastasis.csv"))
colnames(metastasis_gene) <- "gene_name"
metastasis_gene_list <- list(metastasis_hub_gene = metastasis_gene$gene_name)

# data filter -------------------------------------------------------------
bulkRNA_sample_filter <- bulkRNA_sample %>% 
  filter(study == "TCGA" & tissue_type == "tumor")
bulkRNA_fpkm_filter <- bulkRNA_fpkm %>% 
  dplyr::select(c("gene_name",bulkRNA_sample_filter$sample))


# bulkRNA remove duplicate gene -------------------------------------------
bulkRNA_fpkm_filter <- bulkRNA_fpkm_filter[!duplicated(bulkRNA_fpkm_filter$gene_name),]
bulkRNA_fpkm_filter <- column_to_rownames(bulkRNA_fpkm_filter,var = "gene_name")

# GSVA --------------------------------------------------------------------
## run gsva
ssGSEA_matrix <- gsva(expr = as.matrix(bulkRNA_fpkm_filter), 
                      gset.idx.list = metastasis_gene_list,
                      method = 'ssgsea',kcdf = 'Gaussian',abs.ranking = TRUE)


# result combine survival -------------------------------------------------
GSVA_data_transform <- ssGSEA_matrix %>% 
  # column_to_rownames(var = "metastasis_hub_gene") %>% 
  t(.) %>% 
  as.data.frame(.) %>% 
  rownames_to_column(var = "sample") %>% 
  left_join(.,bulkRNA_survival,by = "sample") %>% 
  left_join(.,bulkRNA_sample[,c(1,8)],by = "sample")


# survival analysis -------------------------------------------------------
status_cutoff <- 0.1
for(tumor in unique(GSVA_data_transform$tumor_code)){
  # tumor <- "BRCA"
  # if(tumor == "LAML"){
  #   next
  # }
  # 
  # if(tumor == "PCPG" & module == "Cell.Signaling"){
  #   next
  # }
  # if(tumor == "CHOL" & module == "Coagulation"){
  #   next
  # }
  
  current_GSVA <- GSVA_data_transform %>% 
    filter(tumor_code == tumor)
  
  
  ## best cutoff find
  best_cutoff <- surv_cutpoint(current_GSVA,time = "survival_time", event = "survival_status", variables = "metastasis_hub_gene", minprop = status_cutoff)
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
    # MP_name = module,
    tumor_code = tumor,
    HR = current_cox_result$coefficients[,2],
    L95CI = current_cox_result$conf.int[,3],
    H95CI = current_cox_result$conf.int[,4],
    pvalue = current_cox_result$coefficients[,5]
  )
  all_cox_result <- rbind(all_cox_result,current_result)
  
}

# result plot -------------------------------------------------------------
all_cox_result_transform <- all_cox_result %>% 
  mutate(survival_status = ifelse(HR > 1 & pvalue <= 0.05, "good", 
                                  ifelse(HR < 1 & pvalue <= 0.05, "poor", "other")))
all_cox_result_transform <- all_cox_result_transform[!all_cox_result_transform$tumor_code %in% c("THYM","TGCT","CHOL","DLBC"),]
  # filter(survival_status == "poor")

ggplot(all_cox_result_transform, aes(HR, tumor_code))+
  geom_point(size=4,aes(color=tumor_code))+
  geom_errorbarh(aes(xmax = H95CI, xmin = L95CI,color=tumor_code),size= 1,height = 0) +
  scale_shape_manual(values = c("circle"))+
  scale_color_manual(values = tumor_code_color)+
  scale_x_continuous(limits= c(0, 4))+
  geom_vline(aes(xintercept = 1),color="gray",linetype="dashed", size = 0.8) +
  # geom_text(aes(label=tumor_code,x=-15+3))+
  xlab('ln RR (%) ')+
  ylab(' ')+
  # theme_few()+
  theme(axis.text.x = element_text(size = 24, color = "black"))+
  theme(axis.text.y = element_text(size = 24, color = "black"))+
  theme(title=element_text(size=24))+
  theme_bw()+theme(panel.grid=element_blank())+
  theme(axis.ticks.length=unit(0.1, "cm"),
        axis.text.x = element_text(margin=unit(c(0.1,0.1,0.1,0.1), "cm")),
        axis.text.y = element_text(margin=unit(c(0.1,0.1,0.1,0.1), "cm")) )+theme(legend.position = "none")
function_plot(filename_prefix = str_glue("{output_dir}/bulkRNA_metastasis_hub_gene_survival"), width = 150, height = 300)