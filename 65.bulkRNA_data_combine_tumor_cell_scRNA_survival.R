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
library("scPDtools")
library("ComplexHeatmap")
library("survival")
library("survminer")
(.packages())
setwd(' ')
display.brewer.all()

function_plot <- function(filename_prefix, width, height){
  ggsave(filename = str_c(filename_prefix, '.pdf'), plot = last_plot(), device = 'pdf', width = width, height = height, units = 'mm', dpi = 500)
  ggsave(filename = str_c(filename_prefix, '.svg'), plot = last_plot(), device = 'svg', width = width, height = height, units = 'mm', dpi = 500)
  ggsave(filename = str_c(filename_prefix, '.png'), plot = last_plot(), device = 'png', width = width, height = height, units = 'mm', dpi = 1000)
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
                       axis.text.x = element_text(angle = 45, hjust = 0.5 ),
                       panel.grid=element_blank(),
                       legend.text = element_text(face = "bold",size= 12),
                       legend.title= element_text(face = "bold",size= 12))
survival_plot_theme <- theme(
  panel.background = element_rect(fill = "transparent"),
  panel.border = element_rect(color = "black", size = 2, linetype = "solid", inherit.blank = TRUE, fill='transparent'),
  axis.title.y.left = element_text(size = 20, face = 'bold', vjust=3),
  axis.title.x.bottom = element_text(size = 20, face = 'bold', vjust=-2),
  axis.text = element_text(color="black", size = 16, face = 'bold'),
  axis.ticks = element_line(size = 1), axis.ticks.length = unit("2", "mm"),
  legend.text = element_text(color="black", size = 16, face = 'bold'),
  legend.title = element_text(color="black", size = 16, face = 'bold'),
  plot.margin = unit(c(3, 3, 4, 4),'mm')
)


# data dir ----------------------------------------------------------------
bulkRNA_dir <- " "
GSVA_dir <- " "
output_dir <- " "
module_type <- c("Cell.Cycle.G2M","Cell.Cycle.G1S","Cell.Cycle.HMG","Chromatin","Cell.Cycle.single.nucleus","Stress","Hypoxia",
                 "Stress.in.vitro","Protein.maturation","EMT.I","EMT.II","EMT.VI","MHC","Epithelial.Senescence","MYC",
                 "Respriration","Secreted.I","Secreted.II","Cilia.I","Cilia.II","PDAC.classcial","Alveolar","Colon.related",
                 "PDAC.related","Lymphocyte.Activation","Cell.Signaling","Coagulation","Neuron")

all_cox_result <- data.frame(
  MP_name = as.character(),
  tumor_code = as.character(),
  HR = as.numeric(),
  L95CI = as.numeric(),
  H95CI = as.numeric(),
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
)
survival_status <- c("good"="","poor"="","other"="")



# data read ---------------------------------------------------------------
bulkRNA_survival <- read_csv(str_glue("{bulkRNA_dir}/bulkRNA_survival_filter.csv"))
bulkRNA_sample <- read_csv(str_glue("{bulkRNA_dir}/bulkRNA_smaple_filter.csv"))
GSVA_data <- read_csv(str_glue("{GSVA_dir}/bulkRNA_gsva_sorted.csv"))


# data combine ------------------------------------------------------------
GSVA_data_transform <- GSVA_data %>% 
  column_to_rownames(var = "final_MP") %>% 
  t(.) %>% 
  as.data.frame(.) %>% 
  rownames_to_column(var = "sample") %>% 
  left_join(.,bulkRNA_survival,by = "sample") %>% 
  left_join(.,bulkRNA_sample[,c(1,8)],by = "sample")


# survival analysis -------------------------------------------------------
status_cutoff <- 0.1
for(module in module_type){
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


# result plot -------------------------------------------------------------
all_cox_result_transform <- all_cox_result %>% 
  mutate(survival_status = ifelse(HR > 1 & pvalue <= 0.05, "good", 
                                  ifelse(HR < 1 & pvalue <= 0.05, "poor", "other")))
all_cox_result_transform <- all_cox_result_transform[!is.na(all_cox_result_transform$pvalue),]

write_csv(all_cox_result_transform,str_glue("{output_dir}/all_cox_result_transform.csv"))

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
function_plot(filename_prefix = str_glue("{output_dir}/MP_cor"), width = 420, height = 350)

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
  scale_fill_manual(values = final_MP) +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/MP_survival_poor_stat"), width = 420, height = 250)
