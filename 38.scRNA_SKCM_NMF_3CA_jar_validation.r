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
(.packages())
setwd(' ')
display.brewer.all()

function_plot <- function(filename_prefix, width, height){
  ggsave(filename = str_c(filename_prefix, '.pdf'), plot = last_plot(), device = 'pdf', width = width, height = height, units = 'mm', dpi = 500)
  ggsave(filename = str_c(filename_prefix, '.svg'), plot = last_plot(), device = 'svg', width = width, height = height, units = 'mm', dpi = 500)
  ggsave(filename = str_c(filename_prefix, '.png'), plot = last_plot(), device = 'png', width = width, height = height, units = 'mm', dpi = 500)
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
  theme(axis.text.x = element_text(angle = 90, hjust = 1 )) +
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
jaccard_result <- data.frame(
  MP_3CA_name = as.character(),
  MP_self_name = as.character(),
  jaccard_similarity_score = as.numeric()
)

# data dir ----------------------------------------------------------------
MP_3CA_hub_dir <- " "
input_dir <- " "
output_dir <- " "

# data read ----------------------------------------------------------------
MP_3CA_hub <- read.xlsx(str_glue("{MP_3CA_hub_dir}"),sheetName = "cancer_MPs")
meta_program_top50 <- readRDS(str_glue("{input_dir}/meta_program_top50.RData"))
# scRNA_tumor <- readRDS(str_glue("{input_dir}/scRNA_tumor_NMF.RData"))

# jaccard cor ----------------------------------------------------------------
for(MP_3CA in colnames(MP_3CA_hub)){
  MP_3CA_gene <- head(MP_3CA_hub[[MP_3CA]],50)
  
  for(MP_slef in names(meta_program_top50)){
    
    MP_slef_gene <- head(meta_program_top50[[MP_slef]],50)
    
    ## jaccard
    similarity_score <- jaccard(MP_3CA_gene,MP_slef_gene)
    
    ## result collect
    current_result <- data.frame(
      MP_3CA_name = MP_3CA,
      MP_self_name = MP_slef,
      jaccard_similarity_score = similarity_score
    )
    jaccard_result <- rbind(jaccard_result,current_result)
    
  }
}

# jaccard plot ----------------------------------------------------------------
X_order <- c("IMP23","IMP3","IMP34","IMP14","IMP10","IMP25","IMP30","IMP38","IMP5","IMP26","IMP17","IMP6","IMP19",
             "IMP33","IMP12","IMP8","IMP16","IMP39","IMP2","IMP27","IMP9","IMP31","IMP13","IMP18","IMP24","IMP4","IMP7",
             "IMP1","IMP11","IMP15","IMP20","IMP21","IMP22","IMP28","IMP29","IMP32","IMP35","IMP36","IMP37","IMP40"
             # "IMP18",
             # "IMP1","IMP11","IMP13","IMP19","IMP25","IMP30","IMP36","IMP41","IMP43","IMP46","IMP47","IMP50"
)
jaccard_result$MP_self_name <- factor(jaccard_result$MP_self_name, levels = X_order)

## jat text
ggplot(jaccard_result, aes(x = MP_self_name, y = MP_3CA_name, fill = jaccard_similarity_score)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", jaccard_similarity_score)), vjust = 0.5, size = 2, fontface = "bold") +
  scale_fill_gradient2(low = "#6D9EC1",mid = 'white', high = "#E46726") +
  theme_minimal()+
  scale_y_discrete(limits = colnames(MP_3CA_hub)) +
  theme(plot.title = element_text(face = "bold",size = 12,color="black",hjust = 0.5),
        axis.title = element_text(face = "bold",size = 12,color ="black"), 
        axis.text = element_text(face = "bold",size= 16,color = "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1 ),
        panel.grid=element_blank(),
        legend.text = element_text(face = "bold",size= 12),
        legend.title= element_text(face = "bold",size= 12)
  ) + labs(fill = "pearson") + xlab("") + ylab("")
function_plot(filename_prefix = str_glue("{output_dir}/MP_cor_text"), width = 480, height = 500)

## jat text
ggplot(jaccard_result, aes(x = MP_self_name, y = MP_3CA_name, fill = jaccard_similarity_score)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#6D9EC1",mid = 'white', high = "#E46726") +
  theme_minimal()+
  scale_y_discrete(limits = colnames(MP_3CA_hub)) +
  theme(plot.title = element_text(face = "bold",size = 12,color="black",hjust = 0.5),
        axis.title = element_text(face = "bold",size = 12,color ="black"), 
        axis.text = element_text(face = "bold",size= 16,color = "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1 ),
        panel.grid=element_blank(),
        legend.text = element_text(face = "bold",size= 12),
        legend.title= element_text(face = "bold",size= 12)
  ) + labs(fill = "pearson") + xlab("") + ylab("")
function_plot(filename_prefix = str_glue("{output_dir}/MP_cor_notext"), width = 520, height = 500)

# data save ----------------------------------------------------------------
write_csv(jaccard_result,str_glue("{output_dir}/3CA_jaccard_result.csv"))
