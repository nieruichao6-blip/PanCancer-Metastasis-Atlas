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
bulkRNA_dir <- " "
scRNA_dir <- " "
output_dir <- " "


# data color --------------------------------------------------------------
ano_colors <- list(
  final_MP = c("Cell.Cycle.G2M"="#e0bc58","Cell.Cycle.G1S"="#64abc0","Cell.Cycle.HMG"="#fab37f","Chromatin"="#e98741",
                "Stress"="#8fc0dc","Hypoxia"="#967568","Stress.in.vitro"="#f2d3ca","Protein.maturation"="#eebd85",
                "EMT.I"="#82c785","MHC"="#edeaa4","Epithelial.Senescence"="#cdaa9f","MYC"="#794976",
                "Respriration"="#bcacd3","Secreted.I"="#889b5d","Cilia.I"="#4e9592","PDAC.classcial"="#dbad5f",
                "Alveolar"="#64ae79","PDAC.related"="#ac5092","Lymphocyte.Activation"="#dc8e97","Cell.Signaling"="#e3d1db",
                "Coagulation"="#74a893","Neuron"="#ac9141",
                "Cell.Cycle.single.nucleus"="#5ac6e9","EMT.II"="#ebce8e",
                "EMT.VI"="#e5c06e","Secreted.II"="#7587b1","Cilia.II"="#c7deef","Colon.related"="#e97371"
                # ""="#e1a4c6",""="#916ba6",
                # ""="#cb8f82",""="#7db3af",""="#d2e0ac"
  ),
  tumor_code = c("ACC"="#5A7B8F","PRAD"="#EE4C97","BRCA"="#3D806F","CESC"="#A08634","CRC"="#F37C95","ESCA"="#608541",
                        "HNSC"="#7D4E57","KIRC"="#BC3C29","LUAD"="#958056","LUSC"="#9FAFA3","PCPG"="#6F99AD","PAAD"="#0072B5",
                        "READ"="#CFC59A","SARC"="#E18727","SKCM"="#FFDC91","GC"="#718DAE","TGCT"="#7876B1","THCA"="#F9AC93",
                        "NSCLC"="#E3BC06","THYM"="#2285FE","UVM"="#3E6086","UCEC"="#20854E","UCS"="#4A7985","GBM"="#9F4700",
                        "OV"="#55d151","MESO"="#f69a95","LIHC"="#aa329a","KIRP"="#EDE816","KICH"="#48dc88","DLBC"="#48dc88",
                        "CHOL"="#9ebe71","LGG"="#FDC44D","LAML"="#16A795","STAD"="#6FCDDC","COAD"="#E7ABFD","BLCA"="#CEE769")
  
)
tumor_code = c("ACC"="#5A7B8F","PRAD"="#EE4C97","BRCA"="#3D806F","CESC"="#A08634","CRC"="#F37C95","ESCA"="#608541",
               "HNSC"="#7D4E57","KIRC"="#BC3C29","LUAD"="#958056","LUSC"="#9FAFA3","PCPG"="#6F99AD","PAAD"="#0072B5",
               "READ"="#CFC59A","SARC"="#E18727","SKCM"="#FFDC91","GC"="#718DAE","TGCT"="#7876B1","THCA"="#F9AC93",
               "NSCLC"="#E3BC06","THYM"="#2285FE","UVM"="#3E6086","UCEC"="#20854E","UCS"="#4A7985","GBM"="#9F4700",
               "OV"="#55d151","MESO"="#f69a95","LIHC"="#aa329a","KIRP"="#EDE816","KICH"="#48dc88","DLBC"="#48dc88",
               "CHOL"="#9ebe71","LGG"="#FDC44D","LAML"="#16A795","STAD"="#6FCDDC","COAD"="#E7ABFD","BLCA"="#CEE769")


# data read ---------------------------------------------------------------
## DEG read
module_gene_list <- readRDS(str_glue("{scRNA_dir}/meta_program_top50.RData"))

## bulkRNA fpkm read
bulkRNA_fpkm <- read_csv(str_glue("{bulkRNA_dir}/bulkRNA_fpkm_filter.csv"))

## bulkRNA sample read
bulkRNA_sample <- read_csv(str_glue("{bulkRNA_dir}/bulkRNA_smaple_filter.csv"))


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
colnames(ano_row) <- "final_MP"
ano_col <- as.data.frame(colnames(sorted_mat))
colnames(ano_col) <- "sample"
ano_col <- left_join(ano_col,bulkRNA_sample_filter[,c(1,8)],by = "sample")
rownames(ano_col) <- bulkRNA_sample_filter$sample
ano_col$tumor_code <- factor(ano_col$tumor_code)
ano_col <- ano_col[-1]

pdf(str_glue("{output_dir}/bulkRNA_GSVA.pdf"), width = 19, height = 6)
pheatmap(sorted_mat,
         cluster_row = F, cluster_cols = F, show_rownames = T, 
         annotation_row = ano_row,
         annotation_col = ano_col,
         show_colnames = F, legend = T, annotation_names_row = F, annotation_names_col = F,
         annotation_colors = ano_colors,
         cellwidth = 0.1,cellheight = 14,
         fontsize = 5,
         use_raster = F
)
dev.off()



# hub MP per --------------------------------------------------------------
bulkRNA_gsva_sorted <- as.data.frame(sorted_mat)
bulkRNA_gsva_sorted <- rownames_to_column(bulkRNA_gsva_sorted,var = "final_MP")
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
