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
library("glmGamPoi") 
library("paletteer")
library("NMF")
library("BiocParallel")
library("foreach")
library("doParallel")
(.packages())
# setwd(' ')
display.brewer.all()
nmf.options(maxIter = 1000, verbose = 1)
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

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

# data dir ----------------------------------------------------------------
origin_dir <- " "
cluster_dir <- " "
output_dir <- " "
nor_method <- "SCT"

# data read ----------------------------------------------------------------
## exp data
scRNA_tumor_transform <- read_h5(file = str_glue("{cluster_dir}/scRNA_cell_type_cluster_only.h5"),
                         assay.name = 'RNA', 
                         target.object = 'seurat')
scRNA_tumor_origin <- read_h5(file = str_glue("{origin_dir}/scRNA_gene_filter.h5"),
                         assay.name = 'RNA', 
                         target.object = 'seurat')                        


## copy data
scRNA_tumor_origin@reductions <- scRNA_tumor_transform@reductions

## space free
remove(scRNA_tumor_transform)
gc()

## hig var gene read
# high_gene <- read_csv(str_glue("{cluster_dir}/high_variable_gene_2w.csv"))
# high_gene <- high_gene[!is.na(high_gene$gene_name),]
high_gene <- rownames(scRNA_tumor_origin)[!str_detect(rownames(scRNA_tumor_origin),"^ENSG")]

# create and save loom data ------------------------------------------------

## test filter
scRNA_tumor_origin <- subset(scRNA_tumor_origin,tumor_code %in% c("ESCC"))

## hub filter
counts <- as.matrix(scRNA_tumor_origin@assays$RNA@counts)
# counts <- t(counts)
counts <- t(counts[high_gene,360000:552953])
# colnames(counts) <- colnames(scRNA_tumor_origin)
# rownames(counts) <- rownames(scRNA_tumor_origin)
write.csv(counts,str_glue("{output_dir}/scRNA_tumor_cell_counts_part3.csv"))
write.csv(counts,str_glue("{output_dir}/scRNA_tumor_cell_counts_test.csv"))

remove(counts)
gc()

# sample ID split ----------------------------------------------------------------
scRNA_tumor.list <- SplitObject(scRNA_tumor_origin, split.by = "sample_ID")
sample_list <- names(scRNA_tumor.list)

# exp filter ----------------------------------------------------------------
for(sample_id in sample_list){
    # sample_id <- "HRS1028865"
    print(str_glue("   -------------{sample_id}"))
    current_seurat <- scRNA_tumor.list[[sample_id]]

    ## scale and nor
    if(nor_method == "NormalizeData"){
        current_seurat <- NormalizeData(current_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
        current_seurat <- FindVariableFeatures(current_seurat, selection.method = "vst", nfeatures = 5000)
        current_seurat <- ScaleData(current_seurat, features = rownames(current_seurat), do.scale = FALSE, do.center = TRUE)
        exp_matrix <- current_seurat@assays$RNA@scale.data
    }else if(nor_method == "SCT"){
        current_seurat <- SCTransform(current_seurat,vst.flavor = "v2", method = "glmGamPoi", do.scale = FALSE, do.center = TRUE, return.only.var.genes = FALSE)
        exp_matrix <- current_seurat@assays$SCT@scale.data
    }else if(nor_method == "CPM_origin"){
        raw_counts <- GetAssayData(current_seurat, assay = "RNA", slot = "counts")
        raw_counts <- as.matrix(raw_counts)
        log2_cpm <- log2(t(t(raw_counts)/colSums(raw_counts))*100000+1)
        filter_matrix <- log2_cpm[apply(log2_cpm, 1, function(x) length(which(x > 3.5)) > ncol(log2_cpm)*0.02),]
        exp_matrix <- filter_matrix - rowMeans(filter_matrix)
    }

    ## select matrix and transform
    exp_matrix[exp_matrix < 0] <- 0
    exp_matrix <- exp_matrix[rowSums(exp_matrix) > 0,]

    ## data save
    saveRDS(exp_matrix,str_glue("{output_dir}/01.tumor_cell_exp_filter_SCT/{sample_id}_exp.RData"))

}

## sample list save
saveRDS(sample_list,str_glue("{output_dir}/sample_list.RData"))








