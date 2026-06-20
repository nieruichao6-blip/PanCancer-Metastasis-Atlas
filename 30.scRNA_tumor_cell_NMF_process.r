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
library("tools")
library("tidyverse")
library("optparse")
library("reshape2")
library("ggfortify")
library("ggpubr")
library("scales")
library('glue')
library('funr')
library('xlsx')
library("ggrepel")
library("ggridges")
library("statmod")
library("S4Vectors")
library("stats4")
library("SummarizedExperiment")
library("devtools")
library("Seurat")
library("randomcoloR")
library("patchwork")
library("ggprism")
library("plyr")
library("randomcoloR")
library("future")
library("ggforce")
library("ggsci")
library("doFuture")
library("msigdbr")
library("fgsea")
library("NMF")
library("BiocParallel")
library("foreach")
library("doParallel")

nmf.options(maxIter = 12000, verbose = 1)

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

# set opt and args --------------------------------------------------------
option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = " ", action = "store", help = "exp input dir"),
    make_option(c("-f", "--file_name"), type = "character", default = "GSM", action = "store", help = "exp file name"),
    make_option(c("-n", "--NMF_method"), type = "character", default = "origin_nmf", action = "store", help = "NMF method choose in origin_nmf brunet_nmf snmfr_nmf Rcpp_nmf"),
    make_option(c("-k", "--rank_k"), type = "character", default = "4 5 6 7 8 9", action = "store", help = "rank choose"),
    make_option(c("-g", "--high_gene_num"), type = "double", default = 50, action = "store", help = "rank high gene num"),
    make_option(c("-I", "--Iteration_num"), type = "double", default = 2000, action = "store", help = "Iteration_num"),
    make_option(c("--outputdir"), type = "character", default = " ", action = "store", help = "output dir")
)

opt = parse_args(OptionParser(option_list = option_list, usage = "pipeline for NMF process"))
input = opt$input
output_dir = opt$outputdir
sample_id = opt$file_name
NMF_process = opt$NMF_method
rank = opt$rank_k
hub_gene_filter = opt$high_gene_num
Iteration = opt$Iteration_num


# data read --------------------------------------------------------
scRNA_sample_data <- readRDS(str_glue("{input}/{sample_id}_exp.RData"))

# NMF process --------------------------------------------------------
top50_genes_list_all <- list()
top50_genes_list_k <- list()
rank <- as.numeric(strsplit(rank, " ")[[1]])
for(k in rank){
    print(str_glue("  -----------rank - {k}"))
    if(NMF_process == "origin_nmf"){
        res_nmf <- nmf(scRNA_sample_data, k, seed = 123, method = "snmf/r",.options = "p100v1")
    }else if(NMF_process == brunet_nmf){
        res_nmf <- nmf(scRNA_sample_data, k, seed = "nndsvd", method = "brunet",maxIter = Iteration,.options = "p100v1")
    }else if(NMF_process == snmfr_nmf){
        res_nmf <- nmf(scRNA_sample_data, k, seed = "nndsvd", method = "snmf/r",maxIter = Iteration,.options = "p100v1")
    }else if(NMF_process == "Rcpp_nmf"){
        res_nmf <- RcppML::nmf(scRNA_sample_data, k, L1 = c(0,0), verbose=FALSE, seed=123)
    }

    W <- basis(res_nmf) 
    top50_genes_list_k <- apply(W, 2, function(x) {      
        top_idx <- order(x, decreasing = TRUE)[1:hub_gene_filter]    
        rownames(W)[top_idx]
    })
    colnames(top50_genes_list_k) <- paste0(sample_id, "_Program", k, 1:k) 

    for(program_name in colnames(top50_genes_list_k)){    
        top50_genes_list_all[[ program_name ]] <- top50_genes_list_k[, program_name]
    }
}

  
# data save --------------------------------------------------------
saveRDS(top50_genes_list_all,str_glue("{output_dir}/{sample_id}_NMF.RData"))
saveRDS(res_nmf,str_glue("{output_dir}/{sample_id}_NMF_origin.RData"))
