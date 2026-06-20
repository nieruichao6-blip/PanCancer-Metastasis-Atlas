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
setwd(' ')
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

source(" ")
source(" ")

# data dir ----------------------------------------------------------------
origin_dir <- " "
cluster_dir <- " "
output_dir <- " "

# data read ----------------------------------------------------------------
## high var gene
high_var_gene <- read_csv(str_glue("{cluster_dir}/high_variable_gene.csv"))

## exp data
scRNA_tumor_transform <- read_h5(file = str_glue("{cluster_dir}/scRNA_cell_type_cluster_only.h5"),
                         assay.name = 'RNA', 
                         target.object = 'seurat')
scRNA_tumor_origin <- read_h5(file = str_glue("{origin_dir}/scRNA_gene_filter.h5"),
                         assay.name = 'RNA', 
                         target.object = 'seurat')                        
# scRNA_tumor@assays$RNA@var.features <- high_var_gene$gene_name
# scRNA_tumor <- FindVariableFeatures(scRNA_tumor, nfeatures = 3000)

## copy data
scRNA_tumor_origin@reductions <- scRNA_tumor_transform@reductions

## space free
remove(scRNA_tumor_transform)
gc()



# pre process ----------------------------------------------------------------
## nohub gene
mito.genes <- grep(pattern = "^MT-", x = rownames(scRNA_tumor_origin), value = TRUE)
rbl.genes <- grep(pattern = "^RB-", x = rownames(scRNA_tumor_origin), value = TRUE)
rsl.genes <- grep(pattern = "^RS-", x = rownames(scRNA_tumor_origin), value = TRUE)  
rpl.genes <- grep(pattern = "^RPL-", x = rownames(scRNA_tumor_origin), value = TRUE)  
rbl.genes <- grep(pattern = "^RBL-", x = rownames(scRNA_tumor_origin), value = TRUE)  
rps.genes <- grep(pattern = "^RPS-", x = rownames(scRNA_tumor_origin), value = TRUE)  
rbs.genes <- grep(pattern = "^RBS-", x = rownames(scRNA_tumor_origin), value = TRUE)  
rbl1.genes <- grep(pattern = "^RB", x = rownames(scRNA_tumor_origin), value = TRUE)  
rsl1.genes <- grep(pattern = "^RS", x = rownames(scRNA_tumor_origin), value = TRUE)  
rpl1.genes <- grep(pattern = "^RPL", x = rownames(scRNA_tumor_origin), value = TRUE)  
rbl1.genes <- grep(pattern = "^RBL", x = rownames(scRNA_tumor_origin), value = TRUE)  
rps1.genes <- grep(pattern = "^RPS", x = rownames(scRNA_tumor_origin), value = TRUE)  
rbs1.genes <- grep(pattern = "^RBS", x = rownames(scRNA_tumor_origin), value = TRUE)
hub_gene <- rownames(scRNA_tumor_origin)[!rownames(scRNA_tumor_origin) %in% c(mito.genes,rbl.genes,rsl.genes,rpl.genes,rbl.genes,rps.genes,
                                                  rbs.genes,rbl1.genes,rsl1.genes,rpl1.genes,rbl1.genes,rps1.genes,rbs1.genes)]

## filter hub data
scRNA_tumor_origin <- subset(scRNA_tumor_origin, features = hub_gene)

# sample ID split ----------------------------------------------------------------
scRNA_tumor.list <- SplitObject(scRNA_tumor_origin, split.by = "sample_ID")

all_gene <- rownames(scRNA_tumor_origin)

## space free
remove(scRNA_tumor_origin)
gc()

# run nmf ----------------------------------------------------------------
sample_list <- names(scRNA_tumor.list)[1:50]

register(MulticoreParam(workers = 100, progressbar = TRUE))
bpparam()

for(sample_id in sample_list){
    print(str_glue("   -------------{sample_id}"))
    temp_seurat <- scRNA_tumor.list[[sample_id]]
    gene_mat <- sample_nmf(temp_seurat,nor_method = "CPM_origin",nmf_method = "origin_nmf",rank = 4:9,gene_num = 50,
                           scale = TRUE,center = TRUE,sample = sample_id,output = output_dir,max_num = 1000)
    # all_programs_top50[sample_id] <- gene_mat

    saveRDS(gene_mat,str_glue("{output_dir}/{sample_id}_nmf.RData"))
    # for(program_name in colnames(gene_mat)){
    #     all_programs_top50[[ program_name ]] <- gene_mat[, program_name]  
    # }
}

register(SerialParam())
bpparam()
gc()

# nmf result read ----------------------------------------------------------------
all_programs_top50 <- list()

for(sample_id in sample_list){
    current_nmf <- readRDS(str_glue("{output_dir}/{sample_id}_nmf.RData"))
    all_programs_top50[[sample_id]] <- as.matrix(as.data.frame(current_nmf))

}

# robust NMF analysis -------------------------------------------------------
nmf_filter_ccle <- robust_nmf_programs(all_programs_top50, intra_min = 35, intra_max = 10, inter_filter = TRUE, inter_min = 10)
all_programs_top50 <- lapply(all_programs_top50, function(x) x[, is.element(colnames(x), nmf_filter_ccle), drop = FALSE])
all_programs_top50 <- do.call(cbind, all_programs_top50)

# jar cor first -------------------------------------------------------
nmf_intersect <- apply(all_programs_top50, 2, function(x) apply(all_programs_top50, 2, function(y) length(intersect(x, y))))

nmf_intersect_hc <- hclust(as.dist(50 - nmf_intersect), method = "average")
nmf_intersect_hc <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))
nmf_intersect <- nmf_intersect[order.dendrogram(nmf_intersect_hc),order.dendrogram(nmf_intersect_hc)]
pheatmap(nmf_intersect,
         cluster_rows = T, cluster_cols = T,
         clustering_method = 'average',
         show_colnames = F, show_rownames = F,
         # breaks = c(5, 25),
         color = custom_magma
         )

# jar cor analysis -------------------------------------------------------
## create matrix
program_names <- names(all_programs_top50)
n_prog <- length(program_names)
jaccard_mat <- matrix(0, nrow = n_prog, ncol = n_prog,dimnames = list(program_names, program_names))

## cor process
for(i in seq_len(n_prog)){  
    for(j in seq_len(n_prog)){    
        set1 <- all_programs_top50[[ i ]]    
        set2 <- all_programs_top50[[ j ]]    
        inter_len <- length(intersect(set1, set2))    
        union_len <- length(union(set1, set2))    
        jaccard_mat[i, j] <- inter_len / union_len  
    }
}

## jar cor cluster
dist_mat <- as.dist(1 - jaccard_mat)
hc <- hclust(dist_mat, method = "ward.D2")

## MP num confirm
k_meta <- 12
cluster_cut <- cutree(hc, k = k_meta)
program_cluster_df <- data.frame(  
    Program = program_names,  
    Cluster = cluster_cut)
table(program_cluster_df$Cluster)

## name MP
meta_program_names <- c("MP1",                        
                        "MP2",                        
                        "MP3",                        
                        "MP4",                        
                        "MP5",
                        "MP6",                        
                        "MP7",                        
                        "MP8",                        
                        "MP9",                        
                        "MP10",
                        "MP11",                        
                        "MP12",                        
                        "MP13",                        
                        "MP14",                        
                        "MP15")
program_cluster_df$MetaProgramName <- meta_program_names[program_cluster_df$Cluster]
head(program_cluster_df)

# MP character -------------------------------------------------------
meta_program_top50 <- list()

for(meta in unique(program_cluster_df$MetaProgramName)){    
    progs <- program_cluster_df$Program[program_cluster_df$MetaProgramName == meta]
    genes <- unlist(all_programs_top50[progs])
    gene_freq <- table(genes)
    gene_freq_sorted <- sort(gene_freq, decreasing = TRUE)
    top50 <- names(gene_freq_sorted)[1:50]
    meta_program_top50[[meta]] <- top50
}

# data save -------------------------------------------------------
if(T){
    write_csv(all_programs_top50)
}



