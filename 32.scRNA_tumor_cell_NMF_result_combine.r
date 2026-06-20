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
library("ComplexHeatmap")
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

# data dir ----------------------------------------------------------------
all_gene_dir <- " "
input_dir <- " "
output_dir <- " "
set.seed(123)

# data combine ----------------------------------------------------------------
## all gene choose
scRNA_tumor_origin <- read_h5(file = str_glue("{all_gene_dir}/scRNA_gene_filter.h5"),
                         assay.name = 'RNA', 
                         target.object = 'seurat')
all_gene <- rownames(scRNA_tumor_origin)

remove(scRNA_tumor_origin)
gc()

## sample read
sample_list <- readRDS(str_glue("{input_dir}/sample_list.RData"))
# sample_list <- list.files(
#   path = str_glue("{input_dir}/02.tumor_cell_NMF_process_result/"),
#   pattern = "_NMF\\.RData$"
# )

## data read and combine
all_programs_top50 <- list()
for(sample_id in sample_list){
    # sample_id <- gsub("_NMF\\.RData$", "", sample_id)
    current_nmf <- readRDS(str_glue("{input_dir}/02.tumor_cell_NMF_process_result/{sample_id}_NMF.RData"))
    all_programs_top50[[sample_id]] <- as.matrix(as.data.frame(current_nmf))
}



# robust NMF analysis -------------------------------------------------------
nmf_filter_ccle <- robust_nmf_programs(all_programs_top50, intra_min = 35, intra_max = 10, inter_filter = TRUE, inter_min = 10)
all_programs_top50 <- lapply(all_programs_top50, function(x) x[, is.element(colnames(x), nmf_filter_ccle), drop = FALSE])
all_programs_top50 <- do.call(cbind, all_programs_top50)

# jar cor first -------------------------------------------------------
nmf_intersect <- apply(all_programs_top50, 2, function(x) apply(all_programs_top50, 2, function(y) length(intersect(x, y))))

nmf_intersect_hc <- hclust(as.dist(50 - nmf_intersect), method = "ward.D2")
nmf_intersect_hc <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))
nmf_intersect <- nmf_intersect[order.dendrogram(nmf_intersect_hc),order.dendrogram(nmf_intersect_hc)]
pdf(str_glue("{output_dir}/NMF_count_heatmap.pdf"),width = 10,height = 9)
pheatmap(nmf_intersect,
         cluster_rows = T, cluster_cols = T,
         clustering_method = 'ward.D2',
         show_colnames = F, show_rownames = F,
        #  breaks = seq(0, 30, length.out = 333),
         color = custom_magma,
         use_raster = F
         )
dev.off()


# jar cor analysis -------------------------------------------------------
## create matrix
program_names <- colnames(all_programs_top50)
n_prog <- length(program_names)
jaccard_mat <- matrix(0, nrow = n_prog, ncol = n_prog,dimnames = list(program_names, program_names))

## cor process
for(i in seq_len(n_prog)){  
    for(j in seq_len(n_prog)){    
        set1 <- all_programs_top50[,i]
        set2 <- all_programs_top50[,j]
        inter_len <- length(intersect(set1, set2))    
        union_len <- length(union(set1, set2))    
        jaccard_mat[i, j] <- inter_len / union_len  
    }
}

## jar cor cluster
dist_mat <- as.dist(1 - jaccard_mat)
hc <- hclust(dist_mat, method = "ward.D2")
# hc <- reorder(as.dendrogram(hc), colMeans(jaccard_mat))
# jaccard_mat <- nmf_intersect[order.dendrogram(hc),order.dendrogram(hc)]

## MP num confirm
k_meta <- 50
cluster_cut <- cutree(hc, k = k_meta)
program_cluster_df <- data.frame(  
    Program = program_names,  
    Cluster = cluster_cut)
table(program_cluster_df$Cluster)

## name MP
meta_program_names <- c("IMP1",                      
                        "IMP2",                        
                        "IMP3",                        
                        "IMP4",                        
                        "IMP5",
                        "IMP6",                        
                        "IMP7",                        
                        "IMP8",                        
                        "IMP9",                        
                        "IMP10",
                        "IMP11",                        
                        "IMP12",                        
                        "IMP13",                        
                        "IMP14",                        
                        "IMP15",
                        "IMP16",                        
                        "IMP17",
                        "IMP18",                        
                        "IMP19",                        
                        "IMP20",                        
                        "IMP21",
                        "IMP22",
                        "IMP23",
                        "IMP24",                   
                        "IMP25",                        
                        "IMP26",                        
                        "IMP27",
                        "IMP28",
                        "IMP29",
                        "IMP30",
                        "IMP31",
                        "IMP32",
                        "IMP33",
                        "IMP34",                        
                        "IMP35",                        
                        "IMP36",                        
                        "IMP37",
                        "IMP38",
                        "IMP39",
                        "IMP40",
                        "IMP41",
                        "IMP42",
                        "IMP43",
                        "IMP44",                        
                        "IMP45",                        
                        "IMP46",                        
                        "IMP47",
                        "IMP48",
                        "IMP49",
                        "IMP50"
                        )
program_cluster_df$MetaProgramName <- meta_program_names[program_cluster_df$Cluster]
head(program_cluster_df)

# MP character -------------------------------------------------------
meta_program_top50 <- list()

for(meta in unique(program_cluster_df$MetaProgramName)){    
    # meta <- "MP1"
    progs <- program_cluster_df$Program[program_cluster_df$MetaProgramName == meta]
    genes <- unlist(all_programs_top50[,progs])
    gene_freq <- as.data.frame(table(genes))
    colnames(gene_freq) <- c("gene_name","num")
    gene_freq <- arrange(gene_freq,desc(num))
    cutoff_level <- 0.2 * length(progs)
    gene_freg_filter <- gene_freq %>% 
        filter(num >= cutoff_level) %>% 
        head(50)
    # gene_freq <- table(genes)
    # gene_freq_sorted <- sort(gene_freq, decreasing = TRUE)
    # top50 <- names(gene_freq_sorted)[1:50]
    top50 <- as.character(gene_freg_filter$gene_name)
    if(length(top50) >= 10){
        meta_program_top50[[meta]] <- top50
    }
    
}

program_cluster_df_anno=data.frame(row.names = rownames(program_cluster_df),MetaProgramName=program_cluster_df$MetaProgramName)
labels_to_show <- program_cluster_df_anno$MetaProgramName
row_ha <- rowAnnotation(
  imp_name = anno_text(
    labels_to_show, 
    gp = gpar(fontsize = 0.5),
    just = "left"
  )
)

pdf(str_glue("{output_dir}/NMF_final_heatmap.pdf"),width = 19,height = 19)
p <- ComplexHeatmap::pheatmap(
    jaccard_mat,
    annotation_row = program_cluster_df_anno,
    show_colnames = FALSE,
    show_rownames = FALSE,
    treeheight_row = 20,
    treeheight_col = 20,
    cluster_rows = hc,
    cluster_cols = hc,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    color = custom_magma,
    use_raster = FALSE,
    breaks = seq(0.01, 0.25, length.out = 333)
    
)
p + row_ha
dev.off()

# enrichment analysis -------------------------------------------------------
# meta_program_top50 <- readRDS(" ")
enrich_GO <- lapply(meta_program_top50, function(program){runGSEA(program, universe=all_gene,
                    category = "C5", subcategory = "GO:BP")})
enrich_Hallmarker <- lapply(meta_program_top50, function(program){runGSEA(program, universe=all_gene,
                            category = "H")})
enrich_GO[["IMP1"]]
enrich_Hallmarker[["IMP1"]]

# data save -------------------------------------------------------
if(T){
    # write_csv(as.data.frame(meta_program_top50),str_glue("{output_dir}/meta_program_top50.csv"))
    write_csv(program_cluster_df,str_glue("{output_dir}/program_cluster_df.csv"))
    saveRDS(meta_program_top50,str_glue("{output_dir}/meta_program_top50.RData"))
    saveRDS(all_programs_top50,str_glue("{output_dir}/all_programs_top50.RData"))
    saveRDS(jaccard_mat,str_glue("{output_dir}/jaccard_matrix.RData"))
}









