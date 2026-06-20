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
library("SCP")
library("Startrac")
library("corrplot")
library("Hmisc")
(.packages())
# setwd(' ')
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
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
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5 )) +
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
all_tumor_dir <- " "
MP_dir <- " "
output_dir <- " "

# data color ----------------------------------------------------------------
final_MP <- c("Cell.Cycle.G2M"="#e0bc58","Cell.Cycle.G1S"="#64abc0","Cell.Cycle.HMG"="#fab37f","Chromatin"="#e98741",
              "Stress"="#8fc0dc","Hypoxia"="#967568","Stress.in.vitro"="#f2d3ca","Protein.maturation"="#eebd85",
              "EMT.I"="#82c785","MHC"="#edeaa4","Epithelial.Senescence"="#cdaa9f","MYC"="#794976",
              "Respriration"="#bcacd3","Secreted.I"="#889b5d","Cilia.I"="#4e9592","PDAC.classcial"="#dbad5f",
              "Alveolar"="#64ae79","PDAC.related"="#ac5092","Lymphocyte.Activation"="#dc8e97","Cell.Signaling"="#e3d1db",
              "Coagulation"="#74a893","Neuron"="#ac9141",
              "Cell.Cycle.single.nucleus"="#5ac6e9","EMT.II"="#ebce8e",
              "EMT.VI"="#e5c06e","Secreted.II"="#7587b1","Cilia.II"="#c7deef","Colon.related"="#e97371","Other"="#1CB232"
              # ""="#e1a4c6",""="#916ba6",
              # ""="#cb8f82",""="#7db3af",""="#d2e0ac"
)
tumor_type_color <- c("ACC"="#5A7B8F","ANS"="#EE4C97","BRCA"="#3D806F","CESC"="#A08634","CRC"="#F37C95","ESCC"="#608541",
                      "HNSC"="#7D4E57","KIRC"="#BC3C29","LUAD"="#958056","LUSC"="#9FAFA3","NPC"="#6F99AD","PAAD"="#0072B5",
                      "PNET"="#CFC59A","SARC"="#E18727","SKCM"="#FFDC91","GC"="#718DAE","TGCT"="#7876B1","THCA"="#F9AC93",
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



# data read ---------------------------------------------------------------
## read
MP_sample_stat <- read_csv(str_glue("{all_tumor_dir}/MP_per_sample.csv"))
MP_meta_data <- read_csv(str_glue("{all_tumor_dir}/cell_ano_all.csv"))
meta_program_top50 <- readRDS(str_glue("{all_tumor_dir}/meta_program_top50.RData"))
scRNA_tumor_origin <- readRDS(str_glue("{all_tumor_dir}/scRNA_ano_MP.RData"))
scRNA_tumor_meta <- scRNA_tumor_origin@meta.data %>% 
  rownames_to_column(var = "cell_name")

## matrix filter
# exp_matrix <- as.matrix(scRNA_tumor_origin@assays$RNA@data)


# data plot ---------------------------------------------------------------
for(MP_name in names(meta_program_top50)){
  # MP_name <- "Cell.Cycle.G2M"
  
  ## hub gene filter
  MP_gene <- head(meta_program_top50[[MP_name]],30)
  
  ## cell sort
  scRNA_tumor_meta$sort_MP <- scRNA_tumor_meta[[MP_name]]
  scRNA_tumor_meta <- arrange(scRNA_tumor_meta,desc(sort_MP))
  cell_name_sort <- scRNA_tumor_meta$cell_name
  
  ## matrix filter
  MP_matrix <- FetchData(scRNA_tumor_origin, vars = MP_gene)
  MP_matrix <- t(MP_matrix)
  MP_matrix <- MP_matrix[,cell_name_sort]
  
  ## data plot
  pdf(str_glue("{output_dir}/NMF_marker_exp_{MP_name}.pdf"),width = 6,height = 7)
  # png(str_glue("{output_dir}/NMF_marker_exp_{MP_name}.png"),width = 6,height = 7)
  
  pheatmap(
    MP_matrix,
    show_colnames = FALSE,
    show_rownames = TRUE,
    cluster_rows = F,
    cluster_cols = F,
    cellwidth = 0.0005,cellheight = 12,
    breaks = seq(0, 4, length.out = 333),
    color = custom_magma,
    use_raster = F
    
  )
  dev.off()
  
}

# MP plot -----------------------------------------------------------------

## data filter
color_data <- MP_meta_data[,c("cellid","study_ID","tumor_code","tissue_type","tumor_status")]
  # column_to_rownames(var = "cellid")
MP_data <- MP_meta_data[,c(1,25:52)] %>% 
  column_to_rownames(var = "cellid")

## data transform
## data normal
scale_MP_matrix <- t(scale(MP_data))
scale_MP_matrix[scale_MP_matrix< -2] <- -2
scale_MP_matrix[scale_MP_matrix>2] <- 2
normalization <- function(x){return((x-min(x))/(max(x)-min(x)))} 
nor_MP_matrix <- normalization(scale_MP_matrix)


## new sorted
if(T){
  dominant_clusters <- apply(nor_MP_matrix, 2, function(x) rownames(nor_MP_matrix)[which.max(x)])
  ordered_samples <- character(0)
  for (cluster in rownames(nor_MP_matrix)) {
    candidate_samples <- names(dominant_clusters)[dominant_clusters == cluster]
    candidate_samples <- setdiff(candidate_samples, ordered_samples)
    samples_sorted <- names(sort(nor_MP_matrix[cluster, candidate_samples], decreasing = TRUE))
    
    ordered_samples <- c(ordered_samples, samples_sorted)
  }
  sorted_mat <- nor_MP_matrix[, ordered_samples]
  
}

## result plot
ano_row <-  as.data.frame(rownames(sorted_mat))
rownames(ano_row) <- rownames(sorted_mat)
colnames(ano_row) <- "final_MP"
ano_col <- as.data.frame(colnames(sorted_mat))
colnames(ano_col) <- "cellid"
ano_col <- left_join(ano_col,color_data,by = "cellid")
ano_col <- column_to_rownames(ano_col,var = "cellid")
pdf(str_glue("{output_dir}/scRNA_AUcell_pheatmap.pdf"), width = 19, height = 6)
pheatmap(sorted_mat,
         cluster_row = F, cluster_cols = F, show_rownames = T, 
         annotation_row = ano_row,
         annotation_col = ano_col,
         show_colnames = F, legend = T, annotation_names_row = F, annotation_names_col = F,
         annotation_colors = list(final_MP = c(final_MP),tumor_code = c(tumor_type_color),
                                  tissue_type = c(tissue_type_color),tumor_status = c(tumor_status_color),
                                  study_ID = c(study_database_color)),
         cellwidth = 0.002,cellheight = 14,
         fontsize = 5,
         use_raster = F
)
dev.off()


# MP all marker pheatmap --------------------------------------------------
cell_order <- colnames(sorted_mat)
unique_list <- lapply(seq_along(meta_program_top50),function(i){setdiff(meta_program_top50[[i]],unlist(meta_program_top50[-i]))})


for(MP_name in names(meta_program_top50)){
  # MP_name <- "Stress.in.vitro"
  
  ## hub gene filter
  MP_gene <- head(meta_program_top50[[MP_name]],50)
  
  ## matrix filter
  MP_matrix <- FetchData(scRNA_tumor_origin, vars = MP_gene)
  MP_matrix <- t(MP_matrix)
  MP_matrix <- MP_matrix[,cell_order]
  
  ## data plot
  pdf(str_glue("{output_dir}/all_marker_plot/NMF_all_marker_{MP_name}.pdf"),width = 12,height = 4)
  pheatmap(
    MP_matrix,
    show_colnames = FALSE,
    show_rownames = FALSE,
    cluster_rows = F,
    cluster_cols = F,
    cellwidth = 0.001,cellheight = 6,
    breaks = seq(0, 4, length.out = 333),
    color = custom_magma,
    use_raster = F
    
  )
  dev.off()
  
  remove(MP_matrix)
  gc()
  
}


# NMF other plot ----------------------------------------------------------
## NMF sample stat
MP_sample_stat <- MP_sample_stat %>% 
  mutate(status = ifelse(percent_tumor > 0.5 | percent_tissue > 0.5, "high.in.50%",
                         ifelse(percent_tumor > 0.3 | percent_tissue > 0.3, "high.in.30%","low.percentage")))

ggplot(MP_sample_stat,aes(percent_tumor, percent_tissue)) +
  geom_point(aes(size = 4, color = status)) + 
  geom_vline(xintercept = c(0.3,0.5), lty = 4, col = '#a4b0be', lwd = 1) +
  geom_hline(yintercept = c(0.3,0.5), lty = 4, col = '#a4b0be', lwd = 1) + 
  geom_text_repel(aes(label = MP_name),box.padding = 0.5, min.segment.length = 0,color = 'black',
                  max.overlaps = getOption("ggrepel.max.overlaps",default = 100)) +
  xlim(c(0,1)) + ylim(c(0,1)) +
  nrc_theme + NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/MP_percentage_stat"), width = 200, height = 200)
