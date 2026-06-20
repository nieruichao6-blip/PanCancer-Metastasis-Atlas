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
library("ComplexHeatmap")
library("SCP")
library("SCENIC")
library("SCopeLoomR")
(.packages())
# setwd(' ')
display.brewer.all()
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

function_plot <- function(filename_prefix, width, height){
  ggsave(filename = str_c(filename_prefix, '.pdf'), plot = last_plot(), device = 'pdf', width = width, height = height, units = 'mm', dpi = 500)
  ggsave(filename = str_c(filename_prefix, '.svg'), plot = last_plot(), device = 'svg', width = width, height = height, units = 'mm', dpi = 500)
  ggsave(filename = str_c(filename_prefix, '.png'), plot = last_plot(), device = 'png', width = width, height = height, units = 'mm', dpi = 500)
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


# data dir ----------------------------------------------------------------
scRNA_dir <- " "
scenic_dir <- " "
NMF_dir <- " "
output_dir <- " "

all_cor <- data.frame(
  MP_name = as.character(),
  TF_name = as.character(),
  cor = as.numeric(),
  pvalue = as.numeric()
)

hub_MP <- c("Respriration","Stress.in.vitro","Cell.Cycle.HMG","Hypoxia",
            "Neuron","EMT.I","Lymphocyte.Activation","Stress","Cell.Cycle.G2M")

# data read ---------------------------------------------------------------
loom_part1 <- open_loom(str_glue("{scenic_dir}/02.pyscenic_process_origin_part1/aucell.loom"))
loom_part2 <- open_loom(str_glue("{scenic_dir}/02.pyscenic_process_origin_part2/aucell.loom"))
loom_part3 <- open_loom(str_glue("{scenic_dir}/02.pyscenic_process_origin_part3/aucell.loom"))
scRNA_tumor_cell <- read_h5(file = str_glue("{scRNA_dir}/scRNA_cell_type_cluster_only.h5"),
                            assay.name = 'RNA', 
                            target.object = 'seurat')
NMF_meta <- read_csv(str_glue("{NMF_dir}/addmodulescore_meta.csv"))
gc()

# data select -------------------------------------------------------------

## select
regulons_incidMat_part1 <- get_regulons(loom_part1,column.attr.name = "Regulons")
regulons_part1 <- regulonsToGeneLists(regulons_incidMat_part1)
regulonAUC_part1 <- get_regulons_AUC(loom_part1,column.attr.name = "RegulonsAUC")
regulonAucThresholds_part1 <- get_regulon_thresholds(loom_part1)
regulons_incidMat_part2 <- get_regulons(loom_part2,column.attr.name = "Regulons")
regulons_part2 <- regulonsToGeneLists(regulons_incidMat_part2)
regulonAUC_part2 <- get_regulons_AUC(loom_part2,column.attr.name = "RegulonsAUC")
regulonAucThresholds_part2 <- get_regulon_thresholds(loom_part2)
regulons_incidMat_part3 <- get_regulons(loom_part3,column.attr.name = "Regulons")
regulons_part3 <- regulonsToGeneLists(regulons_incidMat_part3)
regulonAUC_part3 <- get_regulons_AUC(loom_part3,column.attr.name = "RegulonsAUC")
regulonAucThresholds_part3 <- get_regulon_thresholds(loom_part3)

## all result combine
# sub_regulonAUC_part1 <- regulonAUC_part1[,match(colnames(scRNA_tumor_cell),colnames(regulonAUC_part1))]
sub_regulonAUC_part1 <- rownames_to_column(as.data.frame(regulonAUC_part1@assays@data$AUC),var = "motif")
# sub_regulonAUC_part2 <- regulonAUC_part2[,match(colnames(scRNA_tumor_cell),colnames(regulonAUC_part2))]
sub_regulonAUC_part2 <- rownames_to_column(as.data.frame(regulonAUC_part2@assays@data$AUC),var = "motif")
# sub_regulonAUC_part3 <- regulonAUC_part3[,match(colnames(scRNA_tumor_cell),colnames(regulonAUC_part3))]
sub_regulonAUC_part3 <- rownames_to_column(as.data.frame(regulonAUC_part3@assays@data$AUC),var = "motif")
all_regulonAUC <- full_join(sub_regulonAUC_part1,sub_regulonAUC_part3,by = "motif")
all_regulonAUC <- full_join(all_regulonAUC,sub_regulonAUC_part2,by = "motif")
all_regulonAUC <- column_to_rownames(all_regulonAUC,var = "motif")

## data transform
# scRNA_tumor_cell@meta.data <- cbind(scRNA_tumor_cell@meta.data,t(sub_regulonAUC@assays@data$AUC))
# scRNA_tumor_cell <- AddMetaData(scRNA_tumor_cell,t(sub_regulonAUC@assays@data$AUC))
all_motif <- rownames(all_regulonAUC)
scRNA_tumor_cell <- AddMetaData(scRNA_tumor_cell,t(all_regulonAUC))

# scRNA_tumor_cell@meta.data <- scRNA_tumor_cell@meta.data %>% 
#   dplyr::select(c(1:23))

# scenic data filter ------------------------------------------------------
## data filter
NMF_meta <- NMF_meta %>% 
  dplyr::select(c(18,24:51))
MP_name <- colnames(NMF_meta)[-1]
scenic_meta <- scRNA_tumor_cell@meta.data %>% 
  dplyr::select(c("cell_id",all_motif))
TF_name <- colnames(scenic_meta)[-1]

## data combine
combine_meta <- left_join(NMF_meta,scenic_meta,by = "cell_id")


# cor analysis ------------------------------------------------------------
for(MP in MP_name){
  
  for(TF in TF_name){
    
    cor_test <- cor.test(combine_meta[[MP]],combine_meta[[TF]],method = "pearson")
    current_cor <- data.frame(
      MP_name = MP,
      TF_name = TF,
      cor = cor_test$estimate,
      pvalue = cor_test$p.value
    )
    
    ## data rbind
    all_cor <- rbind(all_cor,current_cor)
    
    
  }
}



# cor result filter -------------------------------------------------------
all_cor_filter <- all_cor %>% 
  dplyr::group_by(MP_name) %>%
  dplyr::arrange(desc(cor),.by_group = T) %>% 
  slice_head(n = 5)




# phetmap plot ------------------------------------------------------------
## data
pheatmap_data <- all_cor %>% 
  dcast(MP_name ~ TF_name, value.var = "cor", fun.aggregate = function(x) x[1]) %>% 
  column_to_rownames(var = "MP_name")


if(T){
  dominant_clusters <- apply(pheatmap_data, 2, function(x) rownames(pheatmap_data)[which.max(x)])
  ordered_samples <- character(0)
  for (cluster in rownames(pheatmap_data)) {
    candidate_samples <- names(dominant_clusters)[dominant_clusters == cluster]
    candidate_samples <- setdiff(candidate_samples, ordered_samples)
    samples_sorted <- names(sort(pheatmap_data[cluster, candidate_samples], decreasing = TRUE))
    
    ordered_samples <- c(ordered_samples, samples_sorted)
  }
  sorted_mat <- pheatmap_data[, ordered_samples]
  
}

## all plot
all_cor_filter_3 <- all_cor %>% 
  dplyr::group_by(MP_name) %>%
  dplyr::arrange(desc(cor),.by_group = T) %>% 
  slice_head(n = 3)
ano_col <- data.frame(row.names = unique(all_cor$MP_name),
                                  final_MP=unique(all_cor$MP_name))

pos <- which(rownames(t(sorted_mat)) %in% unique(all_cor_filter_3$TF_name))
label <- rowAnnotation(Zscore = anno_mark(at = pos,labels = unique(all_cor_filter_3$TF_name),labels_gp = gpar(fontSize = 10)))
pdf(str_glue("{output_dir}/all_TF_result.pdf"),width = 12,height = 12)
ComplexHeatmap::pheatmap(
  t(sorted_mat),
  # annotation_row = program_cluster_df_anno,
  annotation_col = ano_col,
  show_colnames = TRUE,
  show_rownames = F,
  treeheight_row = 20,
  treeheight_col = 20,
  cluster_rows = F,
  cluster_cols = F,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  right_annotation = label,
  # color = custom_magma,
  col = colorRampPalette(c("darkblue","white","darkred"))(500),
  annotation_colors = list(final_MP = c(final_MP)),
  # breaks = seq(0.01, 0.25, length.out = 333)
)
dev.off()


# dot plot ----------------------------------------------------------------
all_cor <- read_csv(str_glue("{output_dir}/MP_TF_all_cor.csv"))
high_cor_data <- all_cor %>% 
  filter(MP_name %in% hub_MP) %>% 
  filter(cor >= 0.3)
high_cor_num <- as.data.frame(table(high_cor_data$TF_name))
colnames(high_cor_num) <- c("TF_name","num")
high_cor_num <- arrange(high_cor_num,desc(num))

## hub plot
hub_TF <- c("YBX1(+)","DRAP1(+)","SPDEF(+)",
            "DDIT3(+)","KLF4(+)","FOS(+)",
            "YBX1(+)","DRAP1(+)","SPDEF(+)",
            "MXI1(+)","EMX2(+)","MAF(+)",
            "GLIS3(+)","TEAD1(+)","FOXO3(+)",
            "SPI1(+)","ZNF689(+)","POU3F2(+)",
            "IKZF1(+)","RUNX3(+)","FLI1(+)",
            "ATF3(+)","FOSB(+)","JUN(+)",
            "MYBL2(+)","E2F1(+)","E2F8(+)")
hub_TF_data <- all_cor %>% 
  filter(MP_name %in% hub_MP) %>% 
  filter(TF_name %in% hub_TF)
ggplot(hub_TF_data, aes(x = TF_name, y = MP_name)) +
  geom_point(aes(color = cor),size = 8) +
  # geom_text(aes(label = sprintf("%.2f", cor)), vjust = 0.5, size = 2, fontface = "bold") +
  scale_color_gradient2(low = "#6D9EC1",mid = 'white', high = "#E46726") +
  # scale_size_continuous(range = c(-0.5,0.8),breaks = c(0,20,30,40,50,60,70,80)) +
  # theme_minimal()+
  scale_y_discrete(limits = hub_MP) +
  scale_x_discrete(limits = unique(hub_TF)) +
  nrc_theme + xlab("") + ylab("")
function_plot(filename_prefix = str_glue("{output_dir}/MP_scenic_cor"), width = 280, height = 120)


# data save ---------------------------------------------------------------
write_csv(all_cor,str_glue("{output_dir}/MP_TF_all_cor.csv"))
write_csv(all_cor_filter,str_glue("{output_dir}/MP_TF_all_cor_filter.csv"))




