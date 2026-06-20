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
library("ROGUE")
library("ggrastr")
library("slingshot")
library("CytoTRACE2")
library("ggalluvial")
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
                       axis.text.x = element_text(angle = 0, hjust = 0.5 ),
                       panel.grid=element_blank(),
                       legend.text = element_text(face = "bold",size= 12),
                       legend.title= element_text(face = "bold",size= 12))

# data color --------------------------------------------------------------
final_MP <- c("Cell.Cycle.G2M"="#e0bc58","Cell.Cycle.G1S"="#64abc0","Cell.Cycle.HMG"="#fab37f","Chromatin"="#e98741",
              "Stress"="#8fc0dc","Hypoxia"="#967568","Stress.in.vitro"="#f2d3ca","Protein.maturation"="#eebd85",
              "EMT.I"="#82c785","MHC"="#edeaa4","Epithelial.Senescence"="#cdaa9f","MYC"="#794976",
              "Respriration"="#bcacd3","Secreted.I"="#889b5d","Cilia.I"="#4e9592","PDAC.classcial"="#dbad5f",
              "Alveolar"="#64ae79","PDAC.related"="#ac5092","Lymphocyte.Activation"="#dc8e97","Cell.Signaling"="#e3d1db",
              "Coagulation"="#74a893","Neuron"="#ac9141","Cell.Cycle.single.nucleus"="#5ac6e9","EMT.II"="#ebce8e",
              "EMT.VI"="#e5c06e","Secreted.II"="#7587b1","Cilia.II"="#c7deef","Colon.related"="#e97371")
tumor_type_color <- c("ACC"="#5A7B8F","ANS"="#EE4C97","BRCA"="#3D806F","CESC"="#A08634","CRC"="#F37C95","ESCC"="#608541",
                      "HNSC"="#7D4E57","KIRC"="#BC3C29","LUAD"="#958056","LUSC"="#9FAFA3","NPC"="#6F99AD","PAAD"="#0072B5",
                      "PNET"="#CFC59A","SARC"="#E18727","SKCM"="#FFDC91","GC"="#718DAE","TGCT"="#7876B1","THCA"="#F9AC93",
                      "NSCLC"="#E3BC06","LSCC"="#2285FE","Liver_healthy"="#3E6086","Lymph_node_healthy"="#20854E","Brain_healthy"="#4A7985")
tissue_type_color <- c("Bone_marrow"="#CE1261","Kidney"="#5E4FA2","Lymph_node"="#8CA77B","Brain"="#00441B","Lung"="#ed3d30",
                       "Colon"="#DEDC00","Liver"="#B3DE69","Peritoneal"="#999999","Stomach"="#725e82","Ovary"="#8c4356",
                       "Head_and_neck"="#db5a6b","Pancreas"="#21a675","Vaginal"="#d9b611","Appendix"="#e29c45",
                       "Bone"="#4c8dae","Nasopharynx"="#003371","Salivary_gland"="#a1afc9","Testis"="#88ada6","Cervix"="#c93756",
                       "Skin"="#204B75","Breast"="#588257","Thyroid"="#B6DB7B","Esophagus"="#4285BF","larynx"="#1CB232")
tumor_status_color <- c("Tumor" = "#5BC0EB","Metastasis"="#8DD3C7","NAT" = "#9BC53D","Healthy" = "#C3423F")

# data dir ----------------------------------------------------------------
origin_dir <- " "
input_dir <- " "
output_dir <- " "

# data read ---------------------------------------------------------------
## data read
scRNA_tumor_cell <- read_h5(file = str_glue("{origin_dir}/scRNA_gene_filter.h5"),
                            assay.name = 'RNA',
                            target.object = 'seurat')
scRNA_scTour_meta <- read_csv(str_glue("{input_dir}/scRNA_tumor_cell_metadata_30.csv"))
scRNA_scTour_meta <- column_to_rownames(scRNA_scTour_meta,var = "index")

## meta copy
scRNA_tumor_cell@meta.data <- scRNA_scTour_meta

# data pre process --------------------------------------------------------
scRNA_tumor_cell <- NormalizeData(scRNA_tumor_cell,normalization.method = 'LogNormalize',scale.factor = 10000)


# hub gene exp ------------------------------------------------------------
expression_data.n <- GetAssayData(scRNA_tumor_cell,layer = "data")
expression_data.n.s <- as.data.frame(expression_data.n["LMNA",])
colnames(expression_data.n.s) <- "ARHGDIB"

data <- merge(expression_data.n.s, scRNA_tumor_cell@meta.data,by = 0)
data2 <- merge(data,scRNA_tumor_cell,by = 1,by.y = 0)

# data <- head(data,100000)
ggplot(data, aes(x = ptime, y = ARHGDIB,color = metastasis_tissue,fill = metastasis_tissue)) +
  geom_smooth(method = "gam", se = TRUE, size = 1, span = 0.4, alpha = 0.18) +
  theme_classic(base_size = 5) +
  labs(x = "CytoTRACE2_Relative", y = "ARHGDIB") +
  scale_color_manual(values = tumor_type_color) +
  scale_fill_manual(values = tumor_type_color) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
        axis.line = element_blank(),axis.text = element_text(color = "black", size = 20),
        axis.title = element_text(color = "black", size = 24),
        axis.ticks = element_line(color = "black", linewidth = 1.8),
        axis.ticks.length = unit(0.25, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 22),
        legend.position = "right")
function_plot(filename_prefix = str_glue("{output_dir}/hub_gene"), width = 100, height = 150)

# data process ------------------------------------------------------------
scRNA_scTour_stat <- scRNA_scTour_meta %>% 
  mutate(time_bin = cut(
    ptime,
    breaks = seq(0,1,length.out = 51),
    include.lowest = T,
    right = F
  ))

total_counts <- scRNA_scTour_stat %>%
  dplyr::group_by(tumor_status) %>%
  dplyr::summarise(total_n = n())

scRNA_sctour_result <- scRNA_scTour_stat %>% 
  dplyr::group_by(time_bin, tumor_status) %>% 
  dplyr::summarise(n = n(),.groups = "drop") %>% 
  dplyr::group_by(time_bin) %>% 
  dplyr::mutate(
    total_in_bin = sum(n),
    prop_in_bin = n/total_in_bin
  ) %>% 
  ungroup() %>% 
  left_join(total_counts,by = "tumor_status") %>%
  dplyr::mutate(proportion = prop_in_bin/total_n) %>% 
  dplyr::group_by(time_bin) %>% 
  dplyr::mutate(final_normalized = proportion/sum(proportion)) %>% 
  ungroup()



# result plot -------------------------------------------------------------
ggplot(scRNA_sctour_result, aes(x = time_bin, stratum = tumor_status, alluvium = tumor_status, y = final_normalized, fill = tumor_status)) +
  geom_stratum(width = 0.5, col=NA) +
  geom_alluvium(width = 0.4, alpha = 0.4) +
  scale_fill_manual(values = tumor_status_color) +
  labs(x = '', y = 'Proportion of Cell Type(%)')+
  theme_classic()+ nrc_theme + NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_scTour_stat"), width = 220, height = 80)






