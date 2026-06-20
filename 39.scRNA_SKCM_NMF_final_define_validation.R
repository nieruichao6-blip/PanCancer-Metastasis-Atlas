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

all_MP_per_result <- data.frame(
  MP_name = as.character(),
  tumor_status = as.character(),
  MP_percentage = as.numeric()
)
tissue_MP_per_result <- data.frame(
  meta_tissue = as.character(),
  MP_name = as.character(),
  tumor_status = as.character(),
  MP_percentage = as.numeric()
)

# data dir ----------------------------------------------------------------
MP_3CA_hub_dir <- " "
all_tumor_dir <- " "
group_data_dir <- " "
input_dir <- " "
output_dir <- " "


# data read ----------------------------------------------------------------
MP_3CA_hub <- read.xlsx(str_glue("{MP_3CA_hub_dir}"),sheetName = "cancer_MPs")
jaccard_mat <- readRDS(str_glue("{input_dir}/jaccard_matrix.RData"))
all_programs_top50 <- readRDS(str_glue("{input_dir}/all_programs_top50.RData"))
program_cluster_df <- read_csv(str_glue("{input_dir}/program_cluster_df.csv"))

scRNA_tumor_origin <- read_h5(file = str_glue("{all_tumor_dir}/scRNA_Melanoma_tumor.h5"),
                              assay.name = 'RNA', 
                              target.object = 'seurat')
scRNA_tumor_meta <- scRNA_tumor_origin@meta.data

group_data <- read.xlsx(group_data_dir,sheetName = "scRNA_seq")
colnames(group_data)[2] <- "sample_ID"

# data color ----------------------------------------------------------------
final_MP <- c("Cell.Cycle.G2M"="#e0bc58","Cell.Cycle.G1S"="#64abc0","Cell.Cycle.HMG"="#fab37f","Chromatin"="#e98741",
              "Stress"="#8fc0dc","Hypoxia"="#967568","Stress.in.vitro"="#f2d3ca","Unfolded.Protein.response"="#eebd85",
              "EMT.I"="#82c785","MHC"="#edeaa4","Epithelial.Senescence"="#cdaa9f","MYC"="#794976",
              "Respriration"="#bcacd3","Secreted.I"="#889b5d","Cilia.I"="#4e9592","PDAC.classcial"="#dbad5f",
              "Alveolar"="#64ae79","PDAC.related"="#ac5092","Lymphocyte.Activation"="#dc8e97","Cell.Signaling"="#e3d1db",
              "Coagulation"="#74a893","Neuron"="#ac9141",
              "Cell.Cycle.single.nucleus"="#5ac6e9","EMT"="#ebce8e",
              "Skin.pigmentation"="#e5c06e","Secreted.II"="#7587b1","Translation.initiation"="#c7deef","Colon.related"="#e97371","Other"="#1CB232"
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
                       "Skin"="#204B75","Breast"="#588257","Thyroid"="#B6DB7B","Esophagus"="#4285BF","larynx"="#1CB232")
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
                          "GSE243639"="#6FCDDC")       


# new MP define ----------------------------------------------------------------

program_cluster_df <- program_cluster_df %>%
  mutate(final_MP = case_when(
    MetaProgramName %in% c("IMP23", "IMP3") ~ "Cell.Cycle.G2M",
    MetaProgramName %in% c("IMP34") ~ "Cell.Cycle.G1S",
    MetaProgramName %in% c("IMP14", "IMP10") ~ "Cell.Cycle.HMG",
    MetaProgramName %in% c("IMP25", "IMP30", "IMP38") ~ "Stress",
    MetaProgramName %in% c("IMP5","IMP26") ~ "Hypoxia",
    MetaProgramName %in% c("IMP17") ~ "Unfolded.Protein.response",
    MetaProgramName %in% c("IMP6", "IMP19") ~ "Translation.initiation",
    MetaProgramName %in% c("IMP33","IMP12","IMP16") ~ "EMT",
    MetaProgramName %in% c("IMP27", "IMP9") ~ "MHC",
    # MetaProgramName %in% c("IMP21") ~ "Epithelial.Senescence",
    # MetaProgramName %in% c("IMP31") ~ "MYC",
    MetaProgramName %in% c("IMP31") ~ "Respriration",
    # MetaProgramName %in% c("IMP21","IMP31","IMP36") ~ "Respriration",
    MetaProgramName %in% c("IMP13","IMP18","IMP24","IMP4","IMP7") ~ "Skin.pigmentation",
    # MetaProgramName %in% c("IMP34") ~ "Secreted.I",
    # MetaProgramName %in% c("IMP4", "IMP48") ~ "Secreted.II",
    # MetaProgramName %in% c("IMP49") ~ "Cilia.I",
    # MetaProgramName %in% c("IMP45") ~ "Cilia.II",
    # MetaProgramName %in% c("IMP40", "IMP42") ~ "PDAC.classcial",
    # MetaProgramName %in% c("IMP5") ~ "Alveolar",
    # MetaProgramName %in% c("IMP23") ~ "Colon.related",
    # MetaProgramName %in% c("IMP18") ~ "PDAC.related",
    MetaProgramName %in% c("IMP22", "", "") ~ "Lymphocyte.Activation",
    # MetaProgramName %in% c("IMP41") ~ "Cell.Signaling",
    # MetaProgramName %in% c("IMP43") ~ "Coagulation",
    MetaProgramName %in% c("IMP11","IMP15","IMP20","IMP32") ~ "Neuron",
    TRUE ~ "Other"
  ))



program_cluster_df <- program_cluster_df %>% 
  filter(final_MP != "Other") %>% 
  arrange(final_MP)


# final MP define ----------------------------------------------------------------
## MP gene
meta_program_top50 <- list()
for(meta in unique(program_cluster_df$final_MP)){    
  # meta <- "MP1"
  progs <- program_cluster_df$Program[program_cluster_df$final_MP == meta]
  genes <- unlist(all_programs_top50[,progs])
  # gene_freq <- table(genes)
  # gene_freq_sorted <- sort(gene_freq, decreasing = TRUE)
  # top50 <- names(gene_freq_sorted)[1:50]
  # meta_program_top50[[meta]] <- top50
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


## jar matrix filter
jaccard_mat_filter <- jaccard_mat
jaccard_mat_filter <- jaccard_mat_filter[program_cluster_df$Program,program_cluster_df$Program]
# jaccard_mat_filter <- jaccard_mat_filter %>% 
#   dplyr::select()

## jar cor cluster
dist_mat <- as.dist(1 - jaccard_mat_filter)
hc <- hclust(dist_mat, method = "ward.D2")

## data plot
program_cluster_df_anno <- data.frame(row.names = rownames(program_cluster_df),
                                      final_MP=program_cluster_df$final_MP)
pdf(str_glue("{output_dir}/NMF_final_heatmap.pdf"),width = 13,height = 11)
ComplexHeatmap::pheatmap(
  jaccard_mat_filter,
  annotation_row = program_cluster_df_anno,
  show_colnames = FALSE,
  show_rownames = FALSE,
  treeheight_row = 20,
  treeheight_col = 20,
  cluster_rows = F,
  cluster_cols = F,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  color = custom_magma,
  annotation_colors = list(final_MP = c(final_MP)),
  breaks = seq(0.01, 0.25, length.out = 333)
  
)
dev.off()

# jaccard cor ----------------------------------------------------------------
jaccard_result <- data.frame(
  MP_3CA_name = as.character(),
  MP_self_name = as.character(),
  jaccard_similarity_score = as.numeric()
)

## cor analysis
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


# MP stat -----------------------------------------------------------------
metastasis_MP <- c("Respriration","Cell.Cycle.single.nucleus","Stress.in.vitro","Cell.Cycle.HMG","Hypoxia","Chromatin",
                   "Neuron","EMT","Lymphocyte.Activation","Stress","Cell.Cycle.G2M")

## stat
program_cluster_combine <- program_cluster_df
program_cluster_combine$sample_ID <- sub("_[^_]+$", "", program_cluster_combine$Program)
program_cluster_combine <- left_join(program_cluster_combine,group_data,by = "sample_ID")

scRNA_MP_OR_data <- program_cluster_combine
scRNA_MP_OR_data$MP_code <- str_c(scRNA_MP_OR_data$sample_ID,scRNA_MP_OR_data$final_MP,sep = '_')
scRNA_MP_OR_data <- scRNA_MP_OR_data[!duplicated(scRNA_MP_OR_data$MP_code),]

## primary tissue
primary_tissue_MP_per_result <- data.frame(
  tumor_code = as.character(),
  MP_name = as.character(),
  # tumor_status = as.character(),
  MP_percentage = as.numeric()
)

primary_tissue_OR <- scRNA_MP_OR_data %>% 
  filter(tumor_status == "Tumor")

for(tumor in unique(primary_tissue_OR$tumor_code)){
  # tissue <- "Brain"
  for(MP in unique(primary_tissue_OR$final_MP)){
    # MP <- "Respriration"
    metastasis_tissue_OR_current <- primary_tissue_OR %>% 
      filter(tumor_code == tumor & final_MP == MP)
    
    MP_tissue_per <- length(unique(metastasis_tissue_OR_current$sample_ID)) / 
      length(unique(primary_tissue_OR[primary_tissue_OR$tumor_code == tumor,]$sample_ID))
    
    ## data combine
    current_result <- data.frame(
      tumor_code = tumor,
      MP_name = MP,
      MP_percentage = MP_tissue_per
    )
    primary_tissue_MP_per_result <- rbind(primary_tissue_MP_per_result,current_result)
  }
  
}
primary_tissue_MP_per_result <- primary_tissue_MP_per_result[primary_tissue_MP_per_result$MP_percentage != 0,]
ggplot(primary_tissue_MP_per_result, aes(tumor_code ,MP_percentage, fill = MP_name)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  theme_bw() +
  labs(x = "", y = "Percentage of MP") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = final_MP) +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_status_MP_per_primary_tissue"), width = 180, height = 180)

















