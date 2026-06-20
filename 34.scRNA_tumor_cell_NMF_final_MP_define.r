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
library("ggrastr")
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

scRNA_tumor_origin <- read_h5(file = str_glue("{all_tumor_dir}/scRNA_cell_type_cluster_only.h5"),
                              assay.name = 'RNA', 
                              target.object = 'seurat')
scRNA_tumor_meta <- scRNA_tumor_origin@meta.data
all_gene <- rownames(scRNA_tumor_origin)

group_data <- read.xlsx(group_data_dir,sheetName = "scRNA_seq")
colnames(group_data)[2] <- "sample_ID"


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


# tumor cell character plot -----------------------------------------------
## UMAP plot
DimPlot(scRNA_tumor_origin, reduction = "umap",label = FALSE,raster=FALSE,cols = tumor_status_color,group.by = "tumor_status") +
  cluster_theme + NoLegend() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),axis.ticks.y = element_blank())
function_plot(filename_prefix = str_glue("{output_dir}/tumor_cell_UMAP_status"), width = 400, height = 400)

DimPlot(scRNA_tumor_origin, reduction = "umap",label = FALSE,raster=FALSE,cols = tissue_type_color,group.by = "tissue_type") +
  cluster_theme + NoLegend() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),axis.ticks.y = element_blank())
function_plot(filename_prefix = str_glue("{output_dir}/tumor_cell_UMAP_tissue"), width = 400, height = 400)

DimPlot(scRNA_tumor_origin, reduction = "umap",label = FALSE,raster=FALSE,cols = study_database_color,group.by = "study_ID") +
  cluster_theme + NoLegend() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),axis.ticks.y = element_blank())
function_plot(filename_prefix = str_glue("{output_dir}/tumor_cell_UMAP_study"), width = 400, height = 400)

DimPlot(scRNA_tumor_origin, reduction = "umap",label = FALSE,raster=FALSE,cols = tumor_type_color,group.by = "tumor_code") +
  cluster_theme + NoLegend() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),axis.ticks.y = element_blank())
function_plot(filename_prefix = str_glue("{output_dir}/tumor_cell_UMAP_code"), width = 400, height = 400)

FeatureDimPlot(srt = scRNA_tumor_origin, features = "KRT18",  reduction = "UMAP", theme_use = "theme_blank",raster = F)
# a <- readRDS(" ")

## bar plot
bar_meta_data <- scRNA_tumor_origin@meta.data
bar_meta_data_tumor_code <- as.data.frame(table(bar_meta_data$tumor_code))
colnames(bar_meta_data_tumor_code) <- c("tumor_code","num")
ggplot(bar_meta_data_tumor_code, aes(reorder(tumor_code, - num) ,num, fill = tumor_code)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "", y = "Number of cells") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = tumor_type_color) +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_tumor_cell_num_code"), width = 250, height = 160)
bar_meta_data_tissue <- as.data.frame(table(bar_meta_data$tissue_type))
colnames(bar_meta_data_tissue) <- c("tissue_type","num")
ggplot(bar_meta_data_tissue, aes(reorder(tissue_type, - num) ,num, fill = tissue_type)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "", y = "Number of cells") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = tissue_type_color) +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_tumor_cell_num_tissue"), width = 250, height = 180)

write_csv(bar_meta_data_tumor_code,str_glue("{output_dir}/scRNA_tumor_cell_num_code.csv"))
write_csv(bar_meta_data_tissue,str_glue("{output_dir}/scRNA_tumor_cell_num_tissue.csv"))

# new MP define ----------------------------------------------------------------
program_cluster_df <- program_cluster_df %>%
  mutate(final_MP = case_when(
    MetaProgramName %in% c("IMP31", "IMP32", "IMP28", "IMP10") ~ "Cell.Cycle.G2M",
    MetaProgramName %in% c("IMP35") ~ "Cell.Cycle.G1S",
    MetaProgramName %in% c("IMP29", "IMP12") ~ "Cell.Cycle.HMG",
    MetaProgramName %in% c("IMP16") ~ "Chromatin",
    MetaProgramName %in% c("IMP17") ~ "Cell.Cycle.single.nucleus",
    MetaProgramName %in% c("IMP26", "IMP15") ~ "Stress",
    MetaProgramName %in% c("IMP22") ~ "Hypoxia",
    MetaProgramName %in% c("IMP33") ~ "Stress.in.vitro",
    MetaProgramName %in% c("IMP24") ~ "Protein.maturation",
    MetaProgramName %in% c("IMP9") ~ "EMT.I",
    MetaProgramName %in% c("IMP7") ~ "EMT.II",
    MetaProgramName %in% c("IMP44") ~ "EMT.VI",
    MetaProgramName %in% c("IMP21", "IMP27") ~ "MHC",
    MetaProgramName %in% c("IMP37", "IMP39") ~ "Epithelial.Senescence",
    MetaProgramName %in% c("IMP38") ~ "MYC",
    MetaProgramName %in% c("IMP14") ~ "Respriration",
    MetaProgramName %in% c("IMP34") ~ "Secreted.I",
    MetaProgramName %in% c("IMP4", "IMP48") ~ "Secreted.II",
    MetaProgramName %in% c("IMP49") ~ "Cilia.I",
    MetaProgramName %in% c("IMP45") ~ "Cilia.II",
    MetaProgramName %in% c("IMP40", "IMP42") ~ "PDAC.classcial",
    MetaProgramName %in% c("IMP5") ~ "Alveolar",
    MetaProgramName %in% c("IMP23") ~ "Colon.related",
    MetaProgramName %in% c("IMP18") ~ "PDAC.related",
    MetaProgramName %in% c("IMP20", "IMP25", "IMP36") ~ "Lymphocyte.Activation",
    MetaProgramName %in% c("IMP41") ~ "Cell.Signaling",
    MetaProgramName %in% c("IMP43") ~ "Coagulation",
    MetaProgramName %in% c("IMP1") ~ "Neuron",
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

jaccard_mat_save <- as.data.frame(jaccard_mat_filter) %>% 
  rownames_to_column(var = "sample_MP")
write_csv(jaccard_mat_save,str_glue("{output_dir}/jaccard_mat_save.csv"))

## jar cor cluster
dist_mat <- as.dist(1 - jaccard_mat_filter)
hc <- hclust(dist_mat, method = "ward.D2")

## data plot
program_cluster_df_anno <- data.frame(row.names = rownames(program_cluster_df),
                                      final_MP=program_cluster_df$final_MP)
ano_col <- program_cluster_df
ano_col$sample_ID <- sub("_[^_]+$", "", ano_col$Program)
ano_col <- left_join(ano_col,group_data,by = "sample_ID")
ano_col <- ano_col[,c("study_ID","tumor_code","tissue_type","tumor_status")]
pdf(str_glue("{output_dir}/NMF_final_heatmap.pdf"),width = 15,height = 11)
ComplexHeatmap::pheatmap(
    jaccard_mat_filter,
    annotation_row = program_cluster_df_anno,
    annotation_col = ano_col,
    show_colnames = FALSE,
    show_rownames = FALSE,
    treeheight_row = 20,
    treeheight_col = 20,
    cluster_rows = F,
    cluster_cols = F,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    color = custom_magma,
    annotation_colors = list(final_MP = c(final_MP),tumor_code = c(tumor_type_color),
                             tissue_type = c(tissue_type_color),tumor_status = c(tumor_status_color),study_ID = c(study_database_color)),
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

## result plot
X_order <- c("Cell.Cycle.G2M","Cell.Cycle.G1S","Cell.Cycle.HMG","Chromatin","Cell.Cycle.single.nucleus","Stress","Hypoxia",
             "Stress.in.vitro","Protein.maturation","EMT.I","EMT.II","EMT.VI","MHC","Epithelial.Senescence","MYC",
             "Respriration","Secreted.I","Secreted.II","Cilia.I","Cilia.II","PDAC.classcial","Alveolar","Colon.related",
             "PDAC.related","Lymphocyte.Activation","Cell.Signaling","Coagulation","Neuron")
jaccard_result$MP_self_name <- factor(jaccard_result$MP_self_name, levels = X_order)

jaccard_result <- jaccard_result %>% 
  mutate(hub_cor = ifelse(jaccard_similarity_score == 0, " ", jaccard_similarity_score))
jaccard_result$hub_cor <- as.numeric(jaccard_result$hub_cor)
ggplot(jaccard_result, aes(x = MP_self_name, y = MP_3CA_name, fill = jaccard_similarity_score)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(jaccard_similarity_score >= 0.05, sprintf("%.2f", jaccard_similarity_score)," ")), vjust = 0.5, size = 2, fontface = "bold") +
  scale_fill_gradient2(low = "#6D9EC1",mid = '#D9D9D9', high = "#E46726") +
  theme_minimal()+
   scale_y_discrete(limits = colnames(MP_3CA_hub)) +
  theme(plot.title = element_text(face = "bold",size = 12,color="black",hjust = 0.5),
        axis.title = element_text(face = "bold",size = 12,color ="black"), 
        axis.text = element_text(face = "bold",size= 16,color = "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1 ),
        panel.grid=element_blank(),
        legend.text = element_text(face = "bold",size= 12),
        legend.title= element_text(face = "bold",size= 12)
  ) + labs(fill = "pearson") + xlab("") + ylab("")
function_plot(filename_prefix = str_glue("{output_dir}/MP_cor"), width = 380, height = 550)

write_csv(jaccard_result, str_glue("{output_dir}/jaccard_result.csv"))

## hub plot
hub_3CA <- c("MP1.Cell.Cycle.G2.M","MP2.Cell.Cycle.G1.S","MP3.Cell.Cylce.HMG.rich","MP4.Chromatin","MP5.Cell.cycle.single.nucleus",
             "MP6.Stress.1","MP7.Hypoxia","MP8.Stress.in.vitro","MP12.Protein.maturation","MP15.EMT.I","MP16.EMT.II","MP19.EMT.V",
             "MP22.Interferon.MHC.II.I","MP25.Epithelial.Senescence","MP27.MYC","MP30.Respiration.1","MP33.Secreted.I","MP34.Secreted.II",
             "MP35.Cilia.1","MP36.Cilia.2","MP43.PDAC.classical","MP44.Alveolar","MP47.Colon.related","MP59.PDAC.related.1","MP60.PDAC.related.2")
jaccard_result_filter <- jaccard_result %>% 
  filter(MP_3CA_name %in% hub_3CA)
jaccard_result_filter$MP_self_name <- factor(jaccard_result_filter$MP_self_name, levels = X_order)
jaccard_result_filter$MP_3CA_name <- as.character(jaccard_result_filter$MP_3CA_name)
ggplot(jaccard_result_filter, aes(x = MP_self_name, y = MP_3CA_name, fill = jaccard_similarity_score)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", jaccard_similarity_score)), vjust = 0.5, size = 2, fontface = "bold") +
  scale_fill_gradient2(low = "#6D9EC1",mid = '#D9D9D9', high = "#E46726") +
  theme_minimal()+
  scale_y_discrete(limits = hub_3CA) +
  theme(plot.title = element_text(face = "bold",size = 12,color="black",hjust = 0.5),
        axis.title = element_text(face = "bold",size = 12,color ="black"), 
        axis.text = element_text(face = "bold",size= 16,color = "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1 ),
        panel.grid=element_blank(),
        legend.text = element_text(face = "bold",size= 12),
        legend.title= element_text(face = "bold",size= 12)
  ) + labs(fill = "pearson") + xlab("") + ylab("")
function_plot(filename_prefix = str_glue("{output_dir}/MP_cor_hub"), width = 420, height = 360)


# MP stat ----------------------------------------------------------------
## data
MP_per_sample <- data.frame(
  MP_name = as.character(),
  percent_tumor = as.numeric(),
  percent_tissue = as.numeric(),
  all_sample = as.numeric()
)

## stat
program_cluster_combine <- program_cluster_df
program_cluster_combine$sample_ID <- sub("_[^_]+$", "", program_cluster_combine$Program)
program_cluster_combine <- left_join(program_cluster_combine,group_data,by = "sample_ID")

for(MP in unique(program_cluster_combine$final_MP)){
  # MP <- "Cilia"
  
  current_result <- program_cluster_combine %>% 
    filter(final_MP %in% MP)
  MP_tumor <- length(unique(current_result$tumor_code)) / 18
  MP_tissue <- length(unique(current_result$tissue_type)) / 21
  MP_all <- length(current_result$sample_ID)/1320
  
  ## data rbind
  MP_per_current <- data.frame(
    MP_name = MP,
    percent_tumor = MP_tumor,
    percent_tissue = MP_tissue,
    all_sample = MP_all
  )
  MP_per_sample <- rbind(MP_per_sample,MP_per_current)
  
}


## data plot
ggplot(MP_per_sample, aes(reorder(MP_name, - percent_tumor) ,percent_tumor, fill = MP_name)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "", y = "Number of cells") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = final_MP) +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/MP_per_tumor"), width = 320, height = 284)
ggplot(MP_per_sample, aes(reorder(MP_name, - percent_tissue) ,percent_tissue, fill = MP_name)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "", y = "Number of cells") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = final_MP) +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/MP_per_tissue"), width = 320, height = 280)

## Primary to Metastasis MP plot
ggsankey_plot <- program_cluster_combine %>%
  dplyr::select(tumor_code, final_MP, metastasis_tissue)
ggsankey_plot$metastasis_tissue[ggsankey_plot$metastasis_tissue == "None"] <- NA
  # filter(tumor_status == "Metastasis") %>%
ggsankey_plot <- ggsankey_plot %>% 
  make_long(tumor_code, final_MP, metastasis_tissue)
ggplot(ggsankey_plot, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_manual(values = c(final_MP,tumor_type_color,tissue_type_color)) +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",plot.title = element_text(hjust = .5))
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_ggsankey_plot"), width = 300, height = 300)


# scRNA addmodulescore -----------------------------------------------------
## MP socre
cells_rankings <- AUCell_buildRankings(scRNA_tumor_origin@assays$RNA@data,splitByBlocks=TRUE)
cells_AUC <- AUCell_calcAUC(meta_program_top50, cells_rankings,
                            aucMaxRank=nrow(cells_rankings)*0.1)
auc_score <- as.data.frame(t(getAUC(cells_AUC)[names(meta_program_top50),]))
scRNA_tumor_origin <- AddMetaData(scRNA_tumor_origin,auc_score)
addmodulescore_meta <- scRNA_tumor_origin@meta.data

addmodulescore_meta <- read_csv(str_glue("{output_dir}/addmodulescore_meta.csv"))
## tumor status plot
### violin plot
diff_status <- addmodulescore_meta[,c("tumor_status",names(meta_program_top50))]
diff_status <- melt(diff_status)
colnames(diff_status) <- c("tumor_status","MP_name","score")
ggviolin(diff_status, 
         x="MP_name", y="score", 
         width = 1.5,color = "black",
         fill="tumor_status",
         xlab = F,
         add = "boxplot",
         add.params = list(outlier.shape = NA),
         bxp.errorbar=T, 
         bxp.errorbar.width=0.05,
         size=0.5,
         palette = tumor_status_color, 
         legend = NULL) + nrc_theme +
  stat_compare_means(aes(group=tumor_status), method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", " ")), label = "p.signif") + nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/addmodulescore_diff_status"), width = 300, height = 200)
### boxplot
diff_status_boxplot <- addmodulescore_meta[,c("sample_ID","tumor_status",names(meta_program_top50))]
diff_status_boxplot <- melt(diff_status_boxplot)
colnames(diff_status_boxplot) <- c("sample_ID","tumor_status","MP_name","score")
diff_status_boxplot <- diff_status_boxplot %>%
  dplyr::group_by(sample_ID,tumor_status, MP_name) %>%
  dplyr::summarise(mean_sample_score = mean(score,na.rm = T),.groups = 'drop') %>% 
  # dplyr::group_by(tumor_status,MP_name) %>% 
  # dplyr::summarise(mean_final_socre = mean(mean_sample_score,na.rm = T),sample_count = n()) %>% 
  ungroup()
ggplot(diff_status_boxplot,aes(x = MP_name,y = mean_sample_score)) +
  labs(y = 'Proportion',x = NULL,title = NULL) +
  geom_boxplot(aes(fill = tumor_status),position = position_dodge(0.7),width = 0.5,outlier.alpha = 0) +
  scale_fill_manual(values = tumor_status_color) +
  theme_bw() + nrc_theme +
  stat_compare_means(aes(group = tumor_status),label = "p.signif",method = "wilcox.test",hide.ns = T)

## meta tissue type plot
diff_tissue <- addmodulescore_meta %>% 
  filter(tumor_status == "Metastasis")
diff_tissue <- diff_tissue[,c("tissue_type",names(meta_program_top50))]
diff_tissue <- melt(diff_tissue)
colnames(diff_tissue) <- c("tissue_type","MP_name","score")
ggviolin(diff_tissue, 
         x="MP_name", y="score", 
         width = 1.5,color = "black",
         fill="tissue_type",
         xlab = F,
         add = "boxplot",
         add.params = list(outlier.shape = NA),
         bxp.errorbar=T, 
         bxp.errorbar.width=0.02,
         size=0.5,
         palette = tissue_type_color, 
         legend = NULL) + nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/addmodulescore_diff_tissue_type"), width = 600, height = 200)


# scRNA MP score plot -----------------------------------------------------

## UMAP
for(MP in names(meta_program_top50)){
  FeatureDimPlot(srt = scRNA_tumor_origin, features = c(MP),  reduction = "UMAP", theme_use = "theme_blank",raster = F)
  function_plot(filename_prefix = str_glue("{output_dir}/UMAP_MP_{MP}_socre"), width = 200, height = 200) 
}

## boxplot
for(MP in names(meta_program_top50)){
  MP <- "Secreted.II"
  addmodulescore_meta$MP_socre <- addmodulescore_meta[[MP]]
  ggplot(addmodulescore_meta,aes(x = reorder(leiden,-MP_socre),y = MP_socre))+ 
    geom_boxplot(outlier.alpha = 0)+ 
    theme_classic() +
    # scale_fill_manual(values = final_MP) +
    nrc_theme
  function_plot(filename_prefix = str_glue("{output_dir}/boxplot_MP_{MP}_socre"), width = 300, height = 200) 
}


# scRNA MP define ---------------------------------------------------------
## MP ano
new.cluster.ids <- c(
  "Other", # 0
  "Other", # 1
  "Other", # 2
  "Other", # 3
  "Secreted.I", # 4
  "Secreted.II", # 5
  "Cell.Cycle.HMG", # 6
  "Cell.Cycle.HMG", # 7
  "Other", # 8
  "Respriration", # 9
  "Other", # 10
  "Other", # 11
  "Other", # 12
  "Cell.Signaling", # 13
  "Neuron", # 14
  "Neuron", # 15
  "Other", # 16
  "Respriration", # 17
  "EMT.I", # 18
  "Other", # 19
  "Hypoxia", # 20
  "Cilia.I", # 21
  "Other", # 22
  "Cell.Signaling", # 23
  "Chromatin", # 24
  "Other", # 25
  "Epithelial.Senescence", # 26
  "Other", # 27
  "Other", # 28
  "Coagulation", # 29
  "Other" # 30
  
)



# new.cluster.ids <- c(
#   "Other", # 0
#   "Other", # 1
#   "Other", # 2
#   "Other", # 3
#   "Secreted.I", # 4
#   "Secreted.II", # 5
#   "Cell.Cycle.G2M", # 6
#   "Cell.Cycle.G1S", # 7
#   "Lymphocyte.Activation", # 8
#   "Other", # 9
#   "EMT.VI", # 10
#   "Other", # 11
#   "Alveolar", # 12
#   "Cell.Signaling", # 13
#   "Neuron", # 14
#   "Neuron", # 15
#   "Stress", # 16
#   "MYC", # 17
#   "EMT.I", # 18
#   "Other", # 19
#   "Hypoxia", # 20
#   "Cilia.I", # 21
#   "Other", # 22
#   "Cell.Signaling", # 23
#   "Chromatin", # 24
#   "MHC", # 25
#   "Epithelial.Senescence", # 26
#   "Cell.Cycle.single.nucleus", # 27
#   "Stress.in.vitro", # 28
#   "Coagulation", # 29
#   "Other" # 30
#   
# )
scRNA_ano <- scRNA_tumor_origin
Idents(scRNA_ano) <- scRNA_ano$leiden
names(new.cluster.ids) <- levels(scRNA_ano)
scRNA_ano <- RenameIdents(scRNA_ano, new.cluster.ids)

## join cell type in meta data
new_cell <- rownames_to_column(as.data.frame(new.cluster.ids),var = "leiden")
colnames(new_cell) <- c("leiden","final_MP")
meta_data <- scRNA_ano@meta.data %>% 
  rownames_to_column(var = 'id') %>% 
  left_join(.,new_cell,by = 'leiden') %>% 
  column_to_rownames(var = 'id')
scRNA_ano@meta.data <- meta_data

## scRNA MP plot
DimPlot(scRNA_ano, reduction = "umap",label = FALSE,raster=FALSE,group.by = "final_MP",cols = final_MP) + 
  cluster_theme + NoLegend() + title("") +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),axis.ticks.y = element_blank())
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_MP_plot"), width = 400, height = 410) 


## scRNA MP stat
scRNA_MP_meta <- scRNA_ano@meta.data %>% 
  # filter(final_MP %in% c("Cell.Cycle.G2M","Cell.Cycle.single.nucleus","EMT.I","Neuron","Secreted.II")) %>%
  # filter(final_MP != "Other") %>% 
  dplyr::select(c("sample_ID","tumor_status","final_MP")) %>% 
  dplyr::group_by(sample_ID,tumor_status) %>% 
  dplyr::mutate(total_cells = n()) %>% 
  dplyr::group_by(sample_ID,tumor_status,final_MP,total_cells) %>% 
  dplyr::summarise(count = n(),.groups = 'drop') %>% 
  dplyr::mutate(percentage = (count/total_cells) * 100) %>% 
  dplyr::select(sample_ID,tumor_status,final_MP,percentage)

ggplot(scRNA_MP_meta,aes(x = final_MP,y = percentage)) +
  labs(y = 'Proportion',x = NULL,title = NULL) +
  geom_boxplot(aes(fill = tumor_status),position = position_dodge(0.7),width = 0.5,outlier.alpha = 0) +
  scale_fill_manual(values = tumor_status_color) +
  theme_bw() + nrc_theme +
  stat_compare_means(aes(group = tumor_status),label = "p.signif",method = "wilcox.test",hide.ns = T)


## OR meta
all_meta <- scRNA_tumor_origin@meta.data
scRNA_MP_OR_data <- program_cluster_combine
scRNA_MP_OR_data$MP_code <- str_c(scRNA_MP_OR_data$sample_ID,scRNA_MP_OR_data$final_MP,sep = '_')
scRNA_MP_OR_data <- scRNA_MP_OR_data[!duplicated(scRNA_MP_OR_data$MP_code),]
for(status in c("Tumor","Metastasis")){
  # status <- "Tumor"
  for(MP in X_order){
    # MP <- "Cell.Cycle.G2M"
    current_MP <- scRNA_MP_OR_data %>% 
      filter(tumor_status == status & final_MP == MP)
    MP_per <- nrow(current_MP) / length(unique(all_meta[all_meta$tumor_status == status,]$sample_ID))
    
    ## data collect
    current_result <- data.frame(
      MP_name = MP,
      tumor_status = status,
      MP_percentage = MP_per * 100
    )
    all_MP_per_result <- rbind(all_MP_per_result,current_result)
  }
}
all_MP_per_result$FC <- NA
for(line in 1:nrow(all_MP_per_result)){
  # line <- 1
  all_MP_per_result[line,]$FC <- all_MP_per_result[all_MP_per_result$MP_name == all_MP_per_result[line,]$MP_name & all_MP_per_result$tumor_status == "Metastasis",]$MP_percentage /
    all_MP_per_result[all_MP_per_result$MP_name == all_MP_per_result[line,]$MP_name & all_MP_per_result$tumor_status == "Tumor",]$MP_percentage
                                
  
}
# all_meta <- scRNA_ano@meta.data
for(tissue in c("Brain","Liver","Lymph_node","Other")){
  tissue <- "Brain"
  for(status in c("Tumor","Metastasis")){
    status <- "Metastasis"
    for(MP in X_order){
      MP <- "Cell.Cycle.G2M"
      
      if(tissue == "Other"){
        tumor_list <- scRNA_MP_OR_data %>% 
          filter(metastasis_tissue %in% c("Lung","Ovary","Vaginal","Peritoneal","Bone_marrow"))
      }else{
        tumor_list <- scRNA_MP_OR_data %>% 
          filter(metastasis_tissue == tissue)
      }
      
      if(status == "Tumor"){
        current_MP <- scRNA_MP_OR_data %>% 
          filter(tumor_status == status & final_MP == MP) %>% 
          filter(tumor_code %in% tumor_list$tumor_code)
      }else if(status == "Metastasis"){
        current_MP <- scRNA_MP_OR_data %>% 
          filter(tumor_status == status & final_MP == MP) %>% 
          filter(metastasis_tissue == tissue)
      }

      if(status == "Tumor"){
        current_meta <- all_meta %>% 
          filter(tumor_code %in% tumor_list$tumor_code)
      }else if(status == "Metastasis"){
        current_meta <- all_meta %>% 
          filter(metastasis_tissue == tissue)
      }

      MP_per <- nrow(current_MP) / length(unique(current_meta[current_meta$tumor_status == status,]$sample_ID))
      
      ## data collect
      current_result <- data.frame(
        meta_tissue = tissue,
        MP_name = MP,
        tumor_status = status,
        MP_percentage = MP_per * 100
      )
      tissue_MP_per_result <- rbind(tissue_MP_per_result,current_result)
      
    }
  }
}
write_csv(tissue_MP_per_result,str_glue("{output_dir}/tissue_MP_per_result.csv"))

## OR result plot
### all plot
write_csv(all_MP_per_result,str_glue("{output_dir}/all_MP_per_result.csv"))
ggplot(all_MP_per_result, aes(reorder(MP_name,-FC) ,MP_percentage, fill = tumor_status)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  theme_bw() +
  labs(x = "", y = "Percentage of MP") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = tumor_status_color) +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_status_MP_all"), width = 330, height = 220)
### diff tissue plot
for(tissue in c("Brain","Liver","Lymph_node","Other")){
  # tissue <- "Lymph_node"
  current_plot <- tissue_MP_per_result %>% 
    filter(meta_tissue == tissue)
  ggplot(current_plot, aes(reorder(MP_name,-MP_percentage) ,MP_percentage, fill = tumor_status)) + 
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    theme_bw() +
    labs(x = "", y = "Percentage of MP") +
    scale_y_continuous(expand = c(0.01,0))+
    scale_fill_manual(values = tumor_status_color) +
    nrc_theme
  function_plot(filename_prefix = str_glue("{output_dir}/scRNA_status_MP_{tissue}"), width = 330, height = 220)
}

## meta tissue
tissue_MP_per_result <- data.frame(
  meta_tissue = as.character(),
  MP_name = as.character(),
  # tumor_status = as.character(),
  MP_percentage = as.numeric()
)
# X_order <- c("Cell.Cycle.G1S",
#              "Protein.maturation","EMT.II","EMT.VI","MHC","Epithelial.Senescence","MYC",
#              "Secreted.I","Secreted.II","Cilia.I","Cilia.II","PDAC.classcial","Alveolar","Colon.related",
#              "PDAC.related","Cell.Signaling","Coagulation",)
metastasis_MP <- c("Respriration","Cell.Cycle.single.nucleus","Stress.in.vitro","Cell.Cycle.HMG","Hypoxia","Chromatin",
                   "Neuron","EMT.I","Lymphocyte.Activation","Stress","Cell.Cycle.G2M","MHC")
metastasis_MP <- c(c("Cell.Cycle.G2M","Cell.Cycle.G1S","Cell.Cycle.HMG","Chromatin","Cell.Cycle.single.nucleus","Stress","Hypoxia",
                     "Stress.in.vitro","Protein.maturation","EMT.I","EMT.II","EMT.VI","MHC","Epithelial.Senescence","MYC",
                     "Respriration","Secreted.I","Secreted.II","Cilia.I","Cilia.II","PDAC.classcial","Alveolar","Colon.related",
                     "PDAC.related","Lymphocyte.Activation","Cell.Signaling","Coagulation","Neuron"))


metastasis_tissue_OR <- scRNA_MP_OR_data %>% 
  filter(tumor_status == "Metastasis")
for(tissue in unique(metastasis_tissue_OR$metastasis_tissue)){
  # tissue <- "Brain"
  for(MP in metastasis_MP){
    # MP <- "Respriration"
    metastasis_tissue_OR_current <- metastasis_tissue_OR %>% 
      filter(metastasis_tissue == tissue & final_MP == MP)
    
    MP_tissue_per <- length(unique(metastasis_tissue_OR_current$sample_ID)) / 
      length(unique(metastasis_tissue_OR[metastasis_tissue_OR$metastasis_tissue == tissue,]$sample_ID))
    
    ## data combine
    current_result <- data.frame(
      meta_tissue = tissue,
      MP_name = MP,
      MP_percentage = MP_tissue_per
    )
    tissue_MP_per_result <- rbind(tissue_MP_per_result,current_result)
  }
  
}
write_csv(tissue_MP_per_result,str_glue("{output_dir}/tissue_MP_per_result.csv"))
tissue_MP_per_result[tissue_MP_per_result$MP_percentage == 0,]$MP_percentage <- 0.01
ggplot(tissue_MP_per_result, aes(meta_tissue ,MP_percentage, fill = MP_name)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  theme_bw() +
  labs(x = "", y = "Percentage of MP") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = final_MP) +
  nrc_theme
# function_plot(filename_prefix = str_glue("{output_dir}/scRNA_status_MP_per_meta_tissue"), width = 350, height = 150)
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_status_MP_per_meta_tissue_all"), width = 450, height = 150)

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
  for(MP in metastasis_MP){
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
write_csv(primary_tissue_MP_per_result,str_glue("{output_dir}/primary_tissue_MP_per_result.csv"))
primary_tissue_MP_per_result[primary_tissue_MP_per_result$MP_percentage == 0,]$MP_percentage <- 0.01
# primary_tissue_MP_per_result <- primary_tissue_MP_per_result[primary_tissue_MP_per_result$MP_percentage != 0,]
ggplot(primary_tissue_MP_per_result, aes(tumor_code ,MP_percentage, fill = MP_name)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  theme_bw() +
  labs(x = "", y = "Percentage of MP") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = final_MP) +
  nrc_theme
# function_plot(filename_prefix = str_glue("{output_dir}/scRNA_status_MP_per_primary_tissue"), width = 380, height = 170)
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_status_MP_per_primary_tissue_all"), width = 500, height = 170)

write_csv(tissue_MP_per_result,str_glue("{output_dir}/scRNA_status_MP_per_meta_tissue_all.csv"))
write_csv(primary_tissue_MP_per_result,str_glue("{output_dir}/scRNA_status_MP_per_primary_tissue_all.csv"))

# diff MP cor -------------------------------------------------------------
# addmodulescore_meta <- read_csv(str_glue("{output_dir}/addmodulescore_meta.csv"))
# meta_program_top50 <- readRDS(str_glue("{output_dir}/meta_program_top50.RData"))
addmodulescore_meta_hub <- addmodulescore_meta[,names(meta_program_top50)]
MP_cor_all_result <- rcorr(as.matrix(addmodulescore_meta_hub),type = "pearson")
pdf(height = 10,width = 12,file = str_glue("{output_dir}/MP_cor_all.pdf"))
corrplot(MP_cor_all_result_1$r,type = "upper",
         # order = "hclust",
         tl.col = "black",
         tl.srt = 90)
dev.off()


# NMF top gene enrichment -------------------------------------------------
meta_program_top50 <- readRDS(str_glue("{output_dir}/meta_program_top50.RData"))
enrich_GO_BP <- lapply(meta_program_top50, function(program){runGSEA(program, universe=all_gene,
                                                                  category = "C5", subcategory = "GO:BP")})
enrich_GO_CC <- lapply(meta_program_top50, function(program){runGSEA(program, universe=all_gene,
                                                                     category = "C5", subcategory = "GO:CC")})
enrich_GO_MF <- lapply(meta_program_top50, function(program){runGSEA(program, universe=all_gene,
                                                                     category = "C5", subcategory = "GO:MF")})
enrich_Hallmarker <- lapply(meta_program_top50, function(program){runGSEA(program, universe=all_gene,
                                                                          category = "H")})
enrich_KEGG <- lapply(meta_program_top50, function(program){runGSEA(program, universe=all_gene,
                                                                    category = "C8")})
## enrich save
all_data <- data.frame(
  pathway = as.character(),
  pval = as.numeric(),
  padj = as.numeric(),
  foldEnrichment = as.numeric(),
  overlap = as.numeric(),
  size = as.numeric(),
  overlapGenes = as.character(),
  MP = as.character(),
  type = as.character()
)

hall_data <- data.frame(
  pathway = as.character(),
  pval = as.numeric(),
  padj = as.numeric(),
  foldEnrichment = as.numeric(),
  overlap = as.numeric(),
  size = as.numeric(),
  overlapGenes = as.character(),
  MP = as.character(),
  type = as.character()
)
go_data <- data.frame(
  pathway = as.character(),
  pval = as.numeric(),
  padj = as.numeric(),
  foldEnrichment = as.numeric(),
  overlap = as.numeric(),
  size = as.numeric(),
  overlapGenes = as.character(),
  MP = as.character(),
  type = as.character()
)
for(MP_name in names(meta_program_top50)){
  # MP_name <- "Alveolar"
  # hall_current <- head(enrich_Hallmarker[[MP_name]],10)
  hall_current <- enrich_Hallmarker[[MP_name]]
  hall_current$MP <- MP_name
  hall_current$type <- "HALL"
  hall_data <- rbind(hall_data,hall_current)
  
  # go_current_BP <- head(enrich_GO_BP[[MP_name]],10)
  go_current_BP <- enrich_GO_BP[[MP_name]]
  # go_current_CC <- head(enrich_GO_CC[[MP_name]],10)
  go_current_CC <- enrich_GO_CC[[MP_name]]
  # go_current_MF <- head(enrich_GO_MF[[MP_name]],10)
  go_current_MF <- enrich_GO_MF[[MP_name]]
  go_current_BP$MP <- MP_name
  go_current_BP$type <- "GO_BP"
  go_current_CC$MP <- MP_name
  go_current_CC$type <- "GO_CC"
  go_current_MF$MP <- MP_name
  go_current_MF$type <- "GO_MF"
  go_data <- rbind(go_data,go_current_BP)
  go_data <- rbind(go_data,go_current_CC)
  go_data <- rbind(go_data,go_current_MF)
  
  ## all combine
  all_data <- rbind(all_data,hall_current)
  all_data <- rbind(all_data,go_current_BP)
  all_data <- rbind(all_data,go_current_CC)
  all_data <- rbind(all_data,go_current_MF)
}

## result plot
# all_data <- read_csv(str_glue("{output_dir}/enrichment_all_pathway.csv"))
hub_pathway <- read_csv(str_glue("{output_dir}/enrichment_all_pathway_self_filter.csv"))
all_data_filter <- all_data %>% 
  filter(pathway %in% unique(hub_pathway$pathway))
all_data_filter$log10p <- -log10(all_data_filter$pval)
ggplot(all_data_filter, aes(x = pathway, y = MP)) +
  geom_point(aes(color = type,size = log10p)) +
  scale_size_continuous(range = c(2,8),breaks = c(0,4,8,12,16,20)) +
  scale_y_discrete(limits = unique(all_data_filter$MP)) +
  scale_x_discrete(limits = unique(hub_pathway$pathway)) +
  nrc_theme + scale_color_manual(values = c("GO_BP"="#FDBF6F","GO_CC"="#B2DF8A","GO_MF"="#A6CEE3","HALL"="#33A02C")) +
  xlab("") + ylab("") 
function_plot(filename_prefix = str_glue("{output_dir}/NMF_MP_enrich_dotplot"), width = 500, height = 350)
ggplot(all_data_filter, aes(x = pathway, y = MP)) +
  geom_point(aes(color = type,size = log10p)) +
  scale_size_continuous(range = c(2,8),breaks = c(0,4,8,12,16,20)) +
  scale_y_discrete(limits = unique(all_data_filter$MP)) +
  scale_x_discrete(limits = unique(hub_pathway$pathway)) +
  nrc_theme + scale_color_manual(values = c("GO_BP"="#FDBF6F","GO_CC"="#B2DF8A","GO_MF"="#A6CEE3","HALL"="#33A02C")) +
  xlab("") + ylab("") + coord_flip()
function_plot(filename_prefix = str_glue("{output_dir}/NMF_MP_enrich_dotplot_transform"), width = 380, height = 500)


# data save ----------------------------------------------------------------
scRNA_ano <- readRDS(str_glue("{output_dir}/scRNA_ano_MP.RData"))
if(T){
  saveRDS(meta_program_top50,str_glue("{output_dir}/meta_program_top50.RData"))
  saveRDS(scRNA_ano,str_glue("{output_dir}/scRNA_ano_MP.RData"))
  write_csv(MP_per_sample,str_glue("{output_dir}/MP_per_sample.csv"))
  write_csv(addmodulescore_meta,str_glue("{output_dir}/addmodulescore_meta.csv"))
  write_csv(hall_data,str_glue("{output_dir}/enrichment_MP_hall.csv"))
  write_csv(go_data,str_glue("{output_dir}/enrichment_MP_go.csv"))
  write_csv(all_data,str_glue("{output_dir}/enrichment_all_pathway.csv"))
  write_csv(all_data,str_glue("{output_dir}/enrichment_all_pathway_self_filter.csv"))
}






