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
library("SCP")
library("Startrac")
library("ggrastr")
library("ROGUE")
(.packages())
# setwd(' ')
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
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
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

# source(" ")

# data dir ----------------------------------------------------------------
tumor_cell_dir <- " "
input_dir <- " "
output_dir <- " "

# data color --------------------------------------------------------------
tumor_type_color <- c("ACC"="#5A7B8F","ANS"="#EE4C97","BRCA"="#3D806F","CESC"="#A08634","CRC"="#F37C95","ESCC"="#608541",
                      "HNSC"="#7D4E57","KIRC"="#BC3C29","LUAD"="#958056","LUSC"="#9FAFA3","NPC"="#6F99AD","PAAD"="#0072B5",
                      "PNET"="#CFC59A","SARC"="#E18727","SKCM"="#FFDC91","GC"="#718DAE","TGCT"="#7876B1","THCA"="#F9AC93","MA"="#FDC44D",
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

all_meta <- data.frame(
  sample_ID = character()
)
metadata_summary <- data.frame(
  tumor_code = as.character(),
  new_cell_type = as.character(),
  clust_prop = as.numeric(),
  origin_cell_type = as.character()
)
all_cell_num_stat <- data.frame(
  cell_type = as.character(),
  tumor_code_stat = as.numeric(),
  tissue_type_stat = as.numeric()
)


# none tumor data process -------------------------------------------------
celltype_file <- list.files(input_dir)

## data plot
for(file in 1:length(celltype_file)){
  # file <- 5
  print(str_glue("############################## {celltype_file[file]} Pre Processing ############################"))
  set.seed(1234)
  
  ## data read
  scRNA_annotation <- readRDS(str_glue("{input_dir}/{celltype_file[file]}/scRNA_annotation.RData"))
  
  ## dir create
  output_file <- str_glue("{output_dir}/{celltype_file[file]}/")
  unlink(output_file, recursive = T, force = T, expand = T)
  dir.create(output_file, recursive = T)
  
  ## size confirm
  size <- 20000/nrow(scRNA_annotation@meta.data)
  
  ## cell type umap plot
  CellDimPlot(scRNA_annotation, reduction = "umap",group.by = "new_cell_type",raster=FALSE,pt.size = size,label.size = 5,
              label = T,label_insitu = T,label_repel = T) + 
    NoLegend() + cluster_theme +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),axis.ticks.y = element_blank()) + labs(title = NULL)
  function_plot(filename_prefix = str_glue("{output_file}/UMAP_annotation_label"), width = 200, height = 200)
  CellDimPlot(scRNA_annotation, reduction = "umap",group.by = "new_cell_type",raster=FALSE,pt.size = size,label = F,label_insitu = T,label_repel = T) + 
    cluster_theme + NoLegend() + 
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),axis.ticks.y = element_blank()) + labs(title = NULL)
  function_plot(filename_prefix = str_glue("{output_file}/UMAP_annotation_nolabel"), width = 200, height = 200)
  CellDimPlot(scRNA_annotation, reduction = "umap",group.by = "new_cell_type",raster=FALSE,pt.size = size,label = T,label_insitu = T,label_repel = T) + 
    cluster_theme +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),axis.ticks.y = element_blank()) + labs(title = NULL)
  function_plot(filename_prefix = str_glue("{output_file}/UMAP_annotation_nolabel_legend"), width = 250, height = 200)
  
  ## tumor status plot
  DimPlot(scRNA_annotation, reduction = "umap",group.by = "tumor_status",raster=FALSE,cols = tumor_status_color) +
    cluster_theme + NoLegend() +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),axis.ticks.y = element_blank()) + labs(title = NULL)
  function_plot(filename_prefix = str_glue("{output_file}/UMAP_tumor_status"), width = 150, height = 150)
  
  ## tissue status plot
  DimPlot(scRNA_annotation, reduction = "umap",group.by = "tissue_type",raster=FALSE,cols = tissue_type_color) +
    cluster_theme + NoLegend() +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),axis.ticks.y = element_blank()) + labs(title = NULL)
  function_plot(filename_prefix = str_glue("{output_file}/UMAP_tissue_type"), width = 150, height = 150)
  
  ## tumor type plot
  DimPlot(scRNA_annotation, reduction = "umap",group.by = "tumor_code",raster=FALSE,cols = tumor_type_color) +
    cluster_theme + NoLegend() +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),axis.ticks.y = element_blank()) + labs(title = NULL)
  function_plot(filename_prefix = str_glue("{output_file}/UMAP_tumor_code"), width = 150, height = 150)
  
  ## study database plot
  DimPlot(scRNA_annotation, reduction = "umap",group.by = "study_ID",raster=FALSE,cols = study_database_color) +
    cluster_theme + NoLegend() +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),axis.ticks.y = element_blank()) + labs(title = NULL)
  function_plot(filename_prefix = str_glue("{output_file}/UMAP_study_ID"), width = 150, height = 150)
  
  ## cell type proportion
  proportion_data <- scRNA_annotation@meta.data %>%
    dplyr::group_by(tumor_code, cell_type) %>%
    dplyr::summarise(count = n()) %>%
    ungroup()
  
  ## dotplot per tissue
  current_summary <- scRNA_annotation@meta.data %>%
    group_by(tumor_code, new_cell_type) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::mutate(clust_total = sum(count)) %>%
    dplyr::mutate(clust_prop = count / clust_total * 100)
  current_summary$origin_cell_type <- celltype_file[file]
  metadata_summary <- rbind(metadata_summary,current_summary)
  
  ## OR find and plot
  tumor_status_OR <- calTissueDist(scRNA_annotation@meta.data,byPatient = F,colname.cluster = "new_cell_type",colname.patient = "sample_ID",
                                   colname.tissue = "tumor_status",method = "chisq",min.rowSum = 0)
  tumor_status_OR <- tumor_status_OR[,c("NAT","Tumor","Metastasis","Healthy")[c("NAT","Tumor","Metastasis","Healthy") %in% colnames(tumor_status_OR)]]
  pdf(str_glue("{output_file}/cell_type_OR.pdf"), width = ncol(tumor_status_OR)*1.5+2, height = nrow(tumor_status_OR)*0.4+1)
  pheatmap(tumor_status_OR,
           cluster_rows = TRUE,
           color = c("#FEE8C8", "#FDBB84", "#FC8D59", "#EF6548"),
           breaks = c(0, 1, 1.5, 3, max(tumor_status_OR)),
           cluster_cols = FALSE,
           angle_col = 90,
           fontsize = 18, border_color = "white",
           display_numbers = matrix(ifelse(tumor_status_OR > 3, "+++", ifelse(tumor_status_OR > 1.5, "++", 
                                                                              ifelse(tumor_status_OR > 1, "+", "+/-"))), nrow(tumor_status_OR)),
           number_color = "black"
  )
  dev.off()
  
  
  ## find DEGs
  # new_DEGs <- FindAllMarkers(scRNA_annotation,only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
  
  # ## DEGs marker plot
  # new_DEGs_filter <- new_DEGs %>% 
  #   filter(pct.1 >= 0.3 & avg_log2FC >= 1)
  # pdf(sprintf("%s/Figure1F.pdf", figdir), width = 5, height = 8)
  # AverageHeatmap(object = scRNA_annotation, markerGene = new_DEGs_filter,
  #                gene.order = new_DEGs_filter,htRange = c(-3,0,3),
  #                annoCol = T,myanCol = heatmap_color,cluster.order = unique(new_DEGs_filter$cell_type))
  # dev.off()
  
  
  ## meta collect
  all_meta <- rbind(all_meta,scRNA_annotation@meta.data)
  
  ## data save
  # write_csv(new_DEGs,str_glue("{output_file}/cell_type_DEGs.csv"))
  # write_csv(new_DEGs_filter,str_glue("{output_file}/cell_type_DEGs_filter.csv"))
  write_csv(scRNA_annotation@meta.data,str_glue("{output_file}/cell_type_meta.csv"))
  
  ## space free
  remove(scRNA_annotation)
  gc()

   
}


# tumor data process ------------------------------------------------------
## data read
scRNA_tumor_cell <- read_h5(file = str_glue("{tumor_cell_dir}/scRNA_tumor_cell_in_batch.h5"),
                            assay.name = 'RNA', 
                            target.object = 'seurat')

## ROGUE analysis
expr <- GetAssayData(scRNA_tumor_cell,assay = "RNA", layer = "counts") %>% as.matrix()
rogue.res <- rogue(expr, labels = scRNA_tumor_cell@meta.data$leiden,sample = scRNA_tumor_cell@meta.data$sample_ID,
                   platform = "UMI",span = 0.8, min.genes = 100)
rougue_sample <- rogue.res %>% 
  tidyr::gather(key = clusters, value = ROGUE) %>% 
  filter(!is.na(ROGUE))
sample_num <- as.data.frame(table(rougue_sample$clusters))
colnames(sample_num) <- c("sample","num")
ggplot(rougue_sample, aes(clusters, ROGUE, color = clusters)) +
  geom_boxplot(aes(fill = clusters),outlier.shape = NA,alpha = 0.2) +
  geom_jitter_rast(shape = 16, position = position_jitter(0.2),alpha = 0.5,size = 3) +
  # scale_color_manual(values = new_status_color) +
  # scale_fill_manual(values = new_status_color) +
  nrc_theme + NoLegend() + xlab("")


## umap plot
CellDimPlot(scRNA_tumor_cell, reduction = "umap",group.by = "leiden",raster=FALSE,pt.size = 0.05,label.size = 5,label = T,label_insitu = T,label_repel = T) + 
  cluster_theme + NoLegend() + NoLegend() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),axis.ticks.y = element_blank()) + labs(title = NULL)
function_plot(filename_prefix = str_glue("{output_dir}/Epithelia_tumor/UMAP_leiden_cluster"), width = 200, height = 200)

DimPlot(scRNA_tumor_cell, reduction = "umap",group.by = "tumor_status",raster=FALSE,cols = tumor_status_color) +
  cluster_theme + NoLegend() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),axis.ticks.y = element_blank()) + labs(title = NULL)
function_plot(filename_prefix = str_glue("{output_dir}/Epithelia_tumor/UMAP_tumor_status"), width = 200, height = 200)
DimPlot(scRNA_tumor_cell, reduction = "umap",group.by = "tissue_type",raster=FALSE,cols = tissue_type_color) +
  cluster_theme + NoLegend() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),axis.ticks.y = element_blank()) + labs(title = NULL)
function_plot(filename_prefix = str_glue("{output_dir}/Epithelia_tumor/UMAP_tissue_type"), width = 200, height = 200)
DimPlot(scRNA_tumor_cell, reduction = "umap",group.by = "tumor_code",raster=FALSE,cols = tumor_type_color) +
  cluster_theme + NoLegend() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),axis.ticks.y = element_blank()) + labs(title = NULL)
function_plot(filename_prefix = str_glue("{output_dir}/Epithelia_tumor/UMAP_tumor_code"), width = 200, height = 200)
DimPlot(scRNA_tumor_cell, reduction = "umap",group.by = "study_ID",raster=FALSE,cols = study_database_color) +
  cluster_theme + NoLegend() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),axis.ticks.y = element_blank()) + labs(title = NULL)
function_plot(filename_prefix = str_glue("{output_dir}/Epithelia_tumor/UMAP_study_ID"), width = 200, height = 200)

## leiden stat plot
scRNA_tumor_cell_meta <- scRNA_tumor_cell@meta.data
scRNA_tumor_cell_meta_stat <- scRNA_tumor_cell_meta %>% 
  dplyr::group_by(leiden,tumor_code) %>% 
  dplyr::summarise(cell_count = n(),.groups = 'drop') %>% 
  dplyr::group_by(leiden) %>% 
  dplyr::mutate(total_in_leiden = sum(cell_count),
         percentage = round(cell_count/total_in_leiden * 100,2)) %>% 
  arrange(leiden,desc(cell_count))

## stat plot
x_direct_data <- scRNA_tumor_cell_meta_stat %>% 
  filter(percentage >= 0.1)
x_direct_data <- as.data.frame(table(x_direct_data$leiden))
colnames(x_direct_data) <- c("leiden","num")
x_direct_data <- arrange(x_direct_data,desc(num))
scRNA_tumor_cell_meta_stat$leiden <- factor(scRNA_tumor_cell_meta_stat$leiden,levels = x_direct_data$leiden)
ggplot(scRNA_tumor_cell_meta_stat, aes(x = leiden, y = percentage, fill = tumor_code)) + 
  geom_bar(stat = "identity",position = "fill") +
  theme_bw() +
  labs(x = "", y = "Proportion of Cell Type(%)") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = tumor_type_color) +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/Epithelia_tumor/scRNA_tumor_cell_leiden_stat"), width = 600, height = 150)

## data save
write_csv(scRNA_tumor_cell_meta_stat,str_glue("{output_dir}/Epithelia_tumor/scRNA_tumor_cell_meta_stat.csv"))

# all meta stat -----------------------------------------------------------
## new cell type per tumor code
tumor_count <- metadata_summary %>%
  filter(clust_prop > 1) %>%
  dplyr::group_by(new_cell_type) %>%
  dplyr::summarise(tumor_count = n_distinct(tumor_code))
metadata_summary_sorted <- metadata_summary %>%
  left_join(tumor_count, by = "new_cell_type") %>%
  arrange(origin_cell_type, -tumor_count, new_cell_type) %>% 
  mutate(label = ifelse(clust_prop >= 10, ">10%",
                        ifelse(clust_prop >= 5, "5%-10%",
                               ifelse(clust_prop >= 1, "1%-5%","<1%"))))
metadata_summary_sorted <- metadata_summary_sorted[metadata_summary_sorted$origin_cell_type != "Epithelia",]
metadata_summary_sorted$new_cell_type <- factor(metadata_summary_sorted$new_cell_type)

mycols <- brewer.pal(9, "GnBu")[c(1, 3, 5, 7)]
names(mycols) <- c("<1%", "1%-5%", "5%-10%", ">10%")
ggplot(data = metadata_summary_sorted, mapping = aes_string(x = "tumor_code", y = "new_cell_type")) +
  geom_point_rast(mapping = aes_string(size = "clust_prop", color = "label")) +
  scale_size_continuous(breaks = c(0,1,5,10,15,20)) +
  scale_x_discrete(limits = c("ANS","BRCA","CESC","CRC","ESCC","GC","HNSC","KIRC","LSCC","NSCLC","LUAD","LUSC","MA","NPC","PAAD","PNET",
                              "SARC","TGCT","THCA","ACC","Liver_healthy","Lymph_node_healthy","Brain_healthy")) +
  theme_bw() +
  theme(
    strip.text.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(colour = "grey", linetype = "dashed", size = 0.2)
  ) +
  scale_color_manual(values = mycols)
function_plot(filename_prefix = str_glue("{output_dir}/cell_type_proportion_dotplot"), width = 250, height = 600)

## all OR
all_tumor_status_OR <- calTissueDist(all_meta,byPatient = F,colname.cluster = "new_cell_type",colname.patient = "sample_ID",
                                 colname.tissue = "tumor_status",method = "chisq",min.rowSum = 0)
all_tumor_status_OR <- all_tumor_status_OR[rev(levels(metadata_summary_sorted$new_cell_type)),c("NAT","Tumor","Metastasis","Healthy")]
pdf(str_glue("{output_dir}/cell_type_OR_pheatmap.pdf"), width = ncol(all_tumor_status_OR)*1.5+3.5, height = nrow(all_tumor_status_OR)*0.25+1)
pheatmap(all_tumor_status_OR,
         cluster_rows = FALSE,
         color = c("#FEE8C8", "#FDBB84", "#FC8D59", "#EF6548"),
         breaks = c(0, 1, 1.5, 3, max(all_tumor_status_OR)),
         cluster_cols = FALSE,
         angle_col = 45,
         fontsize = 18, border_color = "white",
         display_numbers = matrix(ifelse(all_tumor_status_OR > 3, "+++", ifelse(all_tumor_status_OR > 1.5, "++", ifelse(all_tumor_status_OR > 1, "+", "+/-"))), nrow(all_tumor_status_OR)),
         number_color = "black"
)
dev.off()



# metastasis tissue stat --------------------------------------------------
metastasis_meta <- all_meta %>% 
  filter(tumor_status == "Metastasis")
metastasis_OR <- calTissueDist(metastasis_meta,byPatient = F,colname.cluster = "new_cell_type",colname.patient = "sample_ID",
                                     colname.tissue = "tissue_type",method = "chisq",min.rowSum = 0)

pdf(sprintf("%s/Figure1F.pdf", figdir), width = nrow(metastasis_OR)*1+2, height = ncol(metastasis_OR)*1+2)
pheatmap(metastasis_OR,
         cluster_rows = FALSE,
         color = c("#FEE8C8", "#FDBB84", "#FC8D59", "#EF6548"),
         breaks = c(0, 1, 1.5, 3, max(metastasis_OR,na.rm = T)),
         cluster_cols = FALSE,
         angle_col = 45,
         fontsize = 18, border_color = "white",
         display_numbers = matrix(ifelse(metastasis_OR > 3, "+++", ifelse(metastasis_OR > 1.5, "++", ifelse(metastasis_OR > 1, "+", "+/-"))), nrow(metastasis_OR)),
         number_color = "black"
)
dev.off()


# all cell type num stat --------------------------------------------------
## proportion in tumor code and tissue type
for(celltype in as.character(unique(all_meta$new_cell_type))){
  # celltype <- "EC_12_ADAMTSL1"
  ## data filter
  current_meta <- all_meta %>% 
    filter(new_cell_type == celltype)
  
  ## data stat
  tumor_code_num <- as.data.frame(table(current_meta$tumor_code))
  colnames(tumor_code_num) <- c("tumor_code","num")
  tumor_code_num <- tumor_code_num %>% 
    filter(num >= 50)
  tissue_type_num <- as.data.frame(table(current_meta$tissue_type))
  colnames(tissue_type_num) <- c("tissue_type","num")
  tissue_type_num <- tissue_type_num %>% 
    filter(num >= 50)
  
  ## data collect
  current_stat <- data.frame(
    cell_type = celltype,
    tumor_code_stat = nrow(tumor_code_num),
    tissue_type_stat = nrow(tissue_type_num)
  )
  all_cell_num_stat <- rbind(all_cell_num_stat,current_stat)
  
}

all_cell_num_stat <- all_cell_num_stat %>% 
  mutate(status = ifelse(tumor_code_stat >= 10 | tissue_type_stat >= 12, "high.in.50%",
                         ifelse(tumor_code_stat >= 20/3 | tissue_type_stat >= 24/3, "high.in.33%","low.percentage")))%>%
  mutate(label = ifelse(status %in% "high.in.50%",cell_type,""))
all_cell_num_stat <- all_cell_num_stat[!grepl("^Ep",all_cell_num_stat$cell_type),]
ggplot(all_cell_num_stat,aes(tumor_code_stat, tissue_type_stat)) +
  geom_point(aes(size = 7,color = status)) + 
  geom_vline(xintercept = c(20/3,10), lty = 4, col = '#a4b0be', lwd = 1) +
  geom_hline(yintercept = c(8,12), lty = 4, col = '#a4b0be', lwd = 1) + 
  geom_text_repel(aes(label = label,color = "red"),box.padding = 1.5, min.segment.length = 0,color = 'black',nudge_x = 0.2, nudge_y = 0.2,
                  point.padding = 1, max.overlaps = getOption("ggrepel.max.overlaps",default = 200),max.time = 10, max.iter = 80000) +
  xlim(c(0,25)) + ylim(c(0,25)) +
  nrc_theme + NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/cell_type_percentage_stat"), width = 320, height = 320)
ggplot(all_cell_num_stat, aes(reorder(cell_type, - tumor_code_stat), tumor_code_stat)) + 
  geom_bar(stat = "identity", fill = "#2285FE") +
  theme_bw() +
  labs(x = "", y = "Number of tumor") +
  scale_y_continuous(expand = c(0.01,0)) +
  nrc_theme + NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/cell_type_num_tumor_code"), width = 700, height = 250)
ggplot(all_cell_num_stat, aes(reorder(cell_type, - tissue_type_stat), tissue_type_stat)) + 
  geom_bar(stat = "identity", fill = "#2285FE") +
  theme_bw() +
  labs(x = "", y = "Number of tissue") +
  scale_y_continuous(expand = c(0.01,0)) +
  nrc_theme + NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/cell_type_num_tissue_type"), width = 700, height = 250)

## cell num stat
all_cell_type_num <- as.data.frame(table(all_meta$new_cell_type))
colnames(all_cell_type_num) <- c("cell_type","num")
all_cell_type_num <- all_cell_type_num[!grepl("^Ep",all_cell_type_num$cell_type),]
ggplot(all_cell_type_num, aes(reorder(cell_type, - num), num)) + 
  geom_bar(stat = "identity", fill = "#2285FE") +
  theme_bw() +
  labs(x = "", y = "Number of cell") +
  scale_y_continuous(expand = c(0.01,0)) +
  # scale_fill_manual(values = study_database_color) +
  nrc_theme + NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/cell_type_num_stat"), width = 700, height = 250)

# all_meta <- read_csv(str_glue("{output_dir}/cell_type_meta_data_all.csv"))
# metadata_summary <- read_csv(str_glue("{output_dir}/cell_type_metadata_summary.csv"))

# data save ---------------------------------------------------------------
if(T){
  write_csv(metadata_summary,str_glue("{output_dir}/cell_type_metadata_summary.csv"))
  write_csv(all_meta,str_glue("{output_dir}/cell_type_meta_data_all.csv"))
  write_csv(all_cell_num_stat,str_glue("{output_dir}/cell_type_tumor_tissue_num.csv"))
  write_csv(all_cell_type_num,str_glue("{output_dir}/cell_type_all_num.csv"))
}

