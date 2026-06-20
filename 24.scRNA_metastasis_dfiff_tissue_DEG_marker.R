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
library("UpSetR")
# library("scPDtools")
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
nohub_gene_dir <- " "
input_dir <- " "
hall_dir <- " "
output_dir <- " "
tumor_code_list <- c("NPC","ACC","KIRC","HNSC","CRC","GC","PAAD","BRCA","CESC","ESCC","PNET",
                     "THCA","NSCLC","LSCC")


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

DEGs_all <- data.frame(
  gene_name = as.character(),
  avg_logFC = as.numeric(),
  pct.1 = as.numeric(),
  pct.2 = as.numeric(),
  p_val_adj = as.numeric(),
  tumor_code = as.character(),
  metastasis_tissue = as.character()
)

metastasis_up_gene_count <- data.frame(
  gene_name = as.character(),
  num = as.numeric(),
  metastasis_tissue = as.character()
)
primary_up_gene_count <- data.frame(
  gene_name = as.character(),
  num = as.numeric(),
  metastasis_tissue = as.character()
)

# data read ---------------------------------------------------------------
nohub_gene <- read.xlsx(nohub_gene_dir,sheetName = "gene_name")
hall_file <- clusterProfiler::read.gmt(str_glue("{hall_dir}/h.all.v2023.1.Hs.symbols.gmt"))
scRNA_tumor_cell <- read_h5(file = str_glue("{input_dir}/scRNA_bbknn.h5"),
                            assay.name = 'RNA', 
                            target.object = 'seurat')
scRNA_tumor_cell@meta.data$tumor_status_code <- str_c(scRNA_tumor_cell@meta.data$metastasis_tissue,
                                                      scRNA_tumor_cell@meta.data$tumor_status,sep = "_")



# HLA gene exp ------------------------------------------------------------
hub_marker <- c(
  "HLA-DRA","HLA-DPA1","HLA-DRB1","HLA-DPB1","HLA-DQA1","HLA-B","HLA-C","CD74","B2M"
)
y_direct <- c("None_Tumor","Brain_Metastasis","Liver_Metastasis","Lymph_node_Metastasis","Peritoneal_Metastasis",
              "Lung_Metastasis","Ovary_Metastasis","Bone_marrow_Metastasis","Vaginal_Metastasis")
my_palette <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(n = 20)
DotPlot(scRNA_tumor_cell, features = hub_marker,group.by = "tumor_status_code") + nrc_theme + ylab("") + xlab("") +
  scale_y_discrete(limits=y_direct) +
  scale_size_continuous(range = c(0.1, 6),breaks = c(0,20,40,60,80,100)) +
  scale_color_gradientn(colors = my_palette)
function_plot(filename_prefix = str_glue("{output_dir}/HLA_gene_marker"), width = 180, height = 100)



# data run ----------------------------------------------------------------
for(tumor in tumor_code_list){
  # tumor <- "KIRC"
  print(str_glue("############################## {tumor} process ##################################"))
  
  if(tumor == "NSCLC"){
    ## data filter
    current_scRNA <- subset(scRNA_tumor_cell,tumor_code %in% c("NSCLC","LUAD","LUSC"))
  }else{
    ## data filter
    current_scRNA <- subset(scRNA_tumor_cell,tumor_code %in% tumor)
  }

  
  for(meta_tissue in c("Other","Lymph_node","Brain","Liver")){
    # meta_tissue <- "Other"
    
    if(meta_tissue %in% current_scRNA@meta.data$metastasis_tissue){
      print(str_glue("           --------------- {meta_tissue} process"))
      ## DEGs run
      DEGs_current <- FindMarkers(current_scRNA,group.by = "tumor_status_code",ident.1 = str_glue("{meta_tissue}_Metastasis"),ident.2 = "None_Tumor")
      
      ## DEGs result correct
      DEGs_current <- DEGs_current %>% 
          rownames_to_column(var = 'gene_name')
      DEGs_current$tumor_code <- tumor
      DEGs_current$metastasis_tissue <- meta_tissue
      
      ## data combine
      DEGs_all <- rbind(DEGs_all,DEGs_current)
      
    }else if(meta_tissue == "Other"){
      print(str_glue("           --------------- {meta_tissue} process"))
      
      hub_list <- unique(current_scRNA@meta.data$metastasis_tissue)[!unique(current_scRNA@meta.data$metastasis_tissue) %in% c("Lymph_node","Brain","Liver")]
      current_scRNA_hub <- subset(current_scRNA,metastasis_tissue %in% hub_list)
      
      if(!"Metastasis" %in% current_scRNA_hub@meta.data$tumor_status | !"Tumor" %in% current_scRNA_hub@meta.data$tumor_status){
        next
      }
      
      ## DEGs run
      
      DEGs_current <- FindMarkers(current_scRNA_hub,group.by = "tumor_status",ident.1 = "Metastasis",ident.2 = "Tumor")
      
      ## DEGs result correct
      DEGs_current <- DEGs_current %>% 
        rownames_to_column(var = 'gene_name')
      DEGs_current$tumor_code <- tumor
      DEGs_current$metastasis_tissue <- meta_tissue
      
      ## data combine
      DEGs_all <- rbind(DEGs_all,DEGs_current)
      
    }
    
  }
  
  remove(current_scRNA)
  gc()

}

# a <- readRDS(" ")
# DEGs filter -------------------------------------------------------------
DEGs_all <- read_csv(" ")
DEGs_filter <- DEGs_all %>% 
  mutate(DEGs_status = ifelse(pct.1 >= 0.3 & avg_log2FC >= log2(1.1),"up.in.metastasis",
                              ifelse(pct.2 >= 0.3 & avg_log2FC <= log2(1/1.1),"up.in.tumor","not"))) %>% 
  filter(DEGs_status != "not")
DEGs_filter <- DEGs_filter[!DEGs_filter$gene_name %in% nohub_gene$gene_name,]



# DEG num plot ------------------------------------------------------------
DEGs_upset <- DEGs_all %>% 
  mutate(DEGs_status = ifelse(pct.1 >= 0.2 & avg_log2FC > log2(1),"up.in.metastasis",
                              ifelse(pct.2 >= 0.2 & avg_log2FC < log2(1/1),"up.in.tumor","not"))) %>% 
  filter(DEGs_status != "not")
DEGs_upset <- DEGs_upset[!DEGs_upset$gene_name %in% nohub_gene$gene_name,]
DEGs_upset$group <- str_c(DEGs_upset$metastasis_tissue,DEGs_upset$tumor_code,sep = "_")
DEGs_upset_metastasis <- DEGs_upset[DEGs_upset$DEGs_status == "up.in.metastasis",]
DEGs_upset_primary <- DEGs_upset[DEGs_upset$DEGs_status == "up.in.tumor",]

## num stat
DEGs_upset_metastasis_num <- as.data.frame(table(DEGs_upset_metastasis$gene_name))
colnames(DEGs_upset_metastasis_num) <- c("gene_name","num")
DEGs_upset_metastasis_num <- arrange(DEGs_upset_metastasis_num,desc(num))
DEGs_upset_primary_num <- as.data.frame(table(DEGs_upset_primary$gene_name))
colnames(DEGs_upset_primary_num) <- c("gene_name","num")
DEGs_upset_primary_num <- arrange(DEGs_upset_primary_num,desc(num))

## data plot
ggplot(DEGs_upset_metastasis_num[1:50,], aes(reorder(gene_name, - num), num)) + 
  geom_bar(stat = "identity", fill = "#2285FE") +
  theme_bw() +
  labs(x = "", y = "Number of tumor code") +
  scale_y_continuous(expand = c(0.01,0)) +
  nrc_theme + NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/all_metastasis_stat"), width = 300, height = 150)
ggplot(DEGs_upset_primary_num[1:50,], aes(reorder(gene_name, - num), num)) + 
  geom_bar(stat = "identity", fill = "#2285FE") +
  theme_bw() +
  labs(x = "", y = "Number of tumor code") +
  scale_y_continuous(expand = c(0.01,0)) +
  nrc_theme + NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/all_priamry_stat"), width = 300, height = 150)

# pancancer marker find ---------------------------------------------------
## all DEG plot
DEGs_all_volcano <- DEGs_all
colnames(DEGs_all_volcano)[1] <- "gene"
DEGs_all_volcano$cluster <- str_c(DEGs_all_volcano$metastasis_tissue,DEGs_all_volcano$tumor_code,sep = '_')
jjVolcano(diffData = DEGs_all_volcano,topGeneN = 0,base_size = 6,tile.col = as.character(tumor_type_color)) + NoLegend()
function_plot(filename_prefix = str_glue("{output_dir}/metastasis_all_tumor_marker"), width = 700, height = 200)



## marker find and plot
high_per_DEG_all <- data.frame(
  gene_name = as.character()
)
for(meta_tissue in c("Lymph_node","Brain","Liver","Other")){
  # meta_tissue <- "Other"
  
  DEGs_filter_meta <- DEGs_filter[DEGs_filter$metastasis_tissue == meta_tissue,]
  DEGs_filter_meta <- DEGs_filter_meta[!str_detect(DEGs_filter_meta$gene_name,"^ENSG"),]
  
  ## marker num
  DEGs_primary_up <- DEGs_filter_meta[DEGs_filter_meta$DEGs_status == "up.in.tumor",]
  DEGs_primary_up_count <- as.data.frame(table(DEGs_primary_up$gene_name))
  colnames(DEGs_primary_up_count) <- c("gene_name","num")
  DEGs_primary_up_count <- arrange(DEGs_primary_up_count,desc(num))
  DEGs_primary_up_count$metastasis_tissue <- meta_tissue
  DEGs_metastasis_up <- DEGs_filter_meta[DEGs_filter_meta$DEGs_status == "up.in.metastasis",]
  DEGs_metastasis_up_count <- as.data.frame(table(DEGs_metastasis_up$gene_name))
  colnames(DEGs_metastasis_up_count) <- c("gene_name","num")
  DEGs_metastasis_up_count <- arrange(DEGs_metastasis_up_count,desc(num))
  DEGs_metastasis_up_count$metastasis_tissue <- meta_tissue
  
  
  ## upset plot
  # upset_plot_data_metastasis <- as.data.frame(table(DEGs_metastasis_up$gene_name,DEGs_metastasis_up$tumor_code) > 0)
  # upset_plot_data_metastasis <- upset_plot_data_metastasis * 1
  # pdf(str_glue("{output_dir}/upset_metastasis_{meta_tissue}.pdf"),width = 12,height = 8)
  # upset(upset_plot_data_metastasis,nsets = 23,nintersects = 100,matrix.color = "black",order.by = "freq")
  # dev.off()
  # upset_plot_data_primary <- as.data.frame(table(DEGs_primary_up$gene_name,DEGs_primary_up$tumor_code) > 0)
  # upset_plot_data_primary <- upset_plot_data_primary * 1
  # pdf(str_glue("{output_dir}/upset_primary_{meta_tissue}.pdf"),width = 12,height = 8)
  # upset(upset_plot_data_primary,nsets = 23,nintersects = 100,matrix.color = "black",order.by = "freq")
  # dev.off()
  
  ## marker filter
  DEGs_plot <- DEGs_filter_meta %>% 
    filter(gene_name %in% c(DEGs_metastasis_up_count$gene_name[1:20],DEGs_primary_up_count$gene_name[1:20]))
  DEGs_plot$gene_name <- factor(DEGs_plot$gene_name,levels = c(DEGs_metastasis_up_count$gene_name[1:20],DEGs_primary_up_count$gene_name[1:20]))
  high_per_DEG_all <- rbind(high_per_DEG_all,DEGs_plot)
  
  ## marker plot
  ggplot(DEGs_plot, aes(x = gene_name, y = avg_log2FC, fill = tumor_code)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = tumor_type_color) +
    labs(x = "", y = "log2FC") + 
    geom_hline(yintercept = log2(1.2),linetype = "dashed",color = "black") +
    geom_hline(yintercept = log2(1/1.2),linetype = "dashed",color = "black") +
    nrc_theme
  function_plot(filename_prefix = str_glue("{output_dir}/metastasis_{meta_tissue}_marker"), width = 320, height = 120)
  
  ## enrich gene filter
  DEGs_primary_high_gene <- DEGs_primary_up_count %>% 
    filter(num >= 4)
  DEGs_metastasis_high_gene <- DEGs_metastasis_up_count %>% 
    filter(num >= 4)
  
  ## enrichment
  DEGs_primary_enrich <- enricher(gene = DEGs_primary_high_gene$gene_name, TERM2GENE = hall_file,
                                  pAdjustMethod = "none", pvalueCutoff = 0.05, qvalueCutoff = 1, minGSSize = 3)
  DEGs_primary_enrich <- clusterProfiler::filter(DEGs_primary_enrich, Count >= 3)
  DEGs_primary_enrich <- as.data.frame(head(DEGs_primary_enrich, nrow(DEGs_primary_enrich)))
  DEGs_primary_enrich <- as.data.frame(head(DEGs_primary_enrich, 10))
  DEGs_primary_enrich$log10pvalue <- -log10(DEGs_primary_enrich$pvalue)
  DEGs_primary_enrich$log10pvalue <- -DEGs_primary_enrich$log10pvalue
  
  DEGs_metastasis_enrich <- enricher(gene = DEGs_metastasis_high_gene$gene_name, TERM2GENE = hall_file,
                                     pAdjustMethod = "none", pvalueCutoff = 0.05, qvalueCutoff = 1, minGSSize = 3)
  DEGs_metastasis_enrich <- clusterProfiler::filter(DEGs_metastasis_enrich, Count >= 3)
  DEGs_metastasis_enrich <- as.data.frame(head(DEGs_metastasis_enrich, nrow(DEGs_metastasis_enrich)))
  DEGs_metastasis_enrich <- as.data.frame(head(DEGs_metastasis_enrich, 10))
  DEGs_metastasis_enrich$log10pvalue <- -log10(DEGs_metastasis_enrich$pvalue)
  
  ## enrichment plot
  combine_enrich <- rbind(DEGs_metastasis_enrich,DEGs_primary_enrich)
  combine_enrich <- arrange(combine_enrich,desc(log10pvalue))
  combine_enrich$Description <- factor(combine_enrich$Description,levels = unique(combine_enrich$Description))
  combine_enrich %>%
    ggplot(aes(x=Description, y=log10pvalue,fill = Count)) + 
    geom_bar(stat = "identity", width = 0.8, colour = 'black') + 
    geom_text(aes(label = Description), y = 0.2,
              fontface = 'bold', size = 4, hjust = 0) + 
    scale_fill_distiller(palette = "YlOrRd", direction = 1) +
    theme(panel.border = element_rect(color = "black", size = 1.5, linetype = "solid", inherit.blank = TRUE, fill='transparent')) +
    theme(panel.background = element_rect(fill = "transparent")) + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
    theme(axis.title = element_text(size = 15, face = 'bold')) +
    theme(axis.text = element_text(color="black", size = 12, face = 'bold')) +
    theme(axis.ticks.y = element_line(size = 1), axis.ticks.length.y = unit("0", "mm")) +
    theme(axis.ticks.x = element_line(size = 1), axis.ticks.length.x = unit("2", "mm")) +
    xlab("Hallmark") +
    scale_y_continuous(expand = c(0.002, 0.002, 0.02, 0.002)) + 
    scale_x_discrete(label=function(x) " ") +
    coord_flip() 
  function_plot(filename_prefix = str_glue("{output_dir}/{meta_tissue}_metastasis_enrich"), width = 200, height = 220)
  
  ## data save
  write_csv(DEGs_primary_up_count,str_glue("{output_dir}/{meta_tissue}_DEGs_tumor_up_count.csv"))
  write_csv(DEGs_metastasis_up_count,str_glue("{output_dir}/{meta_tissue}_DEGs_metastasis_up_count.csv"))
  
  ## high gene rbind
  metastasis_up_gene_count <- rbind(metastasis_up_gene_count,DEGs_metastasis_up_count[DEGs_metastasis_up_count$num >= 5,])
  primary_up_gene_count <- rbind(primary_up_gene_count,DEGs_primary_up_count[DEGs_primary_up_count$num >= 5,])
  
}

## high per DEG enrichment
high_per_DEG_metastasis <- high_per_DEG_all[high_per_DEG_all$DEGs_status=="up.in.metastasis",]
high_per_DEG_metastasis_enrich <- enricher(gene = as.character(unique(high_per_DEG_metastasis$gene_name)), TERM2GENE = hall_file,
                                pAdjustMethod = "none", pvalueCutoff = 0.05, qvalueCutoff = 1, minGSSize = 3)
high_per_DEG_metastasis_enrich <- clusterProfiler::filter(high_per_DEG_metastasis_enrich, Count >= 3)
high_per_DEG_metastasis_enrich <- as.data.frame(head(high_per_DEG_metastasis_enrich, nrow(high_per_DEG_metastasis_enrich)))
high_per_DEG_all_metastasis <- as.data.frame(unique(high_per_DEG_all[high_per_DEG_all$DEGs_status == "up.in.metastasis",]$gene_name))
write_csv(high_per_DEG_all_metastasis,str_glue("{output_dir}/high_per_DEG_all_metastasis.csv"))


# diff metastasis tissue marker find ---------------------------------------------------
scRNA_metastasis_cell <- subset(scRNA_tumor_cell,tumor_status == "Metastasis")
scRNA_metastasis_cell@meta.data$tumor_tissue_type <- str_c(scRNA_metastasis_cell@meta.data$tumor_code,
                                                          scRNA_metastasis_cell@meta.data$tissue_type,
                                                          scRNA_metastasis_cell@meta.data$tumor_status,sep = "_")
Idents(scRNA_metastasis_cell) <- scRNA_metastasis_cell$tumor_tissue_type

## diff metastasis tissue DEG
# diff_metastasis_tissue_DEG <- read_csv(str_glue("{output_dir}/DEGs_diff_metastasis_tissue_all.csv"))
diff_metastasis_tissue_DEG <- FindAllMarkers(scRNA_metastasis_cell,group.by = "tumor_tissue_type")
diff_metastasis_tissue_DEG_filter <- diff_metastasis_tissue_DEG %>% 
  rownames_to_column(var = "gene_name") %>% 
  filter(pct.1 >= 0.3 & avg_log2FC >= log2(1.2)) %>% 
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>% 
  # slice_head(n = 200) %>%
  ungroup()

## diff metastasis tissue high DEGs
brain_DEG <- diff_metastasis_tissue_DEG_filter[grepl("Brain",diff_metastasis_tissue_DEG_filter$cluster),]
brain_DEG_num <- as.data.frame(table(brain_DEG$gene))
colnames(brain_DEG_num) <- c("gene_name","num")
brain_DEG_num <- arrange(brain_DEG_num,desc(num))
liver_DEG <- diff_metastasis_tissue_DEG_filter[grepl("Liver",diff_metastasis_tissue_DEG_filter$cluster),]
liver_DEG_num <- as.data.frame(table(liver_DEG$gene))
colnames(liver_DEG_num) <- c("gene_name","num")
liver_DEG_num <- arrange(liver_DEG_num,desc(num))
lymph_node_DEG <- diff_metastasis_tissue_DEG_filter[grepl("Lymph_node_",diff_metastasis_tissue_DEG_filter$cluster),]
lymph_node_DEG_num <- as.data.frame(table(lymph_node_DEG$gene))
colnames(lymph_node_DEG_num) <- c("gene_name","num")
lymph_node_DEG_num <- arrange(lymph_node_DEG_num,desc(num))
lung_DEG <- diff_metastasis_tissue_DEG_filter[grepl("Lung",diff_metastasis_tissue_DEG_filter$cluster),]
lung_DEG_num <- as.data.frame(table(lung_DEG$gene))
colnames(lung_DEG_num) <- c("gene_name","num")
lung_DEG_num <- arrange(lung_DEG_num,desc(num))
peritoneal_DEG <- diff_metastasis_tissue_DEG_filter[grepl("Peritoneal",diff_metastasis_tissue_DEG_filter$cluster),]
peritoneal_DEG_num <- as.data.frame(table(peritoneal_DEG$gene))
colnames(peritoneal_DEG_num) <- c("gene_name","num")
peritoneal_DEG_num <- arrange(peritoneal_DEG_num,desc(num))
bone_marrow_DEG <- diff_metastasis_tissue_DEG_filter[grepl("Bone_marrow",diff_metastasis_tissue_DEG_filter$cluster),]
vaginal_DEG <- diff_metastasis_tissue_DEG_filter[grepl("Vaginal",diff_metastasis_tissue_DEG_filter$cluster),]
ovary_DEG <- diff_metastasis_tissue_DEG_filter[grepl("Ovary",diff_metastasis_tissue_DEG_filter$cluster),]

## result plot
cell_marker <- c(
                 # "NDUFA4L2","FXYD2","FABP7","CRYAB",
                 "NDUFA4L2","SERPINE1","UGT2B7","UBE2D2","KCNIP4",
                 # "LYZ","ODAM","GPX2","TFF2",
                 "SMPX","SLC40A1","MRPL27","LCN2","TMEM92",
                 # "PCK1","ANPEP","ARFGEF2","PIGR",
                 "MED15","SCN2A","ANPEP","ARFGEF2","PHF20",
                 "ANXA2","AVPI1","CASC15","SIK1",
                 # "KRT17","MYO1E","PERP","PLAUR",
                 # "FGD4","CD2AP","CEMIP2","CYP3A5",
                 "LLGL2","CLMN","INPP4B","JAK1","MAGI1",
                 # "ALDOA","ATP5F1B","COX5A","EEF2",
                 "GABARAP","NDUFA7","NDUFS7","EEF2",
                 # "BST2","HLA-DRA","IFI16","RIMS2",
                 "SNHG25","ATP11A","CNIH4","SELENOM",
                 # "PRDX1","BANF1","ERH","HNRNPA1"
                 "PRDX1","ERH","HNRNPA1","BANF1"
                 )



y_direct <- c("KIRC_Bone_marrow_Metastasis","PAAD_Vaginal_Metastasis","GC_Ovary_Metastasis",
              "PAAD_Lung_Metastasis","ACC_Lung_Metastasis","CRC_Peritoneal_Metastasis","GC_Peritoneal_Metastasis","ANS_Peritoneal_Metastasis",
              "CESC_Liver_Metastasis","CRC_Liver_Metastasis","NPC_Liver_Metastasis","PAAD_Liver_Metastasis","BRCA_Liver_Metastasis",
              "PNET_Liver_Metastasis","GC_Liver_Metastasis",
              "BRCA_Brain_Metastasis","CRC_Brain_Metastasis","NSCLC_Brain_Metastasis","ESCC_Brain_Metastasis","LUAD_Brain_Metastasis","THCA_Brain_Metastasis",
              "LSCC_Lymph_node_Metastasis","TGCT_Lymph_node_Metastasis","BRCA_Lymph_node_Metastasis","ESCC_Lymph_node_Metastasis",
              "HNSC_Lymph_node_Metastasis","LUAD_Lymph_node_Metastasis","LUSC_Lymph_node_Metastasis","THCA_Lymph_node_Metastasis"
)
Idents(scRNA_metastasis_cell) <- scRNA_metastasis_cell$tumor_tissue_type
DotPlot(scRNA_metastasis_cell, features = cell_marker,group.by = "tumor_tissue_type") + nrc_theme + ylab("") + xlab("") +
      scale_y_discrete(limits=y_direct) +
      scale_size_continuous(range = c(0.1, 6),breaks = c(0,20,40,60,80,100)) +
      scale_color_gradient2(high = "#D7301F",mid = "#A1D99B",low = "#6BAED6")
function_plot(filename_prefix = str_glue("{output_dir}/dotplot_diff_metastasis_tissue_marker"), width = 355, height = 185)


# diff tissue marker and pathway find -------------------------------------
Idents(scRNA_metastasis_cell) <- scRNA_metastasis_cell$metastasis_tissue
metastasis_DEG <- FindAllMarkers(scRNA_metastasis_cell,only.pos = T)
metastasis_DEG_filter <- metastasis_DEG %>% 
  filter(pct.1 >= 0.5 & avg_log2FC >= log2(1.2))





# data save ---------------------------------------------------------------
if(T){
  write_csv(DEGs_all,str_glue("{output_dir}/DEGs_all.csv"))
  write_csv(DEGs_filter,str_glue("{output_dir}/DEGs_filter.csv"))
  write_csv(diff_metastasis_tissue_DEG,str_glue("{output_dir}/DEGs_diff_metastasis_tissue_all.csv"))
  write_csv(diff_metastasis_tissue_DEG_filter,str_glue("{output_dir}/DEGs_diff_metastasis_tissue_filter.csv"))

  write_csv(metastasis_DEG,str_glue("{output_dir}/metastasis_tissue_DEG.csv"))
  
  write_csv(brain_DEG_num,str_glue("{output_dir}/tissue_marker_brain_DEG_num.csv"))
  write_csv(liver_DEG_num,str_glue("{output_dir}/tissue_marker_liver_DEG_num.csv"))
  write_csv(lymph_node_DEG_num,str_glue("{output_dir}/tissue_marker_lymph_node_DEG_num.csv"))
  write_csv(lung_DEG_num,str_glue("{output_dir}/tissue_marker_lung_DEG_num.csv"))
  write_csv(peritoneal_DEG_num,str_glue("{output_dir}/tissue_marker_peritoneal_DEG_num.csv"))
  write_csv(bone_marrow_DEG,str_glue("{output_dir}/tissue_marker_bone_marrow_DEG.csv"))
  write_csv(vaginal_DEG,str_glue("{output_dir}/tissue_marker_vaginal_DEG.csv"))
  write_csv(ovary_DEG,str_glue("{output_dir}/tissue_marker_ovary_DEG.csv"))
}


