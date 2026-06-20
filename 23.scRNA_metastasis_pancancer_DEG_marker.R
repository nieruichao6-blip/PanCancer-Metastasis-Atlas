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
library("scPDtools")
(.packages())
setwd(' ')
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
input_dir <- " "
hall_dir <- " "
output_dir <- " "
tumor_code_list <- c("ACC","BRCA","CESC","CRC","ESCC","HNSC","KIRC","NPC","PAAD","PNET","GC",
                     "THCA","NSCLC","LSCC")


# data color --------------------------------------------------------------
tumor_type_color <- c("ACC"="#5A7B8F","ANS"="#EE4C97","BRCA"="#3D806F","CESC"="#A08634","CRC"="#F37C95","ESCC"="#608541",
                      "HNSC"="#7D4E57","KIRC"="#BC3C29","LUAD"="#958056","LUSC"="#9FAFA3","NPC"="#6F99AD","PAAD"="#0072B5",
                      "PNET"="#CFC59A","SARC"="#E18727","SKCM"="#FFDC91","GC"="#718DAE","TGCT"="#7876B1","THCA"="#F9AC93",
                      "NSCLC"="#E3BC06","LSCC"="#2285FE","Liver_healthy"="#3E6086","Lymph_node_healthy"="#20854E","Brain_healthy"="#4A7985")
tissue_type_color <- c("Bone_marrow"="#CE1261","Kidney"="#5E4FA2","Lymph_node"="#8CA77B","Brain"="#00441B","Lung"="#ed3d30",
                       "Colon"="#DEDC00","Liver"="#B3DE69","Peritoneal"="#999999","Stomach"="#725e82","Ovary"="#8c4356",
                       "Head_and_neck"="#db5a6b","Pancreas"="#21a675","Vaginal"="#d9b611","Appendix"="#e29c45",
                       "Bone"="#4c8dae","Nasopharynx"="#003371","Salivary_gland"="#a1afc9","Testis"="#88ada6","Cervix"="#c93756",
                       "Skin"="#204B75","Breast"="#588257","Thyroid"="#B6DB7B","Esophagus"="#4285BF","larynx"="#1CB232")

DEGs_all <- data.frame(
  gene_name = as.character(),
  avg_logFC = as.numeric(),
  pct.1 = as.numeric(),
  pct.2 = as.numeric(),
  p_val_adj = as.numeric()
)

# data read ---------------------------------------------------------------
hall_file <- clusterProfiler::read.gmt(str_glue("{hall_dir}/h.all.v2023.1.Hs.symbols.gmt"))
scRNA_tumor_cell <- read_h5(file = str_glue("{input_dir}/scRNA_bbknn.h5"),
                            assay.name = 'RNA', 
                            target.object = 'seurat')
scRNA_tumor_cell@meta.data$tumor_status_code <- str_c(scRNA_tumor_cell@meta.data$tumor_code,
                                                      scRNA_tumor_cell@meta.data$tumor_status,sep = "_")


# DEGs find ---------------------------------------------------------------
for(tumor in tumor_code_list){
  # tumor <- "LSCC"
  
  print(str_glue("############################## {tumor} process ##################################"))
  
  if(tumor == "NSCLC"){
    ## data filter
    current_scRNA <- subset(scRNA_tumor_cell,tumor_code %in% c("NSCLC","LUAD","LUSC"))
    ## DEGs run
    DEGs_current <- FindMarkers(scRNA_tumor_cell,group.by = "tumor_status",ident.1 = "Tumor",ident.2 = "Metastasis")
  }else{
    ## data filter
    current_scRNA <- subset(scRNA_tumor_cell,tumor_code %in% tumor)
    ## DEGs run
    DEGs_current <- FindMarkers(scRNA_tumor_cell,group.by = "tumor_status_code",ident.1 = str_glue("{tumor}_Tumor"),ident.2 = str_glue("{tumor}_Metastasis"))
  }

  
  ## DEGs result correct
  DEGs_current <- DEGs_current %>% 
    rownames_to_column(var = 'gene_name')
  DEGs_current$tumor_code <- tumor
  
  ## data combine
  DEGs_all <- rbind(DEGs_all,DEGs_current)
  
}


# DEGs filter -------------------------------------------------------------
DEGs_filter <- DEGs_all %>% 
  mutate(DEGs_status = ifelse(pct.1 >= 0.2 & avg_log2FC >= log2(1.001),"up.in.tumor",
                              ifelse(pct.2 >= 0.2 & avg_log2FC <= log2(1/1.001),"up.in.metastasis","not"))) %>% 
  filter(DEGs_status != "not")


DEGs_primary_up <- DEGs_filter[DEGs_filter$DEGs_status == "up.in.tumor",]
DEGs_primary_up_count <- as.data.frame(table(DEGs_primary_up$gene_name))
colnames(DEGs_primary_up_count) <- c("gene_name","num")
DEGs_primary_up_count <- arrange(DEGs_primary_up_count,desc(num))
DEGs_metastasis_up <- DEGs_filter[DEGs_filter$DEGs_status == "up.in.metastasis",]
DEGs_metastasis_up_count <- as.data.frame(table(DEGs_metastasis_up$gene_name))
colnames(DEGs_metastasis_up_count) <- c("gene_name","num")
DEGs_metastasis_up_count <- arrange(DEGs_metastasis_up_count,desc(num))


# data plot ---------------------------------------------------------------
primary_up_list <- c("BAIAP2","BHLHE40","EPHA2","EZR","FRMD4B","ITGA2","ITGAV","LAMC2","LIMA1","LMNA","MAFF","RTN4","TACSTD2","TIPARP","TSPAN1")
metastasis_up_list <- c("ANAPC11","ATOX1","CAMK2N1","CCT3","CIAO2B","CKLF","PHB1","OST4","PALM2AKAP2","ATP5MK","ATP5MG","CWC15","GCSH","H2AZ2","SLC25A6")

DEGs_plot <- DEGs_all %>% 
  # filter(gene_name %in% c(DEGs_metastasis_up_count$gene_name[1:60]))
  filter(gene_name %in% c(primary_up_list,metastasis_up_list))
DEGs_plot$gene_name <- factor(DEGs_plot$gene_name,levels = c(primary_up_list,metastasis_up_list))

ggplot(DEGs_plot, aes(x = gene_name, y = avg_log2FC, fill = tumor_code)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = tumor_type_color) +
  labs(x = "", y = "log2FC") + 
  # geom_hline(yintercept = log2(1.2),linetype = "dashed",color = "black") +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/pancancer_metastasis_marker"), width = 300, height = 220)


# enrichment --------------------------------------------------------------
## gene list filter
DEGs_primary_high_gene <- DEGs_primary_up_count %>% 
  filter(num >= 7)
DEGs_metastasis_high_gene <- DEGs_metastasis_up_count %>% 
  filter(num >= 7)

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
  ggplot(aes(x=Description, y=log10pvalue)) + 
  geom_bar(stat = "identity", width = 0.8, colour = 'black') + 
  geom_text(aes(label = Description), y = 0.2,
            fontface = 'bold', size = 4, hjust = 0) + 
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
function_plot(filename_prefix = str_glue("{output_dir}/pancancer_metastasis_enrich"), width = 200, height = 220)

# data save ---------------------------------------------------------------
if(T){
  write_csv(DEGs_all,str_glue("{output_dir}/DEGs_all.csv"))
  write_csv(DEGs_filter,str_glue("{output_dir}/DEGs_filter.csv"))
  write_csv(DEGs_primary_up_count,str_glue("{output_dir}/DEGs_tumor_up_count.csv"))
  write_csv(DEGs_metastasis_up_count,str_glue("{output_dir}/DEGs_metastasis_up_count.csv"))
}








