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
(.packages())
setwd(' ')

function_plot <- function(filename_prefix, width, height){
  ggsave(filename = str_c(filename_prefix, '.pdf'), plot = last_plot(), device = 'pdf', width = width, height = height, units = 'mm', dpi = 500)
  # ggsave(filename = str_c(filename_prefix, '.svg'), plot = last_plot(), device = 'svg', width = width, height = height, units = 'mm', dpi = 500)
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
stRNA_dir <- " "
output_dir <- " "
hub_MP <- c("Cell.Cycle.G2M","Cell.Cycle.G1S","Cell.Cycle.HMG","Chromatin","Cell.Cycle.single.nucleus","Stress","Hypoxia",
            "Stress.in.vitro","Protein.maturation","EMT.I","EMT.II","EMT.VI","MHC","Epithelial.Senescence","MYC",
            "Respriration","Secreted.I","Secreted.II","Cilia.I","Cilia.II","PDAC.classcial","Alveolar","Colon.related",
            "PDAC.related","Lymphocyte.Activation","Cell.Signaling","Coagulation","Neuron")

# data process ------------------------------------------------------------
st_file <- list.files(stRNA_dir)

for(MP in hub_MP){
  # MP <- "Hypoxia"
  
  MP_file <- str_glue("{output_dir}/{MP}/")
  unlink(MP_file, recursive = T, force = T, expand = T)
  dir.create(MP_file, recursive = T)
  
  for(file in 1:length(st_file)){
    # file <- 5
    print(str_glue("################################ {st_file[file]} Run AUcell ###################################"))
    set.seed(1234)
    scRNA_list <- list()
    input_file <- str_glue("{stRNA_dir}/{st_file[file]}")
    file_name <- list.files(input_file)
    
    ## dir create
    output_file <- str_glue("{MP_file}/{st_file[file]}")
    unlink(output_file, recursive = T, force = T, expand = T)
    dir.create(output_file, recursive = T)
    
    for(i in 1:length(file_name)){
      # i<- 1
      print(str_glue("           -------------- {file_name[i]} Sample Run"))
      ## data read
      stRNA_list <- readRDS(str_glue("{stRNA_dir}/{st_file[file]}/{file_name[i]}"))
      
      if(!MP %in% colnames(stRNA_list@meta.data)){
        next
      }
      
      ## data transform
      names(stRNA_list@images) <- "slice1"
      stRNA_list@images[["slice1"]]@coordinates$tissue <- as.integer(stRNA_list@images[["slice1"]]@coordinates$tissue)
      stRNA_list@images[["slice1"]]@coordinates$row <- as.integer(stRNA_list@images[["slice1"]]@coordinates$row)
      stRNA_list@images[["slice1"]]@coordinates$col <- as.integer(stRNA_list@images[["slice1"]]@coordinates$col)
      stRNA_list@images[["slice1"]]@coordinates$imagerow <- as.integer(stRNA_list@images[["slice1"]]@coordinates$imagerow)
      stRNA_list@images[["slice1"]]@coordinates$imagecol <- as.integer(stRNA_list@images[["slice1"]]@coordinates$imagecol)
      
      ## tumor spot filter
      stRNA_list <- subset(stRNA_list,copykat_result == "malignant")
      
      ## data plot
      SpatialFeaturePlot(stRNA_list,features = MP,alpha = c(0.5, 1),pt.size.factor = nrow(stRNA_list@images$slice1@image)/390) + NoLegend()
      function_plot(filename_prefix = str_glue("{output_file}/{file_name[i]}"), width = 200, height = 200)
      
      
    }
  }
}









