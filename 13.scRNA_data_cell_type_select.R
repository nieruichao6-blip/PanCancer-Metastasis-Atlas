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
library("dior")
library("sceasy")
library("reticulate")
library("ggalluvial")
library("networkD3")
(.packages())
setwd(' ')
display.brewer.all()


function_plot <- function(filename_prefix, width, height){
  ggsave(filename = str_c(filename_prefix, '.pdf'), plot = last_plot(), device = 'pdf', width = width, height = height, units = 'mm', dpi = 500)
  ggsave(filename = str_c(filename_prefix, '.png'), plot = last_plot(), device = 'png', width = width, height = height, units = 'mm', dpi = 2000)
  ggsave(filename = str_c(filename_prefix, '.svg'), plot = last_plot(), device = 'svg', width = width, height = height, units = 'mm', dpi = 500)
  try(ggsave(filename = str_c(filename_prefix, '.wmf'), plot = last_plot(), device = 'wmf', width = width, height = height, units = 'mm', dpi = 500),
      silent = T)
}


## plot theme
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

# data color --------------------------------------------------------------
cell_type_color <- c("Epithelia_tumor"="#E41A1C",
                     "Melanoma_tumor"="#E41A1C",
                     "Osteoblastic_tumor"="#E41A1C",
                     "Epithelia"="#A6CEE3",
                     "Fibroblast"="#FDBF6F",
                     "Endothelia"="#B2DF8A",
                     "Acinar_cell"="#984EA3",
                     "SMC&Pericyte"="#33A02C",
                     "T&NK_cell"="#1F78B4",
                     "B_cell"="#FF7F00",
                     "Plasma_cell"="#FB9A99",
                     "Myeloid_cell"="#CAB2D6",
                     "Mast_cell"="#6A3D9A",
                     "Neutrophils"="#FFFF99",
                     "Hepatocyte"="#9ec9e1",
                     "Melanoma_cell"="#E66F00",
                     "Neuron"="#1d92c0",
                     "Glial_cell"="#B15928",
                     "Osteoblastic_cell"="#fcc5c1",
                     "Alpha_cell"="#42aa5e",
                     "Bela_cell"="#c22b86"
)



# data dir ----------------------------------------------------------------
input_dir <- " "
output_dir <- " "


# hub ref -----------------------------------------------------------------
hub_cell_type <- c("T_cell",
                   "NK_cell"
                   )
# hub_cell_type <- c(
#                    "Epithelia",
#                    "Epithelia_tumor",
#                    "Melanoma_tumor",
#                    "Osteoblastic_tumor",
#                    "Fibroblast",      # "SMC","Pericyte",
#                    "Endothelia",
#                    "T_cell",
#                    "NK_cell",
#                    "B_cell",
#                    "Plasma_cell",
#                    "Myeloid_cell",
#                    "Neutrophils"
#                    )
num_cutoff <- 50



# data process ------------------------------------------------------------
cancer_file <- list.files(input_dir)

##data stat
cell_num_stat <- data.frame(
  tumor_type = character(),
  cell_type = character(),
  num = numeric()
)

## data run
for(file in 1:length(cancer_file)){
  # file <- 4
  print(str_glue("############################## {cancer_file[file]} Pre Processing ############################"))
  set.seed(1234)
  
  ## dir create
  output_file <- str_glue("{output_dir}/{cancer_file[file]}/")
  # unlink(output_file, recursive = T, force = T, expand = T)
  # dir.create(output_file, recursive = T)
  
  ## data read
  scRNA_current <- readRDS(str_glue("{input_dir}/{cancer_file[file]}/scRNA_infercnv_combine.RData"))
  
  ## cell type filter
  for(celltype in hub_cell_type){
    if(celltype %in% scRNA_current@meta.data$new_cell_type){
      print(str_glue("          -------------- {celltype} select"))
      
      ### cell type select
      if(celltype == "Fibroblast"){
        scRNA_hub_cell <- subset(scRNA_current, new_cell_type %in% c("Fibroblast","SMC&Pericyte"))
      }else if(celltype == "B_cell"){
        scRNA_hub_cell <- subset(scRNA_current, new_cell_type %in% c("B_cell"))
      }else if(celltype == "Myeloid_cell"){
        scRNA_hub_cell <- subset(scRNA_current, new_cell_type %in% c("Myeloid_cell","Mast_cell"))
      }else if(celltype == "Epithelia"){
        scRNA_hub_cell <- subset(scRNA_current, new_cell_type %in% c("Epithelia"))
      }else if(celltype == "Epithelia_tumor" & cancer_file[file] != "liver_healthy" & cancer_file[file] != "lymph_node_healthy"){
        scRNA_hub_cell <- subset(scRNA_current, new_cell_type %in% c("Epithelia_tumor"))
      }else if(celltype == "Melanoma_cell"){
        scRNA_hub_cell <- subset(scRNA_current, new_cell_type %in% c("Melanoma_tumor"))
      }else if(celltype == "Osteoblastic_cell"){
        scRNA_hub_cell <- subset(scRNA_current, new_cell_type %in% c("Osteoblastic_tumor"))
      }else if(celltype == "T_cell"){
        scRNA_hub_cell <- subset(scRNA_current, new_cell_type %in% c("T_cell"))
      }else if(celltype == "NK_cell"){
        scRNA_hub_cell <- subset(scRNA_current, new_cell_type %in% c("NK_cell"))
      }else{
        scRNA_hub_cell <- subset(scRNA_current, new_cell_type %in% celltype)
      }
      
      ### sample filter
      sample_num <- as.data.frame(table(scRNA_hub_cell$sample_ID))
      colnames(sample_num) <- c("sample","num")
      sample_num <- sample_num %>% 
        filter(num >= num_cutoff)
      if(nrow(sample_num) >= 1){
        scRNA_hub_cell <- subset(scRNA_hub_cell, sample_ID %in% sample_num$sample)
        
        ### data stat
        cell_num <- nrow(scRNA_hub_cell@meta.data)
        current_stat <- data.frame(
          tumor_type = cancer_file[file],
          cell_type = celltype,
          num = cell_num
        )
        cell_num_stat <- rbind(cell_num_stat,current_stat)
        
        
        ### data save
        saveRDS(scRNA_hub_cell,str_glue("{output_file}/scRNA_{celltype}.RData"))
        dior::write_h5(scRNA_hub_cell, file=str_glue("{output_file}/scRNA_{celltype}.h5"),
                       assay.name = "RNA",object.type = 'seurat')
        
        ### space free
        remove(scRNA_hub_cell)
        gc()
        
        
      }

    }
    
  }
  
  ## space free
  remove(scRNA_current)
  gc()

}


# data save ---------------------------------------------------------------
# cell_num_stat_filter <- cell_num_stat[!(cell_num_stat$cell_type == "Epithelia" & cell_num_stat$tumor_type == "liver_healthy"),]
# cell_num_stat_filter <- cell_num_stat_filter[!(cell_num_stat_filter$cell_type == "Epithelia" & cell_num_stat_filter$tumor_type == "lymph_node_healthy"),]

write_csv(cell_num_stat,str_glue("{output_dir}/cell_num_stat.csv"))

# data plot ---------------------------------------------------------------

ggplot(cell_num_stat, aes(tumor_type, num, fill = cell_type)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "", y = "Number of cells") +
  scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = cell_type_color) +
  nrc_theme
function_plot(filename_prefix = str_glue("{output_dir}/scRNA_useful_cell"), width = 300, height = 200)




