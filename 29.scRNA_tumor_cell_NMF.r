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
library("glmGamPoi") 
library("paletteer")
library("NMF")
library("BiocParallel")
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
origin_dir <- " "
cluster_dir <- " "
output_dir <- " "

# data read ----------------------------------------------------------------
## high var gene
high_var_gene <- read_csv(str_glue("{cluster_dir}/high_variable_gene.csv"))

## exp data
scRNA_tumor_transform <- read_h5(file = str_glue("{cluster_dir}/scRNA_cell_type_cluster_only.h5"),
                         assay.name = 'RNA', 
                         target.object = 'seurat')
scRNA_tumor_origin <- read_h5(file = str_glue("{origin_dir}/scRNA_gene_filter.h5"),
                         assay.name = 'RNA', 
                         target.object = 'seurat')                        
# scRNA_tumor@assays$RNA@var.features <- high_var_gene$gene_name
# scRNA_tumor <- FindVariableFeatures(scRNA_tumor, nfeatures = 3000)

## copy data
scRNA_tumor_origin@reductions <- scRNA_tumor_transform@reductions

## space free
remove(scRNA_tumor_transform)
gc()

# sample ID split ----------------------------------------------------------------
scRNA_tumor.list <- SplitObject(scRNA_tumor_origin, split.by = "sample_ID")

# pre process ----------------------------------------------------------------
## nohub gene
mito.genes <- grep(pattern = "^MT-", x = rownames(scRNA_tumor_origin), value = TRUE)
rbl.genes <- grep(pattern = "^RB-", x = rownames(scRNA_tumor_origin), value = TRUE)
rsl.genes <- grep(pattern = "^RS-", x = rownames(scRNA_tumor_origin), value = TRUE)  
rpl.genes <- grep(pattern = "^RPL-", x = rownames(scRNA_tumor_origin), value = TRUE)  
rbl.genes <- grep(pattern = "^RBL-", x = rownames(scRNA_tumor_origin), value = TRUE)  
rps.genes <- grep(pattern = "^RPS-", x = rownames(scRNA_tumor_origin), value = TRUE)  
rbs.genes <- grep(pattern = "^RBS-", x = rownames(scRNA_tumor_origin), value = TRUE)  
rbl1.genes <- grep(pattern = "^RB", x = rownames(scRNA_tumor_origin), value = TRUE)  
rsl1.genes <- grep(pattern = "^RS", x = rownames(scRNA_tumor_origin), value = TRUE)  
rpl1.genes <- grep(pattern = "^RPL", x = rownames(scRNA_tumor_origin), value = TRUE)  
rbl1.genes <- grep(pattern = "^RBL", x = rownames(scRNA_tumor_origin), value = TRUE)  
rps1.genes <- grep(pattern = "^RPS", x = rownames(scRNA_tumor_origin), value = TRUE)  
rbs1.genes <- grep(pattern = "^RBS", x = rownames(scRNA_tumor_origin), value = TRUE)
hub_gene <- rownames(scRNA_tumor_origin)[!rownames(scRNA_tumor_origin) %in% c(mito.genes,rbl.genes,rsl.genes,rpl.genes,rbl.genes,rps.genes,
                                                  rbs.genes,rbl1.genes,rsl1.genes,rpl1.genes,rbl1.genes,rps1.genes,rbs1.genes)]
## nor transform
for(sample in names(scRNA_tumor.list)){
  # sample <- "HRS1028865"
  print(str_glue("   -------------{sample}"))
  # scRNA_tumor.list[[sample]] <- subset(scRNA_tumor.list[[sample]], features = hub_gene)
  scRNA_tumor.list[[sample]] <- NormalizeData(scRNA_tumor.list[[sample]], normalization.method = "LogNormalize", scale.factor = 10000)
  scRNA_tumor.list[[sample]] <- FindVariableFeatures(scRNA_tumor.list[[sample]], selection.method = "vst", nfeatures = 5000)
  scRNA_tumor.list[[sample]] <- ScaleData(scRNA_tumor.list[[sample]], features = rownames(scRNA_tumor.list[[sample]]),do.center = TRUE)
  # scRNA_tumor.list[[sample]] <- SCTransform(scRNA_tumor.list[[sample]], vst.flavor = "v2", method = "glmGamPoi",verbose = F, 
  #                                           return.only.var.genes = F, do.scale = TRUE, do.center = FALSE)
  
  # exp_matrix <- scRNA_tumor.list[[sample]]@assays$SCT@scale.data
  # exp_matrix <- scRNA_tumor.list[[sample]]@assays$RNA@scale.data
  # exp_matrix[exp_matrix < 0] <- 0
  # exp_matrix <- exp_matrix[rowSums(exp_matrix) > 0,]
  # scRNA_tumor.list[[sample]] <- subset(scRNA_tumor.list[[sample]], features = rownames(exp_matrix))
  # scRNA_tumor.list[[sample]]@assays$RNA@scale.data <- exp_matrix
  

  gc()
}
all_gene <- rownames(scRNA_tumor_origin)

## space free
remove(scRNA_tumor_origin)
gc()

# run NMF ----------------------------------------------------------------
geneNMF.programs <- multiNMF(scRNA_tumor.list, assay="RNA", slot="scale.data", k=4:9,  min.cells.per.sample = 50, nfeatures = 10000,
                            center = FALSE, scale = FALSE, seed = 123)
# geneNMF.programs <- multiNMF(scRNA_tumor.list, assay="SCT", slot="scale.data", k=10, min.cells.per.sample = 50, nfeatures = 10000,
#                             center = FALSE, scale = FALSE, seed = 123)
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

for(MP in c(12:18)){
  geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs, metric = "jaccard", 
                                        nMP=MP, weight.explained = 0.8, max.genes=50,
                                        min.confidence = 0.5)
  geneNMF.metaprograms$programs.similarity <- geneNMF.metaprograms$programs.similarity * 100

  pdf(str_glue("{output_dir}/NMF_confirm_{MP}.pdf"),width = 12,height = 9)
  plotMetaPrograms(geneNMF.metaprograms, similarity.cutoff = c(0.1,0.6))
  dev.off()
}

# NMF confirm ----------------------------------------------------------------
gc()
geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs, metric = "jaccard",
                                        nMP=30, weight.explained = 0.8, max.genes=50,
                                        min.confidence = 0.25)
geneNMF.metaprograms_copy <- geneNMF.metaprograms
geneNMF.metaprograms_copy$programs.similarity <- geneNMF.metaprograms_copy$programs.similarity * 100

# NMF plot ----------------------------------------------------------------
pdf(str_glue("{output_dir}/NMF_final_heatmap.pdf"),width = 12,height = 9)
plotMetaPrograms(geneNMF.metaprograms, 
                 similarity.cutoff = c(0.04, 0.25),downsample = 2000,
                # palette = c('white','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b')
                palette = colorRampPalette(c('white', '#FCFFC9', '#F4E228', '#E61A13', '#5A1833', '#1D0B14'))(50)
                # palette = custom_magma
                )
dev.off()

# NMF module ano ----------------------------------------------------------------
## enrich
enrich_GO <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program){runGSEA(program, universe=all_gene,
                    category = "C5", subcategory = "GO:BP")})
enrich_Hallmarker <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program){runGSEA(program, universe=all_gene,
                            category = "H")})

## scRNA score
mp.genes <- geneNMF.metaprograms$metaprograms.genes
scRNA_tumor <- AddModuleScore_UCell(scRNA_tumor,features = mp.genes,ncores=12,name = "")

# data save ----------------------------------------------------------------
saveRDS(scRNA_tumor.list,str_glue("{output_dir}/scRNA_sample_split.RData"))
saveRDS(geneNMF.metaprograms,str_glue("{output_dir}/NMF_MP_data.RData"))
saveRDS(scRNA_tumor,str_glue("{output_dir}/scRNA_tumor_NMF.RData"))
dior::write_h5(scRNA_tumor, file=str_glue("{output_dir}/scRNA_tumor_NMF.h5"),
                assay.name = "RNA",object.type = 'seurat')


















################################ test
scRNA_tumor.list <- SplitObject(scRNA_tumor_origin, split.by = "sample_ID")
patient <- names(scRNA_tumor.list) # 获取患者ID列表# 对每个患者数据执行NMF分析
topn <- 10000  # 保留变异最大的前10000个基因
rank <- 10  # NMF分解的维度（提取10个特征程序）
for(i in patient){
  # 创建结果存储目录（若不存在则新建）  
  if (!dir.exists(paste(" ", i, sep = "/"))){
    dir.create(paste(" ", i, sep = "/")) }
  # 1. 数据过滤：移除线粒体基因和核糖体基因（避免非特异性表达干扰）  
  snRNA <- scRNA_tumor.list[[i]]  # 获取当前患者数据
  # 定义需要排除的基因（线粒体基因及核糖体相关基因）  
  mito.genes <- grep(pattern = "^MT-", x = rownames(snRNA), value = TRUE)
  rbl.genes <- grep(pattern = "^RB-", x = rownames(snRNA), value = TRUE)
  rsl.genes <- grep(pattern = "^RS-", x = rownames(snRNA), value = TRUE)  
  rpl.genes <- grep(pattern = "^RPL-", x = rownames(snRNA), value = TRUE)  
  rbl.genes <- grep(pattern = "^RBL-", x = rownames(snRNA), value = TRUE)  
  rps.genes <- grep(pattern = "^RPS-", x = rownames(snRNA), value = TRUE)  
  rbs.genes <- grep(pattern = "^RBS-", x = rownames(snRNA), value = TRUE)  
  rbl1.genes <- grep(pattern = "^RB", x = rownames(snRNA), value = TRUE)  
  rsl1.genes <- grep(pattern = "^RS", x = rownames(snRNA), value = TRUE)  
  rpl1.genes <- grep(pattern = "^RPL", x = rownames(snRNA), value = TRUE)  
  rbl1.genes <- grep(pattern = "^RBL", x = rownames(snRNA), value = TRUE)  
  rps1.genes <- grep(pattern = "^RPS", x = rownames(snRNA), value = TRUE)  
  rbs1.genes <- grep(pattern = "^RBS", x = rownames(snRNA), value = TRUE)

  # 过滤基因并更新Seurat对象  
  counts <- snRNA@assays$RNA@counts  
  counts <- counts[-(which(rownames(counts) %in% c(mito.genes,rbl.genes,rsl.genes,rpl.genes,rbl.genes,rps.genes,
                                                  rbs.genes,rbl1.genes,rsl1.genes,rpl1.genes,rbl1.genes,rps1.genes,rbs1.genes))),]
  counts <- counts[rowSums(counts) > 0,]  # 保留至少在一个细胞中表达的基因  
  snRNA <- subset(snRNA, features = rownames(counts))
  # 2. 数据标准化：使用SCTransform（适用于单细胞数据的标准化方法）  
  snRNA <- SCTransform(snRNA, vst.flavor = "v2", method = "glmGamPoi",verbose = F, return.only.var.genes = F)
  # 3. 构建NMF输入矩阵（确保非负性和高信息基因）  
  nmf_mat <- GetAssayData(snRNA, slot = "scale.data")  # 提取标准化后的缩放数据  
  nmf_mat[nmf_mat < 0] <- 0  # NMF要求非负矩阵，负值设为0  
  nmf_mat <- nmf_mat[rowSums(nmf_mat) > 0,]  # 移除无表达的基因  
  gene_sd <- apply(nmf_mat, 1, sd)  # 计算基因表达标准差（衡量变异程度）  
  top_gene <- gene_sd[order(gene_sd, decreasing = T)][1:topn]  # 取变异最大的前10000个基因  
  nmf_mat <- nmf_mat[rownames(nmf_mat) %in% names(top_gene),]  # 筛选高变异基因  
  ina <- which(colSums(is.na(nmf_mat)) == 0)  # 移除含缺失值的细胞  
  nmf_mat <- nmf_mat[,ina]
  # 运行NMF分析（使用brunet算法，nndsvd初始化提高稳定性）  
  res <- nmf(nmf_mat, rank = rank, method = "brunet", seed = "nndsvd", .options = "p20v1")
  # 提取特征基因并保存结果  
  signature <- NMF::basis(res)  # 获取基因特征矩阵（basis矩阵）  
  colnames(signature) <- paste(i, 1:rank, sep = "_")  # 命名特征（患者ID_特征序号）  
  signature <- as.data.frame(signature)
  # 保存特征基因和NMF结果  
  write.table(signature, paste0(" ", i, "/signature_", topn, ".txt"), sep = "\t")
  saveRDS(res, file = paste0(" ", i, "/result_", rank, ".rds"))
  print(paste0("NMF for ", i, " is done!"))  # 提示当前患者处理完成
}

opn <- 10000 # 匹配之前的参数
ranks <- 10   # NMF维度
topRank <- 100  # 每个特征程序保留的top基因数量
programG <- list()  # 存储特征基因列表

# 提取每个患者的特征基因
for (i in seq_along(patient)){ 
  filedir <- paste0("./nmf_res/", patient[i], "/signature_", topn, ".txt")  # 特征基因文件路径  
  geneloading <- read.table(filedir, header = T, sep = "\t")  # 读取基因载荷矩阵  # 标记每个基因的最高载荷对应的特征程序  
  geneloading$maxC <- apply(geneloading, 1, which.max) %>% paste0(patient[i], "_", .)
  # 格式化数据并提取每个程序的top100基因  
  topgenelist <- rownames_to_column(geneloading, var = "gene") %>%   
    pivot_longer(., cols = starts_with(c("P", "R")),
    names_to = "program", values_to = "loading")
  topgenelist <- dplyr::filter(topgenelist, maxC == program) %>%
    group_by(maxC) %>% top_n(n = topRank, wt = loading)  # 按载荷取top100基因  
    topgenelist <- split(topgenelist$gene, topgenelist$maxC)  # 按程序拆分  
    programG <- c(programG, topgenelist)  # 整合到列表
}

# 统一过滤原始数据（与NMF输入一致）
snRNA <- snRNA_tumor  # 原始肿瘤单细胞数据# 再次过滤线粒体和核糖体基因（确保基因集一致）
mito.genes <- grep(pattern = "^MT-", x = rownames(snRNA), value = TRUE)
rbl.genes <- grep(pattern = "^RB-", x = rownames(snRNA), value = TRUE)
rsl.genes <- grep(pattern = "^RS-", x = rownames(snRNA), value = TRUE)
rpl.genes <- grep(pattern = "^RPL-", x = rownames(snRNA), value = TRUE)
rbl.genes <- grep(pattern = "^RBL-", x = rownames(snRNA), value = TRUE)
rps.genes <- grep(pattern = "^RPS-", x = rownames(snRNA), value = TRUE)
rbs.genes <- grep(pattern = "^RBS-", x = rownames(snRNA), value = TRUE)
rbl1.genes <- grep(pattern = "^RB", x = rownames(snRNA), value = TRUE)
rsl1.genes <- grep(pattern = "^RS", x = rownames(snRNA), value = TRUE)
rpl1.genes <- grep(pattern = "^RPL", x = rownames(snRNA), value = TRUE)
rbl1.genes <- grep(pattern = "^RBL", x = rownames(snRNA), value = TRUE)
rps1.genes <- grep(pattern = "^RPS", x = rownames(snRNA), value = TRUE)
rbs1.genes <- grep(pattern = "^RBS", x = rownames(snRNA), value = TRUE)
counts <- snRNA@assays$RNA@countscounts <- counts[-(which(rownames(counts) %in% c(mito.genes,rbl.genes,rsl.genes,rpl.genes,
                                                                                  rbl.genes,rps.genes,rbs.genes,rbl1.genes,
                                                                                  rsl1.genes,rpl1.genes,rbl1.genes,rps1.genes,rbs1.genes))),]counts <- counts[rowSums(counts) > 0,]
snRNA <- subset(snRNA, features = rownames(counts))# 计算每个细胞的特征程序评分（基于特征基因集）
snRNA <- SCTransform(snRNA, vst.flavor = "v2", method = "glmGamPoi", verbose = F, return.only.var.genes = F)
exp_mat <- GetAssayData(snRNA, slot = "data")  # 提取表达矩阵
exp_mat <- as.matrix(exp_mat)# 使用pagoda2的score.cells.puram计算评分
score_list <- list()
score_list <- lapply(programG, function(x){  score <- score.cells.puram(data = t(exp_mat), signature = x) return(score)})

# 整合评分为矩阵（行：细胞，列：特征程序）
score_mat <- sapply(score_list, FUN = function(x) x, simplify = T)
saveRDS(score_mat, file = 'score_mat.rds')  # 保存评分矩阵
