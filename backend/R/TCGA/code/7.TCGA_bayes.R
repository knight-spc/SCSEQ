####Package####
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(maftools)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(stringr)
library(ggpubr)
library(biomaRt)
library(AnnotationDbi)
library(tidyr)
library(Seurat)
library(BayesPrism)
library(pheatmap)
library(RColorBrewer)
library(optparse)

# 定义命令行参数
option_list <- list(
  make_option(c("-w", "--wd"), type="character", default=NULL, 
              help="请输入工作路径"),
  make_option(c("-s", "--sd"), type="character", default=NULL, 
              help="请输入保存路径"),
  make_option(c("-g", "--gene"), type="character", default=NULL, 
              help="选择要观察的基因"),
  make_option(c("-t", "--time"), type="character", default=NULL, 
              help="请输入时间戳") # 新增时间戳参数
)

# 解析命令行参数
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 检查必要参数
if (is.null(opt$wd)) {
  stop("请提供文件路径 (-i/--wd)")
}
print("==================================================================================================================")
print("正在执行7.TCGA_bayes.R")
####WorkingDirection####
# raw_WD <- '/public-SSD/home/sunpc/scseq/backend/R/TCGA'
####Loaddata####

g <- opt$gene
setwd(opt$sd)
HB <- readRDS("seurat_data.rds") #*HB为用户上传的单细胞数据，需要有注释内容
gdata <- read.csv(paste0('TCGA/res/TCGAmrna_',g,"_tpm.csv"),row.names = 1)

setwd(opt$wd)
TCGAinfo <- read.csv('res/TCGAdir.csv',row.names = 1)
TCGAcount <- readRDS('res/TCGAmrna_count.rds')

# 提取单细胞数据的表达矩阵和细胞类型注释
# sc.dat0 <- as.matrix(HB@assays$RNA@layers$counts)  # Seurat4
sc.dat0 <- as.matrix(LayerData(HB, assay = "RNA", layer = "counts"))   # Seurat5
rownames(sc.dat0) <- rownames(HB)
colnames(sc.dat0) <- colnames(HB)

head(rownames(sc.dat0))
head(colnames(sc.dat0))

#加载bulk RNA-seq数据
names(TCGAcount)
#*i为用户选择的癌症类型，范围为names(TCGAcount)中的一个
i <- names(TCGAcount)[16]

# 转换为矩阵
bulk_expr <- as.matrix(TCGAcount[[i]])
head(rownames(bulk_expr))
head(colnames(bulk_expr))

# 处理bulk RNA-seq数据中的负值
bulk_expr[bulk_expr < 0] <- 0
cell.type.labels <- as.character(HB@meta.data[["Annotation"]])#*单细胞数据中注释A
# cell.state.labels <- as.character(HB@meta.data[["AnnotationSubtype"]])#*单细胞数据中注释B----------------------------------------------
cell.state.labels <- cell.type.labels

#QC
common_genes <- intersect(rownames(bulk_expr), rownames(sc.dat0))
bk.dat <- as.matrix(t(bulk_expr[common_genes,]))
sc.dat <- as.matrix(t(sc.dat0[common_genes,]))

# 检查是否存在NA值
any(is.na(bk.dat))
any(is.na(sc.dat))

# 移除包含NA值的行或列
bk.dat <- bk.dat[, colSums(is.na(bk.dat)) == 0]
sc.dat <- sc.dat[, colSums(is.na(sc.dat)) == 0]

dim(bk.dat)
dim(sc.dat)

# 确保两个数据集的列名依然一致
common_genes1 <- intersect(colnames(bk.dat), colnames(sc.dat))
bk.dat <- bk.dat[, common_genes1]
sc.dat <- sc.dat[, common_genes1]
dim(bk.dat)
dim(sc.dat)


setwd(opt$sd)
setwd("TCGA")

#####*以下所有QC的pdf全都自己重新设置地址保存####
pdf(paste('01.QC.pdf',sep = '/'),width = 8,height=8)
plot.cor.phi (input=sc.dat,
              input.labels=cell.state.labels,
              title="cell state correlation",
              #specify pdf.prefix if need to output to pdf
              #pdf.prefix="gbm.cor.cs", 
              cexRow=0.2, cexCol=0.2,
              margins=c(2,2))

plot.cor.phi (input=sc.dat, 
              input.labels=cell.type.labels, 
              title="cell type correlation",
              #specify pdf.prefix if need to output to pdf
              #pdf.prefix="gbm.cor.ct",
              cexRow=0.5, cexCol=0.5,
)
dev.off()

#单细胞数据的离群基因
pdf(paste('02.QC.pdf',sep = '/'),width = 8,height=8)
sc.stat <- plot.scRNA.outlier(
  input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE #return the data used for plotting. 
  #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
)
dev.off()
head(sc.stat)

#bulk数据的离群基因
bk.dat2 <- bk.dat#2^bk.dat#
pdf(paste('03.QC.pdf',sep = '/'),width = 8,height=8)
bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat2,#make sure the colnames are gene symbol or ENSMEBL ID 
  sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE
  #pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
)
dev.off()
#> EMSEMBLE IDs detected.
head(bk.stat)#查看

#过滤异常基因
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs", 
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM",
                                                "MALAT1","chrX","chrY","hb","act"),
                                  exp.cells=5)

#检查不同类型基因表达的一致性
pdf(paste('04.QC.pdf',sep = '/'),width = 8,height=8)
plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = bk.dat2
                 #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
)
dev.off()
#选择相关性最高的组别
sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
                                         gene.type = "protein_coding")#12155


# 创建BayesPrism对象
summary(sc.dat.filtered.pc[,1:5])
summary(bk.dat[,1:5])

myPrism2 <- new.prism(
  reference=sc.dat.filtered.pc, 
  mixture=bk.dat2,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.state.labels,
  key=NULL,# 
  outlier.cut=0.01,
  outlier.fraction=0.1,
)
#> Number of outlier genes filtered from mixture = 6 
#> Aligning reference and mixture... 
#> Nornalizing reference...


# 运行BayesPrism反卷积
bp.res2 <- run.prism(prism = myPrism2, n.cores=50)

saveRDS(bp.res2,file = 'TCGAbayes.rds')

# 查看和分析结果
#提取细胞类型
thetaTcga <- as.data.frame(get.fraction(bp=bp.res2,
                                        which.theta="first",
                                        state.or.type="state"))#*对应56行注释类型的结果
thetaTcga <- as.data.frame(get.fraction(bp=bp.res2,
                                        which.theta="first",
                                        state.or.type="type"))#*对应57行注释类型的结果
#*上述两个表thetaTcga、thetaTcga以图表的形式呈现
#####plot####

x <- thetaTcga#thetaTcga(*选择其中一个)
rownames(gdata) <- gdata$Index
x <- cbind(gdata[rownames(x),g],x)
colnames(x)[1] <- g

####heatmap####

# 设置颜色方案
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)

# 数据标准化 (Z-score 标准化)
thetaT_scaled <- thetaTcga[,]

thetaT_scaled <- t(scale(thetaT_scaled))
thetaT_scaled <- scale(thetaT_scaled)
thetaT_scaled <- t(thetaT_scaled )
# 画热图
# 5. 画热图 
pdf(paste('7.TCGA_bayes.pdf',sep = '/'),width = 8,height=8)
pheatmap(thetaT_scaled, 
         color = heatmap_colors, 
         cluster_rows = T, # 不对行聚类 
         cluster_cols = F, # 只对列聚类 
         show_rownames = F, 
         show_colnames = T, 
         main = "Kuffer Distribution in LIHC",
         angle_col = 45
)####*展示的第一张####
dev.off()


####*展示的第二张图####

# 定义一个函数来计算单个组的相关性
calculate_correlation <- function(subset_data, cg_columns,method,g) {
  results <- lapply(cg_columns, function(cg) {
    # 计算 HSPA1B 和当前 cg 的相关性
    test_result <- cor.test(subset_data[[g]], subset_data[[cg]], method = method)#c('spearman','pearson','kendall)
    data.frame(
      cg = cg,
      Correlation = test_result$estimate,
      P_Value = test_result$p.value
    )
  })
  # 合并结果
  do.call(rbind, results)
}
# 按 Cancer 类型分组并计算相关性
lis <- colnames(x)[-1]
correlation_results <- x %>% na.omit() %>%
  group_modify(~ calculate_correlation(.x, cg_columns=lis,method = 'pearson',g=g)) %>%
  ungroup()

correlation_results$res = case_when(
  correlation_results$P_Value < 0.001 ~ "***",
  correlation_results$P_Value < 0.01  ~ "**",
  correlation_results$P_Value < 0.05  ~ "*",
  is.na(correlation_results$P_Value)  ~ "",
  correlation_results$P_Value > 0.05 ~ ""
)

# 如果需要保存为 CSV 文件
write.csv(correlation_results, paste0('res/',g,"_bayes_imm&TPM_", opt$time, ".csv"), row.names = FALSE)

#####plot#### 

# correlation_results <- read.csv(paste0('res/',g,"_bayes_imm&TPM.csv"))#*先展示表

print("7.TCGA_bayes.R绘图完成")
# 清理所有对象
rm(list = ls())
# 回收内存
gc()
