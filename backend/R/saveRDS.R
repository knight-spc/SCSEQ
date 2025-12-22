library(Seurat)
library(dplyr)
library(optparse)

# 定义命令行参数
option_list <- list(
  make_option(c("-w", "--wd"), type="character", default=NULL, 
              help="请输入工作路径"),
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="请输入umap路径")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 检查必要参数
if (is.null(opt$wd)) {
  stop("请提供文件路径 (-i/--wd)")
}
print('正在读取数据')
setwd(opt$wd)
data <- readRDS("seurat_data.rds")

# 检查HSMM对象中的基因名称
umap_data = read.csv(opt$file)
data@meta.data$celltype<-umap_data$annotation
Idents(data) <- data@meta.data$celltype
saveRDS(data, file = "seurat_data.rds")
print("保存完毕")
rm(list = ls())
gc()