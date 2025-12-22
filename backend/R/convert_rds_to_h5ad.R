# 用于将rds数据转换为h5ad格式，方便后续celltypist进行注释

library(Seurat)
library(fgsea)
library(SeuratDisk)
library(optparse)

# 定义命令行参数
option_list <- list(
  make_option(c("-w", "--wd"), type="character", default=NULL, 
              help="请输入工作路径"),
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="请输入RDS文件名"),
  make_option(c("-s", "--hs"), type="character", default=NULL, 
              help="请输入h5Seurat文件名")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 检查必要参数
if (is.null(opt$wd)) {
  stop("请提供数据路径 (-i/--wd)")
}


convert_rds_to_h5ad <- function(wd_path, file_name, hs_name) {
    setwd(wd_path)
    data <- readRDS(file_name)
    data[["RNA"]] <- as(data[["RNA"]], "Assay")  # 将RNA数据转换为Assay对象，Seurat5版本后需要这一步
    DefaultAssay(data) <- "RNA"
    print(hs_name)

    SaveH5Seurat(data,filename = hs_name)
    Convert(hs_name,dest = 'h5ad')
    rm(data)
    gc()
}

wd_path = opt$wd
file_name = opt$file
hs_name = opt$hs
convert_rds_to_h5ad(wd_path,file_name,hs_name)