library(dplyr) # 1.0.10
library(Seurat)
library(optparse)
library(glue)

# 定义命令行参数
option_list <- list(
  make_option(c("-w", "--wd"), type="character", default=NULL, 
              help="请输入工作路径"),  
  make_option(c("-s", "--se"), type="character", default=NULL, 
              help="请输入工作路径"),
  make_option(c("-c", "--ca"), type="character", default=NULL, 
              help="请输入umap路径"),
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
# 读取 RDS 文件
num_tab <- table( data@meta.data[[opt$se]], data@meta.data[[opt$ca]] )
freq_tab <- num_tab
freq_df <- as.data.frame.matrix(freq_tab)

scale_diagram_path = opt$file
write.csv(freq_df, scale_diagram_path, row.names = TRUE)

print(glue("已保存到{scale_diagram_path}"))
print("保存完毕")
rm(list = ls())
gc()