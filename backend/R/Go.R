# 加载必要的包
library(clusterProfiler)
library(org.Hs.eg.db)  # 这里以人类基因为例
library(glue)
library(dplyr) # 1.0.10

library(optparse)

# 定义命令行参数
option_list <- list(
  make_option(c("-w", "--wd"), type="character", default=NULL, 
              help="请输入工作路径"),
  make_option(c("-c", "--cluster"), type="character", default=NULL, 
              help="请输入注释模型"),
  make_option(c("-m", "--marker"), type="character", default=NULL, 
              help="请输入marker"),
  make_option(c("-t", "--time"), type="character", default=NULL, 
              help="请输入时间戳")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 检查必要参数
if (is.null(opt$wd)) {
  stop("请提供数据路径 (-i/--wd)")
}


getBarPlot <- function(work_path, cluster, timestamp, all_marker) {
    markers = read.csv(all_marker)    
    setwd(work_path)
    # 将基因符号转换为Entrez ID
    gene_symbols <- markers$gene[which(markers$cluster==cluster)]
    gene_entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    # 执行GO富集分析
    go_results <- enrichGO( 
        gene = gene_entrez_ids$ENTREZID,  # 输入的Entrez基因ID
        OrgDb = org.Hs.eg.db,             # 选择基因组数据库
        keyType = "ENTREZID",   # 基因ID的类型
        ont = "BP",             # GO本体：生物过程（BP）、分子功能（MF）或细胞组成（CC）
        pvalueCutoff = 0.05,    # P值的阈值
        pAdjustMethod = "BH",   # P值的校正方法
        qvalueCutoff = 0.2,     # q值的阈值
        minGSSize = 10,         # 富集的最小基因集大小
        maxGSSize = 500,        # 富集的最大基因集大小
        universe = NULL,        # 全基因集合，用于背景分析，默认为NULL
        readable = TRUE,        # 是否将基因ID转换为基因名称
        pool = FALSE           # 是否允许多重检验的池化方法
    )
    go_results = as.data.frame(go_results)
    results = go_results[, c("Count", "Description", "p.adjust")]
    write.csv(results, glue("go_data_{timestamp}.csv"), row.names = FALSE)

    # 清理所有对象
    rm(list = ls())
    # 回收内存
    gc()

    print("完成Go富集分析")

}


work_path = opt$wd
cluster = opt$cluster
timestamp = opt$time
all_marker = opt$marker

getBarPlot(work_path,cluster,timestamp,all_marker)