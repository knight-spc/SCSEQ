print('111111111')
library(Seurat) # 4.2.0
library(RColorBrewer) # 1.1-3
library(dplyr) # 1.0.10
library(harmony) # 1.2.0
library(pals) # 1.7
library(ggplot2) # 3.3.6
library(glue)
library(readxl)
library(SingleR)
library(optparse)

print('111111111')

# 定义命令行参数
option_list <- list(
  make_option(c("-w", "--wd"), type="character", default=NULL, 
              help="输入的工作路径"),
  make_option(c("-t", "--time"), type="character", default=NULL, 
              help="请输入时间戳")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
timestamp = opt$time

# 检查必要参数
if (is.null(opt$wd)) {
  stop("请提供文件路径 (-i/--wd)")
}
setwd(opt$wd)

data <- readRDS("seurat_data.rds")
# 确保细胞类型标签已正确设置


print("当前细胞类型：")
# Idents(data) <- data@meta.data$celltype

print(table(Idents(data)))
# data <- JoinLayers(data)


# 计算marker基因
data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
print(data.markers)

# 提取top基因并确保基因不重复
data.markers.top20 = data.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
data.markers.top50 = data.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

# 添加进度显示
print(sprintf("开始计算%d个基因的表达量", length(unique(data.markers.top50$gene))))

marker_data = data.frame(gene = data.markers.top20$gene, 
                p_val = data.markers.top20$p_val, 
                avg_log2FC = data.markers.top20$avg_log2FC,
                pct.1 = data.markers.top20$pct.1,
                pct.2 = data.markers.top20$pct.2,
                p_val_adj = data.markers.top20$p_val_adj,
                cluster = data.markers.top20$cluster,stringsAsFactors = F) 

# 新增: 计算所有基因在各个细胞类型中的表达情况
expression_data <- data@assays$RNA$data
annotations <- Idents(data)
results <- list()

print('开始计算top基因的表达量')
genes <- unique(data.markers.top50$gene)  # 使用unique避免重复基因
total_genes <- length(genes)

for (i in seq_along(genes)) {
    if (i %% 100 == 0) {
        cat(sprintf("Processing gene %d/%d\n", i, total_genes))
    }
    
    gene_of_interest <- genes[i]
    gene_expression <- expression_data[gene_of_interest, ]
    
    for (annot in unique(annotations)) {
        subset_cells <- which(annotations == annot)
        mean_expr <- mean(gene_expression[subset_cells], na.rm = TRUE)
        expr_threshold <- 0.1
        expr_percentage <- mean(gene_expression[subset_cells] > expr_threshold, na.rm = TRUE)
        
        result_df <- data.frame(
            Annotation = annot,
            Gene = gene_of_interest,
            Mean_Expression = mean_expr,
            Expression_Percentage = expr_percentage,
            stringsAsFactors = FALSE
        )
        results <- c(results, list(result_df))
    }
}


# 计算并保存allmarker基因数据
all_marker_data = data.frame(gene = data.markers$gene, 
                    p_val = data.markers$p_val, 
                    avg_log2FC = data.markers$avg_log2FC,
                    pct.1 = data.markers$pct.1,
                    pct.2 = data.markers$pct.2,
                    p_val_adj = data.markers$p_val_adj,
                    cluster = data.markers$cluster,stringsAsFactors = F) 

# 合并所有结果并保存
all_genes_expression <- do.call(rbind, results)
write.csv(all_genes_expression, glue("all_genes_expression_{timestamp}.csv"), row.names = FALSE)
write.csv(marker_data, glue("marker_data_{timestamp}.csv"), row.names = FALSE)
write.csv(all_marker_data, glue("all_marker_{timestamp}.csv"), row.names = FALSE)

print("计算完成")
# 清理所有对象
rm(list = ls())
# 回收内存
gc()
