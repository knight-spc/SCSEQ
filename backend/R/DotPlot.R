library(dplyr) # 1.0.10

library(Seurat)
library(optparse)

print('111111111')

# 定义命令行参数
option_list <- list(
  make_option(c("-w", "--wd"), type="character", default=NULL, 
              help="输入的工作路径"),
  make_option(c("-m", "--marker"), type="character", default=NULL, 
              help="输入的all_marker路径"),
  make_option(c("-a", "--anno"), type="character", default=NULL, 
              help="输入的annotation路径"),
  make_option(c("-n", "--num"), type="character", default=NULL, 
              help="输入的gene_num"),
  make_option(c("-t", "--time"), type="character", default=NULL, 
              help="请输入时间戳")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

getDotPlot <- function(work_path,all_marker_path,annotation_path,gene_num, timestamp ) {
    setwd(work_path)
    print('开始读取数据')
    # 读取 RDS 文件
    seurat_data <- readRDS("seurat_data.rds")
    expression_data <- seurat_data@assays$RNA$data

    data.markers = read.csv(all_marker_path)
    annotation_data = read.csv(annotation_path)
    # 获取基因列表
    top_markers <-data.markers %>%
        group_by(cluster) %>%
        top_n(n = 5, wt = avg_log2FC)
    gene_list <- top_markers$gene 

    results <- list() # 初始化一个空列表来存储结果    
    annotations <- annotation_data$annotation # 提取Annotation列

    print('开始计算基因表达量')
    total_genes <- length(gene_list)
    total_annotations <- length(unique(annotations))

    # 循环遍历每个基因
    for (i in seq_along(gene_list)) {
        gene_of_interest <- gene_list[i]
        # 提取该基因的表达值
        gene_expression <- expression_data[gene_of_interest, ]
        # cat(sprintf("Processing gene %d/%d: %s\n", i, total_genes, gene_of_interest))

        # 循环遍历每个注释类别并计算平均表达量和表达百分比
        for (j in seq_along(unique(annotations))) {
            annot <- unique(annotations)[j]
            # cat(sprintf("  Processing annotation %d/%d: %s\n", j, total_annotations, annot))
            subset_cells <- which(annotations == annot)
            mean_expr <- mean(gene_expression[subset_cells], na.rm = TRUE)
            expr_threshold <- 0.1  # 计算表达该基因的细胞所占的百分比（假设表达量大于某个阈值，如0.1，为表达）
            expr_percentage <- mean(gene_expression[subset_cells] > expr_threshold, na.rm = TRUE)
            
            # 将结果添加到列表中
            result_df <- data.frame(
            Annotation = annot,
            Gene = gene_of_interest,
            Mean_Expression = mean_expr,
            Expression_Percentage = expr_percentage,
            stringsAsFactors = FALSE # 避免将字符向量转换为因子
            )
            results <- c(results, list(result_df))
        }
    }

    # 将结果转换为数据框并保存到文件
    dotplot_data <- do.call(rbind, results)
    print("计算完成，正在保存数据...")
    
    # 生成带时间戳的文件名
    output_file <- paste0("dotplot_data_", timestamp, ".csv")
    write.csv(dotplot_data, output_file, row.names = FALSE)
    
    print("数据已保存")
    rm()
    gc()
}

timestamp = opt$time
work_path = opt$wd
all_marker_path = opt$marker
annotation_path = opt$anno
gene_num = opt$num

getDotPlot(work_path,all_marker_path,annotation_path,gene_num, timestamp )
