library(Seurat)
library(optparse)

# 命令行参数
option_list = list(
    make_option(c("-w", "--work_path"), type="character", help="工作目录路径"),
    make_option(c("-g", "--gene"), type="character", help="基因名称"),
    make_option(c("-m", "--group_model"), type="character", default="celltype", help="分组模型"),
    make_option(c("-t", "--timestamp"), type="character", help="时间戳")
)
opt = parse_args(OptionParser(option_list=option_list))

# 读取数据
data <- readRDS(file.path(opt$work_path, "seurat_data.rds"))

# 生成PDF图
pdf_path <- file.path(opt$work_path, paste0("marker_genes_", opt$timestamp, ".pdf"))
pdf(pdf_path, width = 6, height = 6)
p <- VlnPlot(data, features = opt$gene, pt.size = 0.1, group.by = opt$group_model)
print(p)
dev.off()

# 获取基因表达数据
expression_data <- data@assays$RNA$data
gene_index <- which(rownames(expression_data) == opt$gene)
if (length(gene_index) != 1) {
    stop("Gene of interest not found in the expression data or multiple matches found.")
}

gene_expression <- expression_data[gene_index, ]
expression_path <- file.path(opt$work_path, paste0("gene_expression_", opt$timestamp, ".csv"))
write.csv(as.data.frame(gene_expression), expression_path)

# 清理内存
rm(list = ls())

gc()