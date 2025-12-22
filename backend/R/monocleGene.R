library(monocle)
library(RColorBrewer)
library(dplyr)
library(optparse)
options(bitmapType='cairo')
# 定义命令行参数
option_list <- list(

  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="输入的 Seurat 对象路径"),
  make_option(c("-p", "--pdfpath"), type="character", default=NULL, 
              help="输入的 Seurat 对象路径"),
  make_option(c("-s", "--SubsetList"), type="character", default=NULL, 
              help="选择要观察的细胞群体")
)
# 解析命令行参数
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
pdf_path = opt$pdfpath
gene = opt$SubsetList
print('正在读取数据')
HSMM_restored <- readRDS(opt$input)
gene_list <- as.character(unlist(strsplit(gene,split=",")))
print("输入的基因列表：")
print(gene_list)
# 检查HSMM对象中的基因名称
all_genes <- fData(HSMM_restored)$gene_short_name
# 检查HSMM对象中所有可用的基因名称
print("HSMM对象中所有基因的前几个示例：")
print(head(all_genes, 20))  # 显示前20个基因名称
print("在HSMM对象中找到的匹配基因：")
print(gene_list[gene_list %in% all_genes])
heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'RdBu')))
example_genes <- row.names(subset(fData(HSMM_restored),
                                gene_short_name %in% gene_list))
print("最终用于绘图的基因行名：")
print(example_genes)

cds_subset <- HSMM_restored[example_genes,]
print("子集维度：")
print(dim(cds_subset))

# # 新增：去除全零基因和含NA的基因
# exprs_mat <- exprs(cds_subset)
# keep_genes <- rowSums(exprs_mat, na.rm=TRUE) > 0 & !apply(exprs_mat, 1, function(x) any(is.na(x) | is.infinite(x)))
# cds_subset <- cds_subset[keep_genes,]
# print("过滤后子集维度：")
# print(dim(cds_subset))

print(pdf_path)
pdf(paste(pdf_path,sep = '/'))

# 计算需要多少页
filtered_genes <- rownames(cds_subset)
n_genes <- length(filtered_genes)
n_pages <- ceiling(n_genes/6)

# 按每页6个基因分组绘图，单基因try-catch
for(i in 1:n_pages) {
    start_idx <- (i-1)*6 + 1
    end_idx <- min(i*6, n_genes)
    current_genes <- filtered_genes[start_idx:end_idx]
    cds_page <- cds_subset[current_genes,]
    tryCatch({
        p <- plot_genes_in_pseudotime(cds_page, 
                                color_by = "Pseudotime",
                                ncol = 1,
                                relative_expr = FALSE) +
            scale_color_gradientn(colours = heatmap_Palette(100))
        print(p)
    }, error=function(e){
        cat("绘制基因", paste(current_genes, collapse=","), "时出错，已跳过。\n")
    })
}
dev.off()
# 清理内存
rm(cds_subset)
rm(HSMM_restored)
gc()