library(infercnv)
library(AnnoProbe)
library(Seurat)
library(optparse)
options(bitmapType='cairo')
# 定义命令行参数
option_list <- list(
  make_option(c("-w", "--wd"), type="character", default=NULL, 
              help="输入的工作路径"),
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="输入的 Seurat 对象路径"),
  make_option(c("-s", "--SubsetList"), type="character", default=NULL, 
              help="选择要观察的细胞群体"),
  make_option(c("-r", "--ReferenceList"), type="character", default=NULL,
              help="选择作为参照的细胞群体")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 检查必要参数
if (is.null(opt$input)) {
  stop("请提供输入文件路径 (-i/--input)")
}

# 读取数据
setwd(opt$wd)
seurat_object <- readRDS(opt$input)
if (!is.null(opt$SubsetList)) {
  
    # 将输入的字符串按逗号分割并去除空格
    cell_types <- strsplit(opt$SubsetList, ",")[[1]]
    print(cell_types)

    cell_types <- trimws(cell_types)  # 去除空格
    print(cell_types)

    Idents(seurat_object) <- seurat_object@meta.data$celltype
    print("Available cell types in object:")
    print(unique(Idents(seurat_object)))
    seurat_object <- subset(seurat_object, idents=cell_types)
    print("daozhemeiwenti")

}

# 表达矩阵文件
counts <- GetAssayData(seurat_object, layer = 'counts')
# 移除NA和Inf
counts[is.na(counts)] <- 0
counts[is.infinite(counts)] <- 0
# 只保留表达量大于0的基因和细胞
counts <- counts[rowSums(counts) > 0, , drop=FALSE]
counts <- counts[, colSums(counts) > 0, drop=FALSE]




# 检查维度
if(nrow(counts) < 10 || ncol(counts) < 10) {
    stop("过滤后基因或细胞数量过少，无法进行PCA。")
}

# 保证anno和counts列名一致
anno <- data.frame(Idents(seurat_object))
anno <- anno[colnames(counts), , drop=FALSE]

# 基因注释文件
gene_order <- "pth/human_gene_pos.txt"
# 新增数据预处理步骤
# counts <- counts[rowSums(counts) > 0, ]  # 过滤全零行
# counts <- log2(counts + 1)               # 添加log标准化
# counts <- scale(counts)                  # 数据标准化
# anno <- anno[colnames(counts), , drop = FALSE]



    # 将输入的字符串按逗号分割并去除空格
    referenceList_types <- strsplit(opt$ReferenceList, ",")[[1]]
    referenceList_types <- trimws(referenceList_types)  # 去除空格
    
    # # 过滤低表达基因（添加额外的有效性检查）
    # counts[is.na(counts)] <- 0  # 处理缺失值
    # genes_to_keep <- rowSums(counts > 0) >= 10  # 至少在10个细胞中表达
    # counts <- counts[genes_to_keep, ]
    # counts <- counts[Matrix::rowSums(counts) > 0, ]  # 确保没有全零行
    # counts <- counts[, Matrix::colSums(counts) > 0]  # 确保没有全零列
    
    # # 新增调试信息输出
    # print(paste("过滤后基因数:", nrow(counts)))
    # print(paste("过滤后细胞数:", ncol(counts)))
    
    # # 新增二次过滤确保有效性
    # if(nrow(counts) == 0 || ncol(counts) == 0) {
    #     stop("数据过滤后矩阵为空，请检查输入参数或数据质量")
    # }
    
    

    
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts,
                                    annotations_file = anno,
                                    delim="\t",
                                    gene_order_file = gene_order,
                                    ref_group_names = referenceList_types)
print("创建完成")
infercnv_obj = infercnv::run(infercnv_obj,
                            cutoff = 0.1,
                            out_dir = ".", 
                            cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE,

                            )
print("infercnv完成")
rm(list = ls())
gc()  # 强制垃圾回收
