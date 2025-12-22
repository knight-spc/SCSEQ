library(SingleR)
library(Seurat)

library(RColorBrewer) # 1.1-3
library(dplyr) # 1.0.10
# library(harmony) # 1.2.0
library(pals) # 1.7
library(ggplot2) # 3.3.6
library(glue)
library(readxl)
print('111111111')

library(optparse)

# 定义命令行参数
option_list <- list(
  make_option(c("-d", "--wd"), type="character", default=NULL, 
              help="请输入工作路径"),
  make_option(c("-m", "--model"), type="character", default=NULL, 
              help="请输入注释模型"),
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


singleR_annotation <- function(data_path, singleR_model, timestamp) {
    ###SingleR###
    setwd(data_path)
    input="seurat_data.rds"
    outdir="./"

    data=readRDS(input)
    # SingleR加载数据集
    print("开始加载SingleR数据集")
    print(singleR_model)
    if (singleR_model =="HumanPrimaryCellAtlasData"){
        load("HumanPrimaryCellAtlas_hpca.se_human.RData")
        print("当前加载模型：HumanPrimaryCellAtlas_hpca.se_human.RData")
        singleR_model <- hpca.se
    } else if (singleR_model =="BlueprintEncodeData"){
        load("BlueprintEncode_bpe.se_human.RData")
        singleR_model <- bpe.se
    } else{
        load(glue("{singleR_model}.Rdata"))
        singleR_model <- ref
    }
    print(singleR_model)

    sce <- data ###赋值给其他变量，避免修改原变量。
    # sce_for_SingleR <- GetAssayData(sce, slot="data.NS_1")
    # sce_for_SingleR

    # 转换数据格式
    expression_matrix <- GetAssayData(sce, layer  = "scale.data", assay = "RNA")
    sce_for_SingleR <- as(expression_matrix, "dgCMatrix")
    clusters=sce@meta.data$seurat_clusters
    # 预测细胞类型
    pred.hesc <- SingleR(test = sce_for_SingleR, ref = singleR_model, labels = singleR_model$label.fine,
                     clusters = clusters,  
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
    table(pred.hesc$labels)


        # 创建celltype数据框并确保seurat_clusters为字符型
    celltype = data.frame(
        seurat_clusters = rownames(pred.hesc), 
        celltype = pred.hesc$labels, 
        stringsAsFactors = FALSE
    ) 
    print(celltype)
if ("celltype" %in% colnames(sce@meta.data)) {
    # 如果存在，先删除旧的celltype列
    sce@meta.data$celltype <- NULL
}
print("1111111")
print(str(sce@meta.data$seurat_clusters))
    sce@meta.data <- left_join(sce@meta.data, celltype, by = "seurat_clusters")
print("2222222")

    # 确保celltype是一个因子或字符向量
    # sce@meta.data$celltype <- as.character(sce@meta.data$celltype)
    print(str(sce@meta.data$celltype))
    print(str(sce@meta.data))
print("333333")
    
   


    print(head(sce@meta.data))
    # 保存结果，拿去画图
    setwd(data_path)

    umap_data = data.frame(seurat_clusters = sce@meta.data$seurat_clusters, orig_ident = sce@meta.data$orig.ident, umap_data = sce@reductions$umap@cell.embeddings,annotation = sce@meta.data$celltype,stringsAsFactors = F) 
    umap_data_path = glue("umap_data_{timestamp}.csv")
    write.csv(umap_data, umap_data_path, row.names = FALSE)

    Idents(sce) <- sce@meta.data$celltype

    rds_file_path <- "seurat_data.rds"
    saveRDS(sce, file = rds_file_path)
    print("保存完毕")
    rm(list = ls())
    gc()

    # return(sce) # 返回一个数据框，注意设置stringsAsFactors = FALSE以避免字符串被转换为因子

}

data_path = opt$wd
singleR_model = opt$model
timestamp = opt$time

singleR_annotation(data_path, singleR_model, timestamp)