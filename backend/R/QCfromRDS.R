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

# 定义命令行参数
option_list <- list(
  make_option(c("-w", "--wd"), type="character", default=NULL, 
            help="输入的工作路径"),

  make_option(c("-s", "--s"), type="character", default=NULL, 
            help="请输入物种信息")

)

# 解析命令行参数
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 检查必要参数
if (is.null(opt$wd)) {
  stop("请提供数据路径 (-i/--wd)")
}

get_QC_figure_final <- function(data_path, species, num_nCount_RNA = 2000, num_nFeature_RNA = 200, num_percent_MT = 25, num_percent_HB = 5) {
    setwd(data_path)
    data <- readRDS("HCC_10_filter_v3.rds")

    # 检查数据结构
    print("检查数据结构：")
    print(dim(data))
    print("细胞数量：")
    print(ncol(data))



    print('111111111')



    # data <- subset(data, subset = nCount_RNA >2000 & nFeature_RNA > 200 & percent.MT < 25 &percent.HB < 5)
    # data <- subset(data, subset = nCount_RNA > num_nCount_RNA & nFeature_RNA > num_nCount_RNA)


    print('222222222')
    data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
    print('333333333')
    data <- FindVariableFeatures(data)
    print('4444444444')
    data <- ScaleData(data, features = rownames(data), assay = "RNA")
    print('5555555555')

    # 运行PCA
    print(data)
    data <- RunPCA(data,npcs = 30,fastpath=FALSE)
    print('6666666666')
    data <- RunHarmony(data, group.by.vars = "orig.ident", dims.use = 1:20, lambda = 1, project.dim = F)
    print('7777777777')
    data <- RunUMAP(data, dim = 1:20, reduction = "harmony")
    print('8888888888')
    data <- FindNeighbors(data, dim = 1:20, reduction = "harmony")
    print('9999999999')
    data <- FindClusters(data, resolution = 0.1)
    print('00000000000')
    data <- RunTSNE(data, dims = 1:20, do.fast = TRUE, reduction = "harmony")
    print(data)
    print("02.umap.pdf")

    pdf("02.umap.pdf", width = 6, height = 6)
    p <- DimPlot(data, group.by = "seurat_clusters", reduction = "umap", label = TRUE, cols = c(cols25(), cols25())) + guides(colour = guide_legend(override.aes = list(size = 3), nrow = 20))
    print(p)
    p <- DimPlot(data, reduction = "umap", group.by = "orig.ident", label = TRUE, cols = c(cols25(), cols25())) + guides(colour = guide_legend(override.aes = list(size = 3), nrow = 20))
    print(p)
    dev.off()

    print("02.umap.pdf绘制完成，路径是：")
    sce <- data ###赋值给其他变量，避免修改原变量。
    umap_data = data.frame(seurat_clusters = sce@meta.data$seurat_clusters, orig_ident = sce@meta.data$orig.ident, umap_data = sce@reductions$umap@cell.embeddings,stringsAsFactors = F) 
    write.csv(umap_data, "umap_data.csv", row.names = FALSE)

    # 计算并保存marker基因数据
    sce.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    sce.markers.top20 = sce.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
    marker_data = data.frame(gene = sce.markers.top20$gene, 
                            p_val = sce.markers.top20$p_val, 
                            avg_log2FC = sce.markers.top20$avg_log2FC,
                            pct.1 = sce.markers.top20$pct.1,
                            pct.2 = sce.markers.top20$pct.2,
                            p_val_adj = sce.markers.top20$p_val_adj,
                            cluster = sce.markers.top20$cluster,stringsAsFactors = F) 
    write.csv(marker_data, "marker_data.csv", row.names = FALSE)
    # 计算并保存allmarker基因数据
    all_marker_data = data.frame(gene = sce.markers$gene, 
                            p_val = sce.markers$p_val, 
                            avg_log2FC = sce.markers$avg_log2FC,
                            pct.1 = sce.markers$pct.1,
                            pct.2 = sce.markers$pct.2,
                            p_val_adj = sce.markers$p_val_adj,
                            cluster = sce.markers$cluster,stringsAsFactors = F) 
    write.csv(all_marker_data, "all_marker.csv", row.names = FALSE)


    rds_file_path <- "seurat_data.rds"
    saveRDS(sce, file = rds_file_path)
    # write.table(data, file.path(figure_path, "data.txt") , sep = "\t", row.names = FALSE, quote = FALSE)
    print("保存完毕")

    rm(list = ls())
    gc()



}


wd_path = opt$wd
species = opt$s

get_QC_figure_final(wd_path, species)
