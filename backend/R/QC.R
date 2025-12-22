print('111111111')
library(Seurat) # 4.2.0
library(RColorBrewer) # 1.1-3
library(dplyr) # 1.0.10
library(harmony) # 1.2.0
library(pals) # 1.7
library(ggplot2) # 3.3.6
library(glue)
library(readxl)
print('111111111')

library(optparse)

# 定义命令行参数
option_list <- list(
  make_option(c("-s", "--sp"), type="character", default=NULL, 
              help="请输入数据路径"),
  make_option(c("-f", "--fp"), type="character", default=NULL, 
              help="请输入图像存储路径"),
  make_option(c("-p", "--s"), type="character", default=NULL, 
              help="请输入物种信息")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 检查必要参数
if (is.null(opt$sp)) {
  stop("请提供数据路径 (-i/--sp)")
}


# 绘制QC图
get_QC_figure_original <- function(data_path, figure_path, species) {
    setwd(data_path)

    folders <- paste(list.files(), "filtered_feature_bc_matrix", sep = "/")
    print(folders)
    skin <- lapply(folders, function(folder) {
    CreateSeuratObject(
        counts = Read10X(folder),
        project = unlist(strsplit(folder, "/"))[1]
    )
    })
    gc()
    names(skin) <- list.files()

    # 计算MT和HB的百分比
    if(species=="human") {
        HBgenes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")

        for(i in 1:length(skin)){
            skin[[i]][["percent.MT"]] <- PercentageFeatureSet(skin[[i]], pattern = "^MT-")
            HBgenes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
            HB_m <- match(HBgenes, rownames(skin[[i]]@assays$RNA)) 
            HBgenes <- rownames(skin[[i]]@assays$RNA)[HB_m] 
            HBgenes <- HBgenes[!is.na(HBgenes)] 
            skin[[i]][["percent.HB"]]<-PercentageFeatureSet(skin[[i]], features=HBgenes) 
        }

    } else if(species=="mouse") {
        for (i in 1:length(skin)) {
            skin[[i]][["percent.MT"]] <- PercentageFeatureSet(skin[[i]], pattern = "^mt-")
            skin[[i]] [["percent.HB"]] <- PercentageFeatureSet(skin[[i]], pattern = "^Hb")
        }
    }

    setwd(figure_path)
    pdf(paste('01.QC.pdf',sep = '/'))

    for(i in 1:length(skin)){
        p <- VlnPlot(object = skin[[i]], features = c("nCount_RNA", "nFeature_RNA", "percent.MT","percent.HB"),ncol = 2,pt.size = 0)
        print(p)
    }
    dev.off()

    ####merge####

    # setwd("F:/6.Hepatoblastoma/")
    setwd(data_path)

    skinlist <- skin 

    skin <- skinlist[[1]]

    skin <- merge(skin, y = skinlist[-1],  
                add.cell.ids = list.files()
                ,project = "skin")

    VlnPlot(object = skin, features = c("nCount_RNA", "nFeature_RNA", "percent.MT","percent.HB"),ncol = 2,pt.size = 0,raster = F)

    skin@meta.data$Tumor <- substr(skin@meta.data$orig.ident,3,3)
    table(skin$Tumor)

    skin@meta.data$Cancer <- rep("skin",dim(skin@meta.data)[1])
    table(skin$Cancer)

    table(skin$orig.ident)


    setwd(figure_path)

    pdf(paste('02.QC.pdf',sep = '/'))

    VlnPlot(object = skin, features = c("nCount_RNA", "nFeature_RNA", "percent.MT","percent.HB"),ncol = 2,pt.size = 0.1,raster = F)
    VlnPlot(object = skin, features = c("nCount_RNA", "nFeature_RNA", "percent.MT","percent.HB"),group.by = "Tumor"
            ,ncol = 2,pt.size = 0.1,raster = F)
    VlnPlot(object = skin, features = c("nCount_RNA", "nFeature_RNA", "percent.MT","percent.HB"),group.by = "Cancer"
            ,ncol = 2,pt.size = 0.1,raster = F)
    dev.off()
    print("创建完毕")


    rm(list = ls())
    gc()
    return 
    # return(skin)
}

get_QC_figure_final <- function(data_path, figure_path, species, num_nCount_RNA, num_nFeature_RNA, num_percent_MT, num_percent_HB) {
    setwd(data_path)

    folders <- paste(list.files(), "filtered_feature_bc_matrix", sep = "/")
    print(folders)
    skin <- lapply(folders, function(folder) {
    CreateSeuratObject(
        counts = Read10X(folder),
        project = unlist(strsplit(folder, "/"))[1]
    )
    })
    gc()
    names(skin) <- list.files()

    # 计算MT和HB的百分比
    if(species=="human") {
        HBgenes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")

        for(i in 1:length(skin)){
            skin[[i]][["percent.MT"]] <- PercentageFeatureSet(skin[[i]], pattern = "^MT-")
            HBgenes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
            HB_m <- match(HBgenes, rownames(skin[[i]]@assays$RNA)) 
            HBgenes <- rownames(skin[[i]]@assays$RNA)[HB_m] 
            HBgenes <- HBgenes[!is.na(HBgenes)] 
            skin[[i]][["percent.HB"]]<-PercentageFeatureSet(skin[[i]], features=HBgenes) 
        }

    } else if(species=="mouse") {
        for (i in 1:length(skin)) {
            skin[[i]][["percent.MT"]] <- PercentageFeatureSet(skin[[i]], pattern = "^mt-")
            skin[[i]] [["percent.HB"]] <- PercentageFeatureSet(skin[[i]], pattern = "^Hb")
        }
    }

    setwd(figure_path)
    # pdf(paste('01.QC.pdf',sep = '/'))

    # for(i in 1:length(skin)){
    #     p <- VlnPlot(object = skin[[i]], features = c("nCount_RNA", "nFeature_RNA", "percent.MT","percent.HB"),ncol = 2,pt.size = 0)
    #     print(p)
    # }
    # dev.off()

    ####merge####

    # setwd("F:/6.Hepatoblastoma/")
    setwd(data_path)

    skinlist <- skin 

    skin <- skinlist[[1]]

    skin <- merge(skin, y = skinlist[-1],  
                add.cell.ids = list.files()
                ,project = "skin")

    # VlnPlot(object = skin, features = c("nCount_RNA", "nFeature_RNA", "percent.MT","percent.HB"),ncol = 2,pt.size = 0,raster = F)

    skin@meta.data$Tumor <- substr(skin@meta.data$orig.ident,3,3)
    table(skin$Tumor)

    skin@meta.data$Cancer <- rep("skin",dim(skin@meta.data)[1])
    table(skin$Cancer)

    table(skin$orig.ident)


    setwd(figure_path)

    # pdf(paste('02.QC.pdf',sep = '/'))

    # VlnPlot(object = skin, features = c("nCount_RNA", "nFeature_RNA", "percent.MT","percent.HB"),ncol = 2,pt.size = 0,raster = F)
    # VlnPlot(object = skin, features = c("nCount_RNA", "nFeature_RNA", "percent.MT","percent.HB"),group.by = "Tumor"
    #         ,ncol = 2,pt.size = 0,raster = F)
    # VlnPlot(object = skin, features = c("nCount_RNA", "nFeature_RNA", "percent.MT","percent.HB"),group.by = "Cancer"
    #         ,ncol = 2,pt.size = 0,raster = F)

    # dev.off()

    #metadata

    # skin <- subset(skin, subset = nCount_RNA >2000 & nFeature_RNA > 200 & percent.MT < 25 &percent.HB < 5)
    skin <- subset(skin, subset = nCount_RNA > num_nCount_RNA & nFeature_RNA > num_nCount_RNA & percent.MT < num_percent_MT & percent.HB< num_percent_HB)

    #QC plot 3

    setwd(figure_path)

    pdf(opt$pathpdf, width = 6, height = 6)

    VlnPlot(object = skin, features = c("nCount_RNA", "nFeature_RNA", "percent.MT","percent.HB"),ncol = 2,pt.size = 0.1,raster = F)
    VlnPlot(object = skin, features = c("nCount_RNA", "nFeature_RNA", "percent.MT","percent.HB"),group.by = "Tumor"
            ,ncol = 2,pt.size = 0.1,raster = F)
    VlnPlot(object = skin, features = c("nCount_RNA", "nFeature_RNA", "percent.MT","percent.HB"),group.by = "Cancer"
            ,ncol = 2,pt.size = 0.1,raster = F)

    dev.off()

    options(future.globals.maxSize = 60000 * 1024^2)  # 增加内存限制

    print('222222222')
    skin <- NormalizeData(skin, normalization.method = "LogNormalize", scale.factor = 10000)
    print('333333333')
    skin <- FindVariableFeatures(skin, selection.method = "vst", nfeatures = 5000)
    print('4444444444')
    skin <- ScaleData(skin, features = rownames(skin), assay = "RNA")
    print('5555555555')
    skin <- RunPCA(skin)
    print('6666666666')
    skin <- RunHarmony(skin, group.by.vars = "orig.ident", dims.use = 1:20, lambda = 1, project.dim = F)
    print('7777777777')
    skin <- RunUMAP(skin, dim = 1:20, reduction = "harmony")
    print('8888888888')
    skin <- FindNeighbors(skin, dim = 1:20, reduction = "harmony")
    print('9999999999')
    skin <- FindClusters(skin, resolution = 0.1)
    print('00000000000')
    skin <- RunTSNE(skin, dims = 1:20, do.fast = TRUE, reduction = "harmony")
    print(skin)
    print("02.umap.pdf")

    pdf("02.umap.pdf", width = 6, height = 6)
    p <- DimPlot(skin, group.by = "seurat_clusters", reduction = "umap", label = TRUE, cols = c(cols25(), cols25())) + guides(colour = guide_legend(override.aes = list(size = 3), nrow = 20))
    print(p)
    p <- DimPlot(skin, reduction = "umap", group.by = "orig.ident", label = TRUE, cols = c(cols25(), cols25())) + guides(colour = guide_legend(override.aes = list(size = 3), nrow = 20))
    print(p)
    dev.off()

    print("02.umap.pdf绘制完成，路径是：")
    print(figure_path)
    sce <- skin ###赋值给其他变量，避免修改原变量。
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
    # write.table(skin, file.path(figure_path, "skin.txt") , sep = "\t", row.names = FALSE, quote = FALSE)
    print("保存完毕")

    rm(list = ls())
    gc()



}

data_path = opt$sp
figure_path = opt$fp
species = opt$s
get_QC_figure_original(data_path, figure_path, species)