
####Working directory####

setwd("F:/6.Hepatoblastoma/")

####Package####

library(Seurat) #4.2.0
library(RColorBrewer) #1.1-3
library(dplyr) #1.0.10
library(harmony) #1.2.0
library(pals) #1.7
library(ggplot2) #3.3.6
library(readxl)

####Load data####

#Clinical information

basicinfo <- read_excel("info/Clinical_information_1to28.xlsx", sheet = "Sheet1") 
basicinfo <- basicinfo[,-c(1:4)]

#Single cell

list.files("sc/self_HB/")

Filepath <- paste("sc/self_HB",list.files("sc/self_HB/"),'filtered_feature_bc_matrix',sep = '/')

HB <- list()

for (i in 1:length(Filepath)) {
  
seurat_object <- CreateSeuratObject(counts = Read10X(Filepath[i]), 
                                    project = list.files("sc/self_HB/")[i] )
HB[[i]] <- seurat_object

}

rm(seurat_object)
gc()

names(HB) = list.files("sc/self_HB/")

####QC####
for(i in 1:length(HB)){
  HB[[i]] [["percent.MT"]] <- PercentageFeatureSet(HB[[i]], pattern = "^MT-")
}

HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")

for(i in 1:length(HB)){
  HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  HB_m <- match(HB.genes, rownames(HB[[i]]@assays$RNA)) 
  HB.genes <- rownames(HB[[i]]@assays$RNA)[HB_m] 
  HB.genes <- HB.genes[!is.na(HB.genes)] 
  HB[[i]][["percent.HB"]]<-PercentageFeatureSet(HB[[i]], features=HB.genes) 
}


#QC plot 1

setwd("F:/6.Hepatoblastoma/fig")

pdf(paste('01.QC.pdf',sep = '/'))

for(i in 1:length(HB)){
  p <- VlnPlot(object = HB[[i]], features = c("nCount_RNA", "nFeature_RNA", "percent.MT","percent.HB"),ncol = 2,pt.size = 0)
  print(p)
}
dev.off()

####merge####

setwd("F:/6.Hepatoblastoma/")

HBlist <- HB 

HB <- HBlist[[1]]

HB <- merge(HB, y = HBlist[-1],  
            add.cell.ids = list.files("sc/self_HB/")
            ,project = "HB")

VlnPlot(object = HB, features = c("nCount_RNA", "nFeature_RNA", "percent.MT","percent.HB"),ncol = 2,pt.size = 0,raster = F)

HB@meta.data$Tumor <- substr(HB@meta.data$orig.ident,3,3)
table(HB$Tumor)

HB@meta.data$Cancer <- rep("HB",dim(HB@meta.data)[1])
table(HB$Cancer)

table(HB$orig.ident)

#QC plot 2

setwd("F:/6.Hepatoblastoma/fig")

pdf(paste('02.QC.pdf',sep = '/'))

VlnPlot(object = HB, features = c("nCount_RNA", "nFeature_RNA", "percent.MT","percent.HB"),ncol = 2,pt.size = 0,raster = F)
VlnPlot(object = HB, features = c("nCount_RNA", "nFeature_RNA", "percent.MT","percent.HB"),group.by = "Tumor"
        ,ncol = 2,pt.size = 0,raster = F)
VlnPlot(object = HB, features = c("nCount_RNA", "nFeature_RNA", "percent.MT","percent.HB"),group.by = "Cancer"
        ,ncol = 2,pt.size = 0,raster = F)

dev.off()

#metadata

HB <- subset(HB, subset = nCount_RNA >2000 & nFeature_RNA > 200 & percent.MT < 25 &percent.HB < 5)

#QC plot 3

setwd("F:/6.Hepatoblastoma/fig")

pdf(paste('03.QC.pdf',sep = '/'))

VlnPlot(object = HB, features = c("nCount_RNA", "nFeature_RNA", "percent.MT","percent.HB"),ncol = 2,pt.size = 0,raster = F)
VlnPlot(object = HB, features = c("nCount_RNA", "nFeature_RNA", "percent.MT","percent.HB"),group.by = "Tumor"
        ,ncol = 2,pt.size = 0,raster = F)
VlnPlot(object = HB, features = c("nCount_RNA", "nFeature_RNA", "percent.MT","percent.HB"),group.by = "Cancer"
        ,ncol = 2,pt.size = 0,raster = F)

dev.off()

####Clustering####

HB <- NormalizeData(HB, normalization.method = "LogNormalize", scale.factor = 10000)
HB <- FindVariableFeatures(HB, selection.method = "vst", nfeatures = 5000)
HB <- ScaleData(HB,features = rownames(HB), assay = "RNA")
HB <- RunPCA(HB)
HB <- RunHarmony(HB, group.by.vars = "orig.ident", dims.use = 1:20, lambda = 1,project.dim = F)
HB <- RunUMAP(HB, dim = 1:20, reduction = "harmony")
HB <- FindNeighbors(HB, dim = 1:20, reduction = "harmony")
HB <- FindClusters(HB, resolution = 1)
HB <- RunTSNE(HB, dims = 1:20,do.fast = TRUE, reduction = "harmony")


table(HB$seurat_clusters)

HB.markers <- FindAllMarkers(HB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
HB.markers.top01 = HB.markers %>% group_by(cluster) %>% top_n(n = 01, wt = avg_log2FC)
HB.markers.top02 = HB.markers %>% group_by(cluster) %>% top_n(n = 02, wt = avg_log2FC)
HB.markers.top10 = HB.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
HB.markers.top20 = HB.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

DimPlot(HB, reduction = "umap", label = TRUE,raster=FALSE)
DimPlot(HB, reduction = "tsne", label = TRUE,raster=FALSE)


















