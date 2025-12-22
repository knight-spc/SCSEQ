


HCC_kuffer <- NormalizeData(HCC_kuffer, normalization.method = "LogNormalize", scale.factor = 10000)
HCC_kuffer <- FindVariableFeatures(HCC_kuffer, selection.method = "vst", nfeatures = 5000)
HCC_kuffer <- ScaleData(HCC_kuffer,features = rownames(HCC_kuffer), assay = "RNA")
HCC_kuffer <- RunPCA(HCC_kuffer)
HCC_kuffer <- RunHarmony(HCC_kuffer, group.by.vars = "Patient_Part", dims.use = 1:20, lambda = 1,project.dim = F)
HCC_kuffer <- RunUMAP(HCC_kuffer, dim = 1:20, reduction = "harmony")
HCC_kuffer <- FindNeighbors(HCC_kuffer, dim = 1:20, reduction = "harmony")
HCC_kuffer <- FindClusters(HCC_kuffer, resolution = 0.5)
HCC_kuffer <- RunTSNE(HCC_kuffer, dims = 1:20,do.fast = TRUE, reduction = "harmony")


table(HCC_kuffer$seurat_clusters)

####marker volcano####

####load####
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

####access####

scRNA = HCC_kuffer
scRNA@meta.data$seurat_clusters = Idents(scRNA)
ctys = levels(scRNA)
ctys

VlnPlot(HCC_kuffer,group.by = "seurat_clusters",features=c('PTPRC'),ncol=1,pt.size = 0)
VlnPlot(HCC_kuffer,group.by = "seurat_clusters",features=c('APOE',"C1QC","C1QB","C1QA"),ncol=1,pt.size = 0)

scRNA.markers <- FindAllMarkers(scRNA, min.pct = 0.25, 
                                logfc.threshold = 0.25)
head(scRNA.markers)

colnames(scRNA.markers)[6] = "seurat_clusters"
k = scRNA.markers$p_val_adj<0.05;table(k)

scRNA.markers = scRNA.markers[k,]

#上下调
scRNA.markers$label <- ifelse(scRNA.markers$avg_log2FC<0,"sigDown","sigUp")
topgene <- scRNA.markers %>%
  group_by(seurat_clusters) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  bind_rows(group_by(scRNA.markers, seurat_clusters) %>%
              top_n(n = 10, wt = -avg_log2FC))
head(topgene)

#根据log2FC范围确定背景柱长度：
dfbar = scRNA.markers %>%
  group_by(seurat_clusters) %>%
  summarise(low = round(min(avg_log2FC)-0.5),
            up = round(max(avg_log2FC)+0.5))

#绘制背景柱和散点图：
p1 <- ggplot()+
  geom_col(aes(x = seurat_clusters ,y = low),dfbar,
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(aes(x = seurat_clusters ,y = up),dfbar,
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(aes(x = seurat_clusters, y = avg_log2FC, color = label),scRNA.markers,
              width =0.4,size = 1)+
  scale_color_manual(values = c("#0077c0","#c72d2e"))+
  theme_classic()
p1

#X轴的色块标签：
library(RColorBrewer)
mycol <- colorRampPalette(rev(brewer.pal(n = 7, name ="Set1")))(length(ctys))
p2 <- p1 + 
  geom_tile(aes(x = ctys,y = 0),
            height = 0.5,fill = mycol, show.legend = F)+
  geom_text(aes(x= ctys, y = 0, label = ctys),
            size = 3,fontface = "bold")
p2

library(ggrepel)
#给每种细胞类型的top基因加上标签,调整细节：
p3 <- p2 + 
  geom_text_repel(aes(x = seurat_clusters,y = avg_log2FC,label = gene),
                  topgene,size = 3 )+
  labs(x = "CellType",y = "Average log2FoldChange",
       title = "Differential expression genes")+
  theme(
    plot.title = element_text(size = 14,color = "black",face = "bold"),
    axis.title = element_text(size = 12,color = "black",face = "bold"),
    axis.line.y = element_line(color = "black",linewidth = 0.8),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position.inside  = c(0.98,0.96),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 12)
  )+
  guides(color = guide_legend(override.aes = list(size = 4)))  
p3

####删除了俩样本后，重新聚类####

Idents(HCC_n) <- HCC_n$Patient
HCC_s <- subset(HCC_n,idents = c("P03","P06","P09"))
HCC_s0 <- subset(HCC_n,idents = c("P03","P06","P09"))

HCC_s <- NormalizeData(HCC_s, normalization.method = "LogNormalize", scale.factor = 10000)
HCC_s <- FindVariableFeatures(HCC_s, selection.method = "vst", nfeatures = 5000)
HCC_s <- ScaleData(HCC_s,features = rownames(HCC_s), assay = "RNA")
HCC_s <- RunPCA(HCC_s)
HCC_s <- RunHarmony(HCC_s, group.by.vars = "Patient_Part", dims.use = 1:20, lambda = 1,project.dim = F)
HCC_s <- RunUMAP(HCC_s, dim = 1:20, reduction = "harmony")
HCC_s <- FindNeighbors(HCC_s, dim = 1:20, reduction = "harmony")
HCC_s <- FindClusters(HCC_s, resolution = 0.5)
HCC_s <- RunTSNE(HCC_s, dims = 1:20,do.fast = TRUE, reduction = "harmony")

table(HCC_s$seurat_clusters)

#####singleR####

#ref

#HCC_s

hpca.se <- readRDS("G:/SingleR_Celldex_ref/Human_General_HumanPrimaryCellAtlasData.rds")
hpca.se


sce_for_SingleR <- GetAssayData(HCC_s, layer ="data")
sce_for_SingleR

sce <- HCC_s
clusters=sce@meta.data$seurat_clusters
pred.hesc <- SingleR(test = sce_for_SingleR, ref = hpca.se, labels = hpca.se$label.fine,
                     clusters = clusters,  
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
table(pred.hesc$labels)
celltype_HCC_s = data.frame(ClusterID=rownames(pred.hesc), celltype=pred.hesc$labels, stringsAsFactors = F) 

@meta.data$labels <-pred.hesc$label

gc()

####annotation####

VlnPlot(HCC_s,group.by = "seurat_clusters",features=c('PTPRC'),ncol=1,pt.size = 0)
VlnPlot(HCC_s,group.by = "seurat_clusters",features=c('APOE',"C1QC","C1QB","C1QA"),ncol=1,pt.size = 0) #kuffer 3,8
VlnPlot(HCC_s,group.by = "seurat_clusters",features=c('CST3',"HLA-DRA","CLEC10A"),ncol=1,pt.size = 0) #DC 8
VlnPlot(HCC_s,group.by = "seurat_clusters",features=c('MYL9',"TAGLN","COL1A2","ACTA2"),ncol=1,pt.size = 0) #fibroblast 12
VlnPlot(HCC_s,group.by = "seurat_clusters",features=c('APOB',"FGB","APOA1","APOE"),ncol=1,pt.size = 0) #hepatocyte 3
VlnPlot(HCC_s,group.by = "seurat_clusters",features=c("PECAM1","AQP1","VWF","PLVAP"),ncol=1,pt.size = 0) #endothelial 6
VlnPlot(HCC_s,group.by = "seurat_clusters",features=c("CD3D","CD3E","CD4","CD8A","CD8B"),ncol=1,pt.size = 0) #T cell 0,1,10,11,14
VlnPlot(HCC_s,group.by = "seurat_clusters",features=c("CD79A","IGHM"),ncol=1,pt.size = 0) #B cell 4
VlnPlot(HCC_s,group.by = "seurat_clusters",features=c("CCL4","PRF1","KLRD1"),ncol=1,pt.size = 0) #NK 1,11
VlnPlot(HCC_s,group.by = "seurat_clusters",features=c("LYZ","FPR1","FCN1"),ncol=1,pt.size = 0) #monocyte 5
VlnPlot(HCC_s,group.by = "seurat_clusters",features=c("HDC","KIT","CPA3","TPSAB1"),ncol=1,pt.size = 0) #mast cell
VlnPlot(HCC_s,group.by = "seurat_clusters",features=c("HBA2","STMN1"),ncol=1,pt.size = 0)


table(HCC_s$seurat_clusters)

HCC_s@meta.data$Annotation <- rep(NA, dim(HCC_s@meta.data)[1])

HCC_s$Annotation[HCC_s$seurat_clusters %in% c(3)] <- "Hepatocyte" #'APOB',"FGB","APOA1","APOE",
HCC_s$Annotation[HCC_s$seurat_clusters %in% c(12)] <- "Fibroblast" #'MYL9',"TAGLN","COL1A2","ACTA2",
HCC_s$Annotation[HCC_s$seurat_clusters %in% c(6)] <- "Endothelial" #"PECAM1","AQP1","VWF","PLVAP",

HCC_s$Annotation[!(HCC_s$seurat_clusters %in% c(3,6,12))] <- "Immune" #"PTPRC",

table(HCC_s$Annotation)

HCC_s@meta.data$Annotation1 <- HCC_s@meta.data$Annotation  
HCC_s$Annotation1[(HCC_s$seurat_clusters %in% c(4))] <- "B cell" #"CD79A",
HCC_s$Annotation1[(HCC_s$seurat_clusters %in% c(4))] <- "B cell" #"CD79A",



####immune####

Idents(HCC_s) <- HCC_s$Annotation

HCC_s_im <- subset(HCC_s,idents = c("Immune"))
#####harmony####
HCC_s_im <- NormalizeData(HCC_s_im, normalization.method = "LogNormalize", scale.factor = 10000)
HCC_s_im <- FindVariableFeatures(HCC_s_im, selection.method = "vst", nfeatures = 5000)
HCC_s_im <- ScaleData(HCC_s_im,features = rownames(HCC_s_im), assay = "RNA")
HCC_s_im <- RunPCA(HCC_s_im)
HCC_s_im <- RunHarmony(HCC_s_im, group.by.vars = "Patient_Part", dims.use = 1:20, lambda = 1,project.dim = F)
HCC_s_im <- RunUMAP(HCC_s_im, dim = 1:20, reduction = "harmony")
HCC_s_im <- FindNeighbors(HCC_s_im, dim = 1:20, reduction = "harmony")
HCC_s_im <- FindClusters(HCC_s_im, resolution = 0.2)
HCC_s_im <- RunTSNE(HCC_s_im, dims = 1:20,do.fast = TRUE, reduction = "harmony")

#####cca####
# 使用CCA代替Harmony
HCC_s_im <- NormalizeData(HCC_s_im, normalization.method = "LogNormalize", scale.factor = 10000)
HCC_s_im <- FindVariableFeatures(HCC_s_im, selection.method = "vst", nfeatures = 5000)
HCC_s_im <- ScaleData(HCC_s_im, features = rownames(HCC_s_im), assay = "RNA")
HCC_s_im <- RunPCA(HCC_s_im)

# 使用CCA进行数据整合
HCC_s_im <- RunCCA(HCC_s_im, features = rownames(HCC_s_im))

# 进行UMAP降维
HCC_s_im <- RunUMAP(HCC_s_im, dims = 1:20, reduction = "cca")

# 找邻居和进行聚类
HCC_s_im <- FindNeighbors(HCC_s_im, dims = 1:20, reduction = "cca")
HCC_s_im <- FindClusters(HCC_s_im, resolution = 0.2)

# 使用tSNE进行可视化
HCC_s_im <- RunTSNE(HCC_s_im, dims = 1:20, do.fast = TRUE, reduction = "cca")

#####cca2####

# 1. 按照 Patient 字段拆分对象
HCC_list <- SplitObject(HCC_s_im, split.by = "Patient")

# 2. 对每个子对象进行标准化
HCC_list <- lapply(HCC_list, NormalizeData)

# 3. 找到每个子对象的变量特征
HCC_list <- lapply(HCC_list, FindVariableFeatures, selection.method = "vst", nfeatures = 2000)

# 4. 选择共同特征
features_to_integrate <- SelectIntegrationFeatures(object.list = HCC_list, nfeatures = 2000)

# 5. 找到整合锚点
HCC_anchors <- FindIntegrationAnchors(object.list = HCC_list, anchor.features = features_to_integrate)

# 6. 整合数据
HCC_combined <- IntegrateData(anchorset = HCC_anchors)

# 7. 运行 PCA（对整合后的数据进行主成分分析）
HCC_combined <- RunPCA(HCC_combined)

# 8. 进一步分析与降维
HCC_combined <- RunUMAP(HCC_combined, dims = 1:20)

# 9. 运行 t-SNE
HCC_combined <- RunTSNE(HCC_combined, dims = 1:20, do.fast = TRUE)

# 10. 聚类（设置resolution）
HCC_combined <- FindNeighbors(HCC_combined, dims = 1:20)  # 你可以选择适合的dim范围
HCC_combined <- FindClusters(HCC_combined, resolution = 0.5)  # 设置 resolution，控制聚类结果的细致度
#####res####
table(HCC_s_im$seurat_clusters)
table(HCC_s_im$RNA_snn_res.0.5)
table(HCC_s_im$seurat_clusters,HCC_s_im$RNA_snn_res.0.5)


table(HCC_s$seurat_clusters)

gc()

####annotation####

VlnPlot(HCC_s_im,group.by = "seurat_clusters",features=c('PTPRC'),ncol=1,pt.size = 0)

VlnPlot(HCC_s_im,group.by = "seurat_clusters",features=c("CD3D","CD3E","CD4","CD8A","CD8B"),ncol=1,pt.size = 0) #T cell
VlnPlot(HCC_s_im,group.by = "seurat_clusters",features=c("CCL4","PRF1","KLRD1"),ncol=1,pt.size = 0) #NK 

VlnPlot(HCC_s_im,group.by = "seurat_clusters",features=c("CD3D","CD3E","CD4","CD8A","CD8B",
                                                         "CCL4","PRF1","KLRD1"),ncol=1,pt.size = 0) #T/NK cell


VlnPlot(HCC_s_im,group.by = "seurat_clusters",features=c("CD79A","IGHM"),ncol=1,pt.size = 0) #B cell 

VlnPlot(HCC_s_im,group.by = "seurat_clusters",features=c('APOE',"C1QC","C1QB","C1QA"),ncol=1,pt.size = 0) #kuffer 
VlnPlot(HCC_s_im,group.by = "seurat_clusters",features=c('CST3',"HLA-DRA","CLEC10A"),ncol=1,pt.size = 0) #DC 


VlnPlot(HCC_s_im,group.by = "seurat_clusters",features=c("LYZ","FPR1","FCN1","CD68"),ncol=1,pt.size = 0) #monocyte 
VlnPlot(HCC_s_im,group.by = "seurat_clusters",features=c("HDC","KIT","CPA3","TPSAB1"),ncol=1,pt.size = 0) #mast cell
VlnPlot(HCC_s_im,group.by = "seurat_clusters",features=c("HBA2","STMN1"),ncol=1,pt.size = 0)

####marker####
HCC_s_im.markers <- FindAllMarkers(HCC_s_im, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
HCC_s_im.markers.top01 = HCC_s_im.markers %>% group_by(cluster) %>% top_n(n = 01, wt = avg_log2FC)
HCC_s_im.markers.top02 = HCC_s_im.markers %>% group_by(cluster) %>% top_n(n = 02, wt = avg_log2FC)
HCC_s_im.markers.top10 = HCC_s_im.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
HCC_s_im.markers.top20 = HCC_s_im.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

HCC_s_im@meta.data$Annotation1 <- rep(NA, dim(HCC_s_im@meta.data)[1])

HCC_s_im$Annotation1[HCC_s_im$seurat_clusters %in% c(0,2,4,6)] <- "T cell" #'APOB',"FGB","APOA1","APOE",
HCC_s_im$Annotation1[HCC_s_im$seurat_clusters %in% c(3,9)] <- "B cell" #'MYL9',"TAGLN","COL1A2","ACTA2",
HCC_s_im$Annotation1[HCC_s_im$seurat_clusters %in% c(1)] <- "Kuffer cell" #'APOE',"C1QC","C1QB","C1QA",
HCC_s_im$Annotation1[HCC_s_im$seurat_clusters %in% c(5)] <- "NK cell" ##"CCL4","PRF1","KLRD1",
HCC_s_im$Annotation1[HCC_s_im$seurat_clusters %in% c(7,8,10)] <- "Dendritic cell" #'CST3',"HLA-DRA","CLEC10A",

####metadata送还####

# 提取HCC_s_immune的Annotation1
immune_annotations <- HCC_s_im@meta.data$Annotation1
names(immune_annotations) <- rownames(HCC_s_im@meta.data)

# 创建一个空向量，初始值为HCC_s的Annotation
HCC_s$Annotation1 <- HCC_s$Annotation

# 用HCC_s_immune的Annotation1覆盖对应细胞
immune_cells <- intersect(rownames(HCC_s@meta.data), rownames(HCC_s_im@meta.data))
HCC_s$Annotation1[immune_cells] <- immune_annotations[immune_cells]

# 检查结果
table(HCC_s$Annotation1)


#####singleR####

#ref

#HCC_s

hpca.se <- readRDS("G:/SingleR_Celldex_ref/Human_General_HumanPrimaryCellAtlasData.rds")
hpca.se


sce_for_SingleR <- GetAssayData(HCC_s_im, layer ="data")
sce_for_SingleR

sce <- HCC_s_im
clusters=sce@meta.data$seurat_clusters
pred.hesc <- SingleR(test = sce_for_SingleR, ref = hpca.se, labels = hpca.se$label.fine,
                     clusters = clusters,  
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
table(pred.hesc$labels)
celltype_HCC_s_im = data.frame(ClusterID=rownames(pred.hesc), celltype=pred.hesc$labels, stringsAsFactors = F) 

gc()

HCC_s@meta.data$Annotation3 <- factor(
  HCC_s@meta.data$Annotation1, 
  levels = c("T cell","NK cell", "B cell","Dendritic cell","Kuffer cell","Hepatocyte","Endothelial","Fibroblast") # 按需调整为你的实际顺序
)

setwd("F:/4.panc_cancer_fig/1126")

saveRDS(HCC_s, file = "F:/4.panc_cancer_data/HCC_10_filter_1126.rds")

#####01####
pdf(paste('01.annodim.pdf',sep = '/'),width = 10,height=10)
DimPlot(HCC_s, reduction = "tsne", label = TRUE,raster=FALSE,group.by = "Annotation3")
dev.off()


#####02####
pdf(paste('02.annodim.pdf',sep = '/'),width = 10,height=10)

DimPlot(HCC_s, reduction = "tsne", label = TRUE,raster=FALSE,group.by = "HBV")
DimPlot(HCC_s, reduction = "tsne", label = TRUE,raster=FALSE,group.by = "Tumor")
DimPlot(HCC_s, reduction = "tsne", label = TRUE,raster=FALSE,group.by = "Patient")

dev.off()

#####03####
pdf(paste('03.annodot.pdf',sep = '/'),width = 15,height=8)
DotPlot(HCC_n, scale = F,
        features = c(
          "PTPRC",
          "CD3D","CD3E",
          "CCL4","PRF1","KLRD1",
          "CD79A","IGHM","MS4A1",
          'CST3',"HLA-DRA","CLEC10A",
          "C1QC","C1QB","C1QA",'APOE',
          'APOB',"FGB","APOA1",
          "PECAM1","AQP1","VWF","PLVAP",
          'MYL9',"TAGLN","COL1A2","ACTA2"
        ),
        group.by = "Annotation3",
        dot.scale = 8,
)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, color = "black"))+
  theme(axis.text.y = element_text(size = 10, color = "black"))+
  scale_color_gradientn(values = seq(0,1,0.5), colours = c("#FFFFFF","#FF0000","#CC0000"))

dev.off()

#####04####

HCC_s0 <- HCC_s

table(Idents(HCC_s0))
HCC_s <- subset(HCC_s0,idents = c("Fibroblast","Hepatocyte","Endothelial"))

scRNA = HCC_s
Idents(scRNA) <- scRNA@meta.data$Annotation1
ctys = levels(scRNA)
ctys

mycol <- colorRampPalette(rev(brewer.pal(n = 8, name ="Set1")))(length(ctys))

id <- HCC_s$HBV
clusters <- HCC_s$Annotation3
stim <- HCC_s$Patient_Part
df <- data.frame(stim = stim, clusters = clusters, id = id)
df$stim <- factor(df$stim,levels = c("P03_T","P03_N","P09_T","P09_N","P06_T","P06_N"))
df2 <- df %>% group_by(stim, clusters) %>% summarise(count = n())
# 调一下df2的levels顺序

# 为 HBV 信息添加颜色条所需数据
hbv_info <- df %>% distinct(stim, id)

# 自定义 HBV 条的颜色，设定两个颜色值供修改
hbv_colors <- c("Pos" = "red", "Neg" = "blue")

p2 <- ggplot(df2, aes(stim, count, fill = clusters, group = clusters)) +
  geom_bar(stat = 'identity', position = 'fill') +
  labs(x = '', y = 'Percentage') +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, color = 'black'), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    values = cols25()
  )

p2

#####bar####
# 创建主图表
p2 <- ggplot(df2, aes(stim, count, fill = clusters, group = clusters)) +
  geom_bar(stat = 'identity', position = 'fill') +
  labs(x = '', y = 'Percentage') +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, color = 'black'), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    values = cols25()#mycol
  )

# 创建 HBV 条图表
hbv_bar <- ggplot(hbv_info, aes(x = stim, y = 1, fill = id)) +
  geom_tile() +
  scale_fill_manual(values = hbv_colors) +
  theme_void() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# 合并主图和颜色条
library(patchwork)

final_plot <- p2 / hbv_bar + 
  plot_layout(heights = c(10, 0.5))  # 设置主图和颜色条的高度比例

# 输出最终图表
final_plot

pdf(paste('04.annopropor_immune.pdf',sep = '/'),width = 10,height=10)

ggplot(df2, aes(stim, count, fill = clusters, group = clusters)) +
  geom_bar(stat = 'identity', position = 'fill') +
  labs(x = '', y = 'Percentage') +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, color = 'black'), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    values = cols25()
  )

ggplot(hbv_info, aes(x = stim, y = 1, fill = id)) +
  geom_tile() +
  scale_fill_manual(values = hbv_colors) +
  theme_void() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

#final_plot
dev.off()

#####05####
pdf(paste('05.annopropor_all.pdf',sep = '/'),width = 10,height=10)
final_plot
dev.off()

#####06####
pdf(paste('06.annopropor_nonimmune.pdf',sep = '/'),width = 10,height=10)
ggplot(df2, aes(stim, count, fill = clusters, group = clusters)) +
  geom_bar(stat = 'identity', position = 'fill') +
  labs(x = '', y = 'Percentage') +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, color = 'black'), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    values = cols25()
  )

ggplot(hbv_info, aes(x = stim, y = 1, fill = id)) +
  geom_tile() +
  scale_fill_manual(values = hbv_colors) +
  theme_void() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
dev.off()

####kuffer+dc####
gc()

Idents(HCC_s0) <- HCC_s0$Annotation1
table(Idents(HCC_s0))
HCC_M <- subset(HCC_s0,idents=c("Kuffer cell","Dendritic cell"))

HCC_M <- NormalizeData(HCC_M, normalization.method = "LogNormalize", scale.factor = 10000)
HCC_M <- FindVariableFeatures(HCC_M, selection.method = "vst", nfeatures = 5000)
HCC_M <- ScaleData(HCC_M,features = rownames(HCC_M), assay = "RNA")
HCC_M <- RunPCA(HCC_M)
HCC_M <- RunHarmony(HCC_M, group.by.vars = "Patient_Part", dims.use = 1:20, lambda = 1,project.dim = F)
HCC_M <- RunUMAP(HCC_M, dim = 1:20, reduction = "harmony")
HCC_M <- FindNeighbors(HCC_M, dim = 1:20, reduction = "harmony")
HCC_M <- FindClusters(HCC_M, resolution = 0.2)
HCC_M <- RunTSNE(HCC_M, dims = 1:20,do.fast = TRUE, reduction = "harmony")

HCC_M

table(HCC_M$seurat_clusters)

HCC_M.markers <- FindAllMarkers(HCC_M, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
HCC_M.markers.top01 = HCC_M.markers %>% group_by(cluster) %>% top_n(n = 01, wt = avg_log2FC)
HCC_M.markers.top02 = HCC_M.markers %>% group_by(cluster) %>% top_n(n = 02, wt = avg_log2FC)
HCC_M.markers.top10 = HCC_M.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
HCC_M.markers.top20 = HCC_M.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

#####volcano access####

scRNA = HCC_M
scRNA@meta.data$seurat_clusters = Idents(scRNA)
ctys = levels(scRNA)
ctys

VlnPlot(HCC_M,group.by = "seurat_clusters",features=c('PTPRC'),ncol=1,pt.size = 0)
VlnPlot(HCC_M,group.by = "seurat_clusters",features=c('APOE',"C1QC","C1QB","C1QA"),ncol=1,pt.size = 0)

table(HCC_M$Annotation4,HCC_M$Annotation2)
VlnPlot(HCC_M,group.by = "seurat_clusters",features=c('APOE',"C1QC","C1QB","C1QA",
                                                      'CST3',"HLA-DRA","CLEC10A"),ncol=4,pt.size = 0) #kuffer + DC
VlnPlot(HCC_M,group.by = "Annotation4",features=c('APOE',"C1QC","C1QB","C1QA",
                                                      'CST3',"HLA-DRA","CLEC10A"),ncol=4,pt.size = 0) #kuffer + DC
VlnPlot(HCC_M,group.by = "Annotation4",features=c("ITGAX","CD83","THBD",
                                                  "CD1C",#cDC2 2 6
                                                  "CLEC9A",#cDC1 3
                                                  "IL3RA",#pDC(质粒样DC) 5
                                                  "XCR1"),ncol=4,pt.size = 0) #kuffer + DC

scRNA.markers <- FindAllMarkers(scRNA, min.pct = 0.25, 
                                logfc.threshold = 0.25)
head(scRNA.markers)

colnames(scRNA.markers)[6] = "seurat_clusters"
k = scRNA.markers$p_val_adj<0.05;table(k)

scRNA.markers = scRNA.markers[k,]

#上下调
scRNA.markers$label <- ifelse(scRNA.markers$avg_log2FC<0,"sigDown","sigUp")
topgene <- scRNA.markers %>%
  group_by(seurat_clusters) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  bind_rows(group_by(scRNA.markers, seurat_clusters) %>%
              top_n(n = 10, wt = -avg_log2FC))
head(topgene)

#根据log2FC范围确定背景柱长度：
dfbar = scRNA.markers %>%
  group_by(seurat_clusters) %>%
  summarise(low = round(min(avg_log2FC)-0.5),
            up = round(max(avg_log2FC)+0.5))

#绘制背景柱和散点图：
p1 <- ggplot()+
  geom_col(aes(x = seurat_clusters ,y = low),dfbar,
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(aes(x = seurat_clusters ,y = up),dfbar,
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(aes(x = seurat_clusters, y = avg_log2FC, color = label),scRNA.markers,
              width =0.4,size = 1)+
  scale_color_manual(values = c("#0077c0","#c72d2e"))+
  theme_classic()
p1

#X轴的色块标签：
library(RColorBrewer)
mycol <- colorRampPalette(rev(brewer.pal(n = 7, name ="Set1")))(length(ctys))
p2 <- p1 + 
  geom_tile(aes(x = ctys,y = 0),
            height = 0.5,fill = mycol, show.legend = F)+
  geom_text(aes(x= ctys, y = 0, label = ctys),
            size = 3,fontface = "bold")
p2

library(ggrepel)
#给每种细胞类型的top基因加上标签,调整细节：
p3 <- p2 + 
  geom_text_repel(aes(x = seurat_clusters,y = avg_log2FC,label = gene),
                  topgene,size = 3 )+
  labs(x = "CellType",y = "Average log2FoldChange",
       title = "Differential expression genes")+
  theme(
    plot.title = element_text(size = 14,color = "black",face = "bold"),
    axis.title = element_text(size = 12,color = "black",face = "bold"),
    axis.line.y = element_line(color = "black",linewidth = 0.8),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position.inside  = c(0.98,0.96),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 12)
  )+
  guides(color = guide_legend(override.aes = list(size = 4)))  
p3

#####07####
pdf(paste('07.annopropor_DCkuffer.pdf',sep = '/'),width = 10,height=10)
p3
dev.off()

#####08####
HCC_s <- HCC_M#subset(HCC_s0,idents = c("Fibroblast","Hepatocyte","Endothelial"))

scRNA = HCC_s
Idents(scRNA) <- scRNA@meta.data$seurat_clusters
ctys = levels(scRNA)
ctys

mycol <- colorRampPalette(rev(brewer.pal(n = 8, name ="Set1")))(length(ctys))

id <- HCC_s$HBV
clusters <- HCC_s$seurat_clusters
stim <- HCC_s$Patient_Part
df <- data.frame(stim = stim, clusters = clusters, id = id)
df$stim <- factor(df$stim,levels = c("P03_T","P03_N","P09_T","P09_N","P06_T","P06_N"))
df2 <- df %>% group_by(stim, clusters) %>% summarise(count = n())
# 调一下df2的levels顺序

# 为 HBV 信息添加颜色条所需数据
hbv_info <- df %>% distinct(stim, id)

# 自定义 HBV 条的颜色，设定两个颜色值供修改
hbv_colors <- c("Pos" = "red", "Neg" = "blue")

p2 <- ggplot(df2, aes(stim, count, fill = clusters, group = clusters)) +
  geom_bar(stat = 'identity', position = 'fill') +
  labs(x = '', y = 'Percentage') +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, color = 'black'), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    values = cols25()
  )

p2

#####bar####
# 创建主图表
p2 <- ggplot(df2, aes(stim, count, fill = clusters, group = clusters)) +
  geom_bar(stat = 'identity', position = 'fill') +
  labs(x = '', y = 'Percentage') +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, color = 'black'), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    values = cols25()#mycol
  )

# 创建 HBV 条图表
hbv_bar <- ggplot(hbv_info, aes(x = stim, y = 1, fill = id)) +
  geom_tile() +
  scale_fill_manual(values = hbv_colors) +
  theme_void() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# 合并主图和颜色条
library(patchwork)

final_plot <- p2 / hbv_bar + 
  plot_layout(heights = c(10, 0.5))  # 设置主图和颜色条的高度比例

# 输出最终图表
final_plot

pdf(paste('08.annopropor_DCkuffer.pdf',sep = '/'),width = 10,height=10)

ggplot(df2, aes(stim, count, fill = clusters, group = clusters)) +
  geom_bar(stat = 'identity', position = 'fill') +
  labs(x = '', y = 'Percentage') +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, color = 'black'), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    values = cols25()
  )

ggplot(hbv_info, aes(x = stim, y = 1, fill = id)) +
  geom_tile() +
  scale_fill_manual(values = hbv_colors) +
  theme_void() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

final_plot

#final_plot
dev.off()

####bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb####
####Go####

k = HCC_M.markers$p_val_adj<0.05;table(k)

HCC_M.markers = HCC_M.markers[k,]


table(HCC_M$seurat_clusters)
table(HCC_M.markers$cluster)

# 加载必要的包
library(clusterProfiler)
library(org.Hs.eg.db)  # 这里以人类基因为例

levels(HCC_M.markers$cluster)

length(levels(HCC_M.markers$cluster))

go_res <- list()

gc()

for (i in 1:length(levels(HCC_M.markers$cluster))) {
  
  # 将基因符号转换为Entrez ID
  gene_symbols <- HCC_M.markers$gene[which(HCC_M.markers$cluster==levels(HCC_M.markers$cluster)[i])]
  gene_entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # 执行GO富集分析
  go_results <- enrichGO( 
    gene = gene_entrez_ids$ENTREZID,  # 输入的Entrez基因ID
    OrgDb = org.Hs.eg.db,             # 选择基因组数据库
    keyType = "ENTREZID",    # 基因ID的类型
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
  
  go_res[[i]] <- as.data.frame(go_results)
}

names(go_res) <- c(0:(length(go_res)-1))
names(go_res) <- paste0("Myeloid_",names(go_res))
names(go_res)

# 使用barplot可视化GO富集结果
barplot(go_results, showCategory = 10)

# 使用dotplot可视化GO富集结果
dotplot(go_results, showCategory = 10)

####save####
# 安装并加载 writexl 包（如果尚未安装）
#install.packages("writexl")
#install.packages("openxlsx")
library(openxlsx)
library(writexl)

# 创建一个新的工作簿
wb <- createWorkbook()

# 将每个 data.frame 添加为工作表
for (name in names(go_res)) {
  addWorksheet(wb, name)
  writeData(wb, name, go_res[[name]])
}

# 保存为一个单一的 XLSX 文件
#saveWorkbook(wb, "F:\\4.panc_cancer_fig\\1126\\combined_data.xlsx", overwrite = TRUE)
write_xlsx(go_res, path = "F:\\4.panc_cancer_fig\\1126\\combined_data.xlsx")

####T/NK####
gc()

Idents(HCC_s0) <- HCC_s0$Annotation1
table(Idents(HCC_s0))
HCC_T <- subset(HCC_s0,idents=c("NK cell","T cell"))

HCC_T <- NormalizeData(HCC_T, normalization.method = "LogNormalize", scale.factor = 10000)
HCC_T <- FindVariableFeatures(HCC_T, selection.method = "vst", nfeatures = 5000)
HCC_T <- ScaleData(HCC_T,features = rownames(HCC_T), assay = "RNA")
HCC_T <- RunPCA(HCC_T)
HCC_T <- RunHarmony(HCC_T, group.by.vars = "Patient_Part", dims.use = 1:20, lambda = 1,project.dim = F)
HCC_T <- RunUMAP(HCC_T, dim = 1:20, reduction = "harmony")
HCC_T <- FindNeighbors(HCC_T, dim = 1:20, reduction = "harmony")
HCC_T <- FindClusters(HCC_T, resolution = 0.5)
HCC_T <- RunTSNE(HCC_T, dims = 1:20,do.fast = TRUE, reduction = "harmony")

HCC_T

table(HCC_T$seurat_clusters)

Idents(HCC_T) <- HCC_T$Annotation4
HCC_T.markers <- FindAllMarkers(HCC_T, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
HCC_T.markers.top01 = HCC_T.markers %>% group_by(cluster) %>% top_n(n = 01, wt = avg_log2FC)
HCC_T.markers.top02 = HCC_T.markers %>% group_by(cluster) %>% top_n(n = 02, wt = avg_log2FC)
HCC_T.markers.top10 = HCC_T.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
HCC_T.markers.top20 = HCC_T.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

VlnPlot(HCC_T,group.by = "seurat_clusters",features=c("CD3D","CD3E","CD4","CD8A","CD8B"),ncol=1,pt.size = 0) #T cell 0,1,10,11,14

VlnPlot(HCC_T,group.by = "seurat_clusters",features=c("CCL4","PRF1","KLRD1"),ncol=1,pt.size = 0) #NK 1,11

VlnPlot(HCC_T,group.by = "seurat_clusters",features=c("CD3D","CD3E","CD4","CD8A",
                                                      "CD8B","CCL4","PRF1","KLRD1"),
        ncol=2,pt.size = 0) #T cell 0,1,10,11,14

VlnPlot(HCC_T,group.by = "seurat_clusters",features=c("CCL4","PRF1","KLRD1"),ncol=1,pt.size = 0) #NK 1,11

#####09####
HCC_s <- HCC_T#subset(HCC_s0,idents = c("Fibroblast","Hepatocyte","Endothelial"))

scRNA = HCC_s
Idents(scRNA) <- scRNA@meta.data$seurat_clusters
ctys = levels(scRNA)
ctys

mycol <- colorRampPalette(rev(brewer.pal(n = 8, name ="Set1")))(length(ctys))

id <- HCC_s$HBV
clusters <- HCC_s$seurat_clusters
stim <- HCC_s$Patient_Part
df <- data.frame(stim = stim, clusters = clusters, id = id)
df$stim <- factor(df$stim,levels = c("P03_T","P03_N","P09_T","P09_N","P06_T","P06_N"))
df2 <- df %>% group_by(stim, clusters) %>% summarise(count = n())
# 调一下df2的levels顺序

# 为 HBV 信息添加颜色条所需数据
hbv_info <- df %>% distinct(stim, id)

# 自定义 HBV 条的颜色，设定两个颜色值供修改
hbv_colors <- c("Pos" = "red", "Neg" = "blue")

p2 <- ggplot(df2, aes(stim, count, fill = clusters, group = clusters)) +
  geom_bar(stat = 'identity', position = 'fill') +
  labs(x = '', y = 'Percentage') +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, color = 'black'), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    values = cols25()
  )

p2

#####bar####
# 创建主图表
p2 <- ggplot(df2, aes(stim, count, fill = clusters, group = clusters)) +
  geom_bar(stat = 'identity', position = 'fill') +
  labs(x = '', y = 'Percentage') +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, color = 'black'), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    values = cols25()#mycol
  )

# 创建 HBV 条图表
hbv_bar <- ggplot(hbv_info, aes(x = stim, y = 1, fill = id)) +
  geom_tile() +
  scale_fill_manual(values = hbv_colors) +
  theme_void() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# 合并主图和颜色条
library(patchwork)

final_plot <- p2 / hbv_bar + 
  plot_layout(heights = c(10, 0.5))  # 设置主图和颜色条的高度比例

# 输出最终图表
final_plot

pdf(paste('09.annopropor_DCkuffer.pdf',sep = '/'),width = 10,height=10)

ggplot(df2, aes(stim, count, fill = clusters, group = clusters)) +
  geom_bar(stat = 'identity', position = 'fill') +
  labs(x = '', y = 'Percentage') +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, color = 'black'), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    values = cols25()
  )

ggplot(hbv_info, aes(x = stim, y = 1, fill = id)) +
  geom_tile() +
  scale_fill_manual(values = hbv_colors) +
  theme_void() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

final_plot

VlnPlot(HCC_T,group.by = "seurat_clusters",features=c("CD3D","CD3E","CD4","CD8A",
                                                      "CD8B","CCL4","PRF1","KLRD1"),
        ncol=2,pt.size = 0) #T cell 0,1,10,11,14

#final_plot
dev.off()

#####annotation####

HCC_T@meta.data$Annotation4 <- rep(NA, dim(HCC_T@meta.data)[1])

HCC_T$Annotation4[HCC_T$seurat_clusters %in% c(0)] <- "T cell 1" #'APOB',"FGB","APOA1","APOE",

HCC_T$Annotation4[HCC_T$seurat_clusters %in% c(6)] <- "T cell 2" #'APOB',"FGB","APOA1","APOE",

HCC_T$Annotation4[HCC_T$seurat_clusters %in% c(2)] <- "CD4+ T cell" #'MYL9',"TAGLN","COL1A2","ACTA2",

HCC_T$Annotation4[HCC_T$seurat_clusters %in% c(1)] <- "CD8+ T cell 1" #"PECAM1","AQP1","VWF","PLVAP",
HCC_T$Annotation4[HCC_T$seurat_clusters %in% c(3)] <- "CD8+ T cell 2" #"PECAM1","AQP1","VWF","PLVAP",
HCC_T$Annotation4[HCC_T$seurat_clusters %in% c(5)] <- "CD8+ T cell 3" #"PECAM1","AQP1","VWF","PLVAP",
HCC_T$Annotation4[HCC_T$seurat_clusters %in% c(7)] <- "CD8+ T cell 4" #"PECAM1","AQP1","VWF","PLVAP",

HCC_T$Annotation4[HCC_T$seurat_clusters %in% c(4)] <- "NK cell"

#####10####
HCC_s <- HCC_T#subset(HCC_s0,idents = c("Fibroblast","Hepatocyte","Endothelial"))
HCC_s$seurat_clusters <- HCC_s$Annotation4
HCC_s$seurat_clusters <- factor(HCC_s$seurat_clusters,levels = c("NK cell","T cell 1","T cell 2",
                                                                 "CD4+ T cell","CD8+ T cell 1","CD8+ T cell 2",
                                                                 "CD8+ T cell 3","CD8+ T cell 4"
                                                                 ))


scRNA = HCC_s
Idents(scRNA) <- scRNA@meta.data$seurat_clusters
ctys = levels(scRNA)
ctys

my_colors <- c('#b1d2aa', '#86b880', '#067b41', '#898f32', '#6a984c', '#1b5228',
               '#e8e842', '#94c54e', '#a2d6f0', '#0785ba', '#404999', '#c12740')
mycol <- colorRampPalette(rev(brewer.pal(n = 8, name ="Set1")))(length(ctys))

id <- HCC_s$HBV
clusters <- HCC_s$seurat_clusters
stim <- HCC_s$Patient_Part
df <- data.frame(stim = stim, clusters = clusters, id = id)
df$stim <- factor(df$stim,levels = c("P03_T","P03_N","P09_T","P09_N","P06_T","P06_N"))
df2 <- df %>% group_by(stim, clusters) %>% summarise(count = n())
# 调一下df2的levels顺序

# 为 HBV 信息添加颜色条所需数据
hbv_info <- df %>% distinct(stim, id)

# 自定义 HBV 条的颜色，设定两个颜色值供修改
hbv_colors <- c("Pos" = "#fbcf41", "Neg" = "#5966b1")

p2 <- ggplot(df2, aes(stim, count, fill = clusters, group = clusters)) +
  geom_bar(stat = 'identity', position = 'fill') +
  labs(x = '', y = 'Percentage') +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, color = 'black'), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    values = my_colors[-c(4,5,7,6)]
  )

p2

#####bar####
# 创建主图表
p2 <- ggplot(df2, aes(stim, count, fill = clusters, group = clusters)) +
  geom_bar(stat = 'identity', position = 'fill') +
  labs(x = '', y = 'Percentage') +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, color = 'black'), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    values = my_colors[-c(4,5,7,6)]#mycol
  )

# 创建 HBV 条图表
hbv_bar <- ggplot(hbv_info, aes(x = stim, y = 1, fill = id)) +
  geom_tile() +
  scale_fill_manual(values = hbv_colors) +
  theme_void() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# 合并主图和颜色条
library(patchwork)

final_plot <- p2 / hbv_bar + 
  plot_layout(heights = c(10, 0.5))  # 设置主图和颜色条的高度比例

# 输出最终图表
final_plot

pdf(paste('10.annopropor_DCkuffer.pdf',sep = '/'),width = 10,height=10)


final_plot

VlnPlot(HCC_T,group.by = "Annotation4",features=c("CD3D","CD3E","CD4","CD8A",
                                                      "CD8B","CCL4","PRF1","KLRD1"),
        ncol=2,pt.size = 0) #T cell 0,1,10,11,14

#final_plot
dev.off()

#####volcano####

scRNA = HCC_T
scRNA@meta.data$Annotation4 <- factor(scRNA@meta.data$Annotation4,levels = c("NK cell","T cell 1","T cell 2",
                                                                               "CD4+ T cell","CD8+ T cell 1","CD8+ T cell 2",
                                                                               "CD8+ T cell 3","CD8+ T cell 4"
))
Idents(scRNA) <- scRNA@meta.data$Annotation4 
ctys <-  levels(scRNA$Annotation4)
ctys

#VlnPlot(HCC_T,group.by = "seurat_clusters",features=c('PTPRC'),ncol=1,pt.size = 0)
#VlnPlot(HCC_T,group.by = "seurat_clusters",features=c('APOE',"C1QC","C1QB","C1QA"),ncol=1,pt.size = 0)

scRNA.markers <- FindAllMarkers(scRNA, min.pct = 0.25, 
                                logfc.threshold = 0.25)
head(scRNA.markers)

colnames(scRNA.markers)[6] = "seurat_clusters"
k = scRNA.markers$p_val_adj<0.05;table(k)

scRNA.markers = scRNA.markers[k,]

#上下调
scRNA.markers$label <- ifelse(scRNA.markers$avg_log2FC<0,"sigDown","sigUp")
topgene <- scRNA.markers %>%
  group_by(seurat_clusters) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  bind_rows(group_by(scRNA.markers, seurat_clusters) %>%
              top_n(n = 10, wt = -avg_log2FC))
head(topgene)

#根据log2FC范围确定背景柱长度：
dfbar = scRNA.markers %>%
  group_by(seurat_clusters) %>%
  summarise(low = round(min(avg_log2FC)-0.5),
            up = round(max(avg_log2FC)+0.5))

#绘制背景柱和散点图：
p1 <- ggplot()+
  geom_col(aes(x = seurat_clusters ,y = low),dfbar,
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(aes(x = seurat_clusters ,y = up),dfbar,
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(aes(x = seurat_clusters, y = avg_log2FC, color = label),scRNA.markers,
              width =0.4,size = 1)+
  scale_color_manual(values = c("#0077c0","#c72d2e"))+
  theme_classic()
p1

#X轴的色块标签：
library(RColorBrewer)
mycol <- colorRampPalette(rev(brewer.pal(n = 7, name ="Set1")))(length(ctys))
p2 <- p1 + 
  geom_tile(aes(x = ctys,y = 0),
            height = 0.5,fill = mycol, show.legend = F)+
  geom_text(aes(x= ctys, y = 0, label = ctys),
            size = 3,fontface = "bold")
p2

library(ggrepel)
#给每种细胞类型的top基因加上标签,调整细节：
p3 <- p2 + 
  geom_text_repel(aes(x = seurat_clusters,y = avg_log2FC,label = gene),
                  topgene,size = 3 )+
  labs(x = "CellType",y = "Average log2FoldChange",
       title = "Differential expression genes")+
  theme(
    plot.title = element_text(size = 14,color = "black",face = "bold"),
    axis.title = element_text(size = 12,color = "black",face = "bold"),
    axis.line.y = element_line(color = "black",linewidth = 0.8),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position.inside  = c(0.98,0.96),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 12)
  )+
  guides(color = guide_legend(override.aes = list(size = 4)))  
p3

#####11####
pdf(paste('11.annopropor_TNK.pdf',sep = '/'),width = 10,height=10)
p3
dev.off()

####GO T/NK####

k = HCC_T.markers$p_val_adj<0.05;table(k)

HCC_T.markers = HCC_T.markers[k,]


table(HCC_T$Annotation4)
table(HCC_T.markers$cluster)

# 加载必要的包
library(clusterProfiler)
library(org.Hs.eg.db)  # 这里以人类基因为例

levels(HCC_T.markers$cluster)

length(levels(HCC_T.markers$cluster))

go_res <- list()

gc()

for (i in 1:length(levels(HCC_T.markers$cluster))) {
  
  # 将基因符号转换为Entrez ID
  gene_symbols <- HCC_T.markers$gene[which(HCC_T.markers$cluster==levels(HCC_T.markers$cluster)[i])]
  gene_entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # 执行GO富集分析
  go_results <- enrichGO( 
    gene = gene_entrez_ids$ENTREZID,  # 输入的Entrez基因ID
    OrgDb = org.Hs.eg.db,             # 选择基因组数据库
    keyType = "ENTREZID",    # 基因ID的类型
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
  
  go_res[[i]] <- as.data.frame(go_results)
}

names(go_res) <- levels(HCC_T.markers$cluster)

# 使用barplot可视化GO富集结果
barplot(go_results, showCategory = 10)

# 使用dotplot可视化GO富集结果
dotplot(go_results, showCategory = 10)

####save####
# 安装并加载 writexl 包（如果尚未安装）
#install.packages("writexl")
#install.packages("openxlsx")
library(openxlsx)
library(writexl)

# 创建一个新的工作簿
wb <- createWorkbook()

# 将每个 data.frame 添加为工作表
for (name in names(go_res)) {
  addWorksheet(wb, name)
  writeData(wb, name, go_res[[name]])
}

# 保存为一个单一的 XLSX 文件
#saveWorkbook(wb, "F:\\4.panc_cancer_fig\\1126\\combined_data.xlsx", overwrite = TRUE)
write_xlsx(go_res, path = "F:\\4.panc_cancer_fig\\1126\\combined_data_TNK.xlsx")

####保存，准备去做细胞通讯####

table(HCC_s0$Annotation1)
table(HCC_M$seurat_clusters)
table(HCC_T$Annotation4)

VlnPlot(HCC_M,group.by = "Annotation4",features=c('APOE',"C1QC","C1QB","C1QA",
                                                  'CST3',"HLA-DRA","CLEC10A"),ncol=4,pt.size = 0) #kuffer + DC

table(HCC_M$Annotation4,HCC_M$Annotation1)

VlnPlot(HCC_M,group.by = "Annotation4",features=c("ITGAX","CD83","THBD",
                                                  "CD1C",#cDC2 2 6
                                                  "CLEC9A",#cDC1 3
                                                  "IL3RA",#pDC(质粒样DC) 5
                                                  "XCR1"),ncol=4,pt.size = 0) #kuffer + DC

table(HCC_M$Annotation)
table(HCC_M$Annotation1)
table(HCC_M$Annotation2)
table(HCC_M$Annotation3)
table(HCC_M$Annotation4)

HCC_M$Annotation4 <- paste0("Myeloid_",HCC_M$seurat_clusters)

HCC_M@meta.data$Annotation3 <- rep(NA, dim(HCC_M@meta.data)[1])

HCC_M$Annotation3[HCC_M$seurat_clusters %in% c(5)] <- "pDC" #'IL3RA",CD123
HCC_M$Annotation3[HCC_M$seurat_clusters %in% c(3)] <- "cDC1" #'CLEC9A",
HCC_M$Annotation3[HCC_M$seurat_clusters %in% c(2,6)] <- "cDC2" #"CD1C",
HCC_M$Annotation3[HCC_M$seurat_clusters %in% c(0)] <- "Kuffer 1" 
HCC_M$Annotation3[HCC_M$seurat_clusters %in% c(1)] <- "Kuffer 2" 
HCC_M$Annotation3[HCC_M$seurat_clusters %in% c(4)] <- "Kuffer 3" 

table(HCC_M$Annotation4,HCC_M$Annotation3)

table(HCC_s0$Annotation)#Raw_type
table(HCC_s0$Annotation1)#Major_type
table(HCC_s0$Annotation2)#原始注释
table(HCC_s0$Annotation3)#Minor_factor

HCC_s0$Annotation2 <- HCC_s0$Annotation1

####metadata送还####

HCC_s <- HCC_s0

table(HCC_s$Annotation1)

# 检查结果
table(HCC_s$Annotation)
table(HCC_s$Annotation1)
table(HCC_s$Annotation2)
table(HCC_s$Annotation3)
table(HCC_M$Annotation4)
table(HCC_M$Annotation3)
table(HCC_T$Annotation4)

# 提取HCC_s_immune的Annotation1
T_annotations <- HCC_T@meta.data$Annotation4
names(T_annotations) <- rownames(HCC_T@meta.data)

myeloid_annotations <- HCC_M@meta.data$Annotation3
names(myeloid_annotations) <- rownames(HCC_M@meta.data)

# 创建一个空向量，初始值为HCC_s的Annotation
HCC_s$Annotation2 <- HCC_s$Annotation1

# 用HCC_s_immune的Annotation1覆盖对应细胞
T_cells <- intersect(rownames(HCC_s@meta.data), rownames(HCC_T@meta.data))
HCC_s$Annotation2[T_cells] <- T_annotations[T_cells]

M_cells <- intersect(rownames(HCC_s@meta.data), rownames(HCC_M@meta.data))
HCC_s$Annotation2[M_cells] <- myeloid_annotations[M_cells]

table(HCC_s$Annotation)#major
table(HCC_s$Annotation1)#minor
table(HCC_s$Annotation2)#subtype
table(HCC_s$Annotation3)#minor_raw

HCC_s$Annotation1[which(HCC_s$Annotation1 %in% c("T cell","NK cell"))] <- "T/NK cell"
HCC_s$Annotation1[which(HCC_s$Annotation1 %in% c("Kuffer cell","Dendritic cell"))] <- "Myeloid cell"

table(HCC_s$Annotation1,HCC_s$Annotation3)

####save####

###figure
###cellchat
###TF
###monocle

table(HCC_s$Annotation)#major
table(HCC_s$Annotation1)#minor
table(HCC_s$Annotation2)#subtype
table(HCC_s$Annotation3)#minor_raw

colnames(HCC_s@meta.data)
colnames(HCC_s@meta.data)[16:19] <- c("Major_celltype","Minor_celltype","Subtype_celltype","Raw_anno")
colnames(HCC_s@meta.data)[16:19] 

saveRDS(HCC_s, file = "F:/4.panc_cancer_data/HCC_10_filter_1127.rds")


































