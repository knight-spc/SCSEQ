library(dplyr)
library(monocle)
library(pals)
library(Seurat)
library(RColorBrewer)
library(AnnoProbe)
library(optparse)

# 定义命令行参数
option_list <- list(
  make_option(c("-w", "--wd"), type="character", default=NULL, 
              help="输入的工作路径"),
  make_option(c("-s", "--SubsetList"), type="character", default=NULL, 
              help="选择要观察的细胞群体")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 检查必要参数
if (is.null(opt$wd)) {
  stop("请提供文件路径 (-i/--wd)")
}

# 读取数据
setwd(opt$wd)
heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'RdBu')))

stdata=readRDS('seurat_data.rds')

stdata$annot <- stdata$celltype
table(stdata@meta.data$annot)


if (!is.null(opt$SubsetList)) {
    # 将输入的字符串按逗号分割并去除空格
    selected_clusters <- strsplit(opt$SubsetList, ",")[[1]]
    selected_clusters <- trimws(selected_clusters)  # 去除空格
    # seurat_object <- subset(seurat_object, idents=selected_clusters)
}
print(111111111)
##！！！！！！建议你给他写个输入细胞判定， >5000细胞，就subset downsample 5000；
## <=5000细胞数可以直接做；细胞数最少的那一群有潜在报错可能

##！！！！！！也可以下面代码自选细胞数
#stdata_new_0=subset(stdata,(annot=='0'))
#stdata_new_0=subset(stdata_new_0,cells=sample(colnames(stdata_new_0) ,300))
#stdata_new_1...........
#stdata_new_2............
#stdata_new_MAC = merge(stdata_new_0,y=c(stdata_new_1,stdata_new_3,stdata_new_4,stdata_new_5))

# # 为每个选定的cluster随机选择5000个细胞
# stdata_subset_list <- lapply(selected_clusters, function(cluster) {
#   subset_data <- subset(stdata, annot == cluster)
#   # 如果该cluster的细胞数小于5000，则全部保留
#   cell_count <- min(5000, ncol(subset_data))
#   subset(subset_data, cells = sample(colnames(subset_data), cell_count))
# })
# # 合并所有选定的细胞
# stdata_new_MAC <- merge(stdata_subset_list[[1]], 
#                        y = stdata_subset_list[2:length(stdata_subset_list)])

# selected_clusters <- c('B cells', 'T cells', 'Mig.cDCs')  # 这里填入您想要的cluster
# 为每个选定的cluster随机选择细胞
min_cells_per_cluster <- 50  # 设置每个群体最小细胞数

stdata_subset_list <- lapply(selected_clusters, function(cluster) {
  subset_data <- subset(stdata, annot == cluster)
  if (ncol(subset_data) == 0) {
    warning(paste("未找到群体:", cluster))
    return(NULL)
  }
  if (ncol(subset_data) < min_cells_per_cluster) {
    warning(paste("群体", cluster, "细胞数不足", min_cells_per_cluster, "个，将被跳过"))
    return(NULL)
  }
  # 如果该cluster的细胞数小于5000，则全部保留
  cell_count <- min(5000, ncol(subset_data))
  subset(subset_data, cells = sample(colnames(subset_data), cell_count))
})

# 移除NULL值
stdata_subset_list <- Filter(Negate(is.null), stdata_subset_list)

# 检查是否有有效的数据集
if (length(stdata_subset_list) == 0) {
  stop("没有找到符合条件的细胞")
}

# 合并所有选定的细胞
if (length(stdata_subset_list) == 1) {
  stdata_new_MAC <- stdata_subset_list[[1]]
} else {
  stdata_new_MAC <- merge(stdata_subset_list[[1]], 
                         y = stdata_subset_list[2:length(stdata_subset_list)])
}

#直接用投入数据
# stdata_new_MAC=stdata

# ============== 3. 准备Monocle输入数据 ==============
# 获取所有层的名称
layer_names <- grep("^counts", names(stdata_new_MAC@assays$RNA@layers), value = TRUE)

# 提取并合并所有层的数据
data_list <- lapply(layer_names, function(layer) {
  as(as.matrix(LayerData(stdata_new_MAC, assay = "RNA", layer = layer)), 'sparseMatrix')
})

# 合并所有层的数据
data <- do.call(cbind, data_list)
print(dim(data))  # 打印合并后的数据维度

# 在构建CellDataSet之前，先过滤低表达基因
# 保留在至少10个细胞中表达的基因
keep_genes <- rowSums(data > 0) >= 10
data <- data[keep_genes, ]
pd <- new('AnnotatedDataFrame', data = stdata_new_MAC@meta.data[colnames(data), ])
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
print(dim(data))  # 打印合并后的数据维度
# ============== 4. 构建CellDataSet对象 ==============
# 使用负二项分布创建CellDataSet对象
HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())
print(222222222)

# ============== 5. 数据质控和标准化 ==============
# 估计细胞大小因子
HSMM <- estimateSizeFactors(HSMM)
# 可以考虑在estimateDispersions之前添加数据过滤和预处理
# 过滤低质量基因和细胞
HSMM <- detectGenes(HSMM, min_expr = 3)  # 检测基因表达----------------------------------------------------------------
# HSMM <- HSMM[fData(HSMM)$num_cells_expressed >= 20,  # 保留在至少20个细胞中表达的基因
#              pData(HSMM)$num_genes_expressed >= 200]  # 保留表达至少200个基因的细胞

# 估计基因表达的离散度
HSMM <- estimateDispersions(HSMM)  # 使用local拟合
# HSMM <- detectGenes(HSMM, min_expr = 3) 
print(head(fData(HSMM)))
# 获取表达基因列表
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))

# ============== 6. 差异基因分析 ==============
# 设置细胞分组信息
Idents(stdata_new_MAC) <- stdata_new_MAC@meta.data$annot
DefaultAssay(stdata_new_MAC)='RNA'
# all.markers <- FindAllMarkers(object = stdata_new_MAC)######一般不用这个函数找marker,网上找其他教程
all.markers = read.csv("all_marker.csv")  ### 这里选择直接读


###！！！！！！！我一般这一步检查有多少个marker，然后下一步diff_markers通过avg_log2FC调整数目，能影响最终轨迹
# all.markers

# 筛选显著差异基因
diff_markers=subset(all.markers,p_val_adj<0.1 & avg_log2FC>0.5)

# ============== 7. 轨迹分析 ==============
# 设置用于排序的基因
ordering_genes=row.names(diff_markers)
HSMM <- setOrderingFilter(HSMM, ordering_genes)
# 可视化排序基因
plot_ordering_genes(HSMM)

# 使用DDRTree方法进行降维
# HSMM <- reduceDimension(HSMM, 
#                        max_components = 2,
#                        num_dim = 30,
#                        method = 'DDRTree',
#                        auto_param_selection = TRUE,  # 自动参数选择
#                        verbose = TRUE)  # 显示进度
HSMM <- reduceDimension(HSMM, max_components = 2,
                        num_dim = 20,method = 'DDRTree',norm_method = 'log',) # DDRTree方式！！！！！！！可改ICA

# 对细胞进行排序
HSMM <- orderCells(HSMM, root_state = NULL, num_paths = NULL, reverse = FALSE)

# ============== 8. 可视化 ==============
# 绘制伪时间轨迹图
plot1<-plot_cell_trajectory(HSMM, color_by = "Pseudotime",cell_size=2,show_backbone =FALSE,show_branch_points=FALSE)+
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        axis.title = element_text(size=20),
        legend.key.size = unit(1, 'cm'),legend.text=element_text(size=10),
        legend.title=element_blank())+scale_color_gradientn(colours = heatmap_Palette(100))
plot1

# 绘制按注释着色的轨迹图
# 获取实际的cluster数量
n_clusters <- length(unique(HSMM$annot))

# 生成足够的颜色
color_palette <- c('#FFD700','#8B0000','#7B68EE','#FF5733','#5DCE7E','#3F48CC',
                  '#FF69B4','#4B0082','#00FF7F','#FF4500','#1E90FF',
                  '#FF1493','#32CD32','#FF8C00','#4169E1')

plot2 <- plot_cell_trajectory(HSMM, color_by = "annot", cell_size=2,
                            show_backbone = FALSE, show_branch_points = FALSE) +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank(),
        axis.title = element_text(size=20),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size=10),
        legend.title = element_blank()) +
  scale_color_manual(values = color_palette[1:n_clusters]) +
  guides(colour = guide_legend(override.aes = list(size=8)))
plot2

# 绘制分面的轨迹图
plot3 <- plot_cell_trajectory(HSMM, color_by = "annot", cell_size=2,
                            show_backbone = FALSE, show_branch_points = FALSE) +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank(),
        axis.title = element_text(size=20),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size=10),
        legend.title = element_blank()) +
  scale_color_manual(values = color_palette[1:n_clusters]) +
  guides(colour = guide_legend(override.aes = list(size=8))) +
  facet_wrap(~annot)
plot3

# 保存图片
ggsave(plot1,file='001.trajectory.pdf',width=6,height=6)
ggsave(plot2,file='001.annot_trajectory.pdf',width=6,height=6)
ggsave(plot3,file='001.annot_trajectory_split.pdf',width=12,height=12)
print("图片保存完成")

# # ============== 9. 差异基因分析 ==============
# # 确保基因名称存在于HSMM对象中
# valid_genes <- intersect(unique(diff_markers$gene), rownames(HSMM))
# if(length(valid_genes) == 0) {
#     stop("没有找到匹配的基因名称")
# }

# # 使用有效的基因名称进行差异检验
# diff_test_res <- differentialGeneTest(HSMM[valid_genes,],
#                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
# sig_gene_names <- row.names(subset(diff_test_res, qval < 0.000001))
# length(sig_gene_names)

# # ============== 10. 热图可视化 ==============
# # 绘制伪时间热图
# res1=plot_pseudotime_heatmap(HSMM[sig_gene_names,],
#                              num_clusters = 3,
#                              cores = 1,
#                              show_rownames = T,return_heatma=T,hmcols = heatmap_Palette(100))
# # 保存热图
# pdf(paste('002.pseudotime_heatmap.pdf',sep = '/'), width = 8, height = 30)
# res1
# dev.off()

# saveRDS(HSMM, file = "HSMM_object.rds")

# ============== 11. 后续分析建议 ==============
##!!!!!!!!!!!!!后面还有大量个性化的内容，目前我没做


###比如每个cluster的结果做GO


##基因在轨迹上变化！！！！！！！！这个比较个性化，也比较常用
# example_genes <- row.names(subset(fData(HSMM),
#                           gene_short_name %in% c("SPP1",'LGMN','A2M','GPNMB','TIMP1','APOE')))
# cds_subset <- HSMM[example_genes,]
# plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime",ncol = 1,relative_expr=FALSE)+scale_color_gradientn(colours = heatmap_Palette(100))


execution_time <- system.time({
  HSMM_restored <- readRDS("HSMM_object.rds")
  example_genes <- row.names(subset(fData(HSMM_restored),
                                  gene_short_name %in% c("SPP1",'LGMN','A2M','GPNMB','TIMP1','APOE')))
  cds_subset <- HSMM_restored[example_genes,]
  plot_genes_in_pseudotime(cds_subset, 
                          color_by = "Pseudotime",
                          ncol = 1,
                          relative_expr = FALSE) +
    scale_color_gradientn(colours = heatmap_Palette(100))
})

print(execution_time)  # 输出用户时间、系统时间和实际时间