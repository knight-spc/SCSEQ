print('111111111')
library(Seurat)
library(ggplot2)
library(pheatmap)
library(grid)
library(dplyr) 
library(optparse)

print('111111111')

# 定义命令行参数
option_list <- list(
  make_option(c("-w", "--wd"), type="character", default=NULL, 
              help="输入的工作路径"),
  make_option(c("-t", "--time"), type="character", default=NULL, 
              help="请输入时间戳"),
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
setwd(opt$wd)
data <- readRDS("seurat_data.rds")

time <- opt$time 
tf <- opt$SubsetList
tf_list <- as.character(unlist(strsplit(tf,split=",")))
print("输入的TF列表：")
print(tf_list)


rownames(data@meta.data) <- Cells(data)
colnames(data) <- Cells(data)

setwd("SCENIC")
auc <- read.csv("auc_mtx.csv", row.names=1)
data[["SCENIC"]] <- CreateAssayObject(
  data = t(as.matrix(auc))
)

markers <- FindAllMarkers(data, assay="SCENIC", only.pos = TRUE)
markers.top = markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)
top_regulons <- markers.top %>%
  pull(gene) %>%
  unique()


DefaultAssay(data) <- "SCENIC"

# 1. UMAP
p1 <- FeaturePlot(
  data,
  features = tf_list, # c("IRF7","CEBPD","TCF7")
  reduction = "umap"
) & theme_bw()

ggsave(
  filename = paste0("1.SCENIC_UMAP_", time, ".pdf"),
  plot = p1,
  width = 8,
  height = 6
)


# 2. DotPlot
top_regulons <- unique(top_regulons)

p2 <- DotPlot(
  data,
  assay = "SCENIC",
  features = top_regulons,
  group.by = "celltype"
) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

ggsave(
  filename = paste0("2.SCENIC_DotPlot_", time, ".pdf"),
  plot = p2,
  width = 18,
  height = 7
)

# 3. Heatmap
library(pheatmap)

avg <- AverageExpression(
  data,
  assays = "SCENIC",
  group.by = "celltype"
)$SCENIC

valid_regulons <- intersect(top_regulons, rownames(avg))
mat <- avg[valid_regulons, , drop = FALSE]

p <- pheatmap(mat)

ggsave(
  filename = paste0("3.SCENIC_Heatmap_", time, ".pdf"),
  plot = p$gtable,
  width = 8,
  height = 10
)

# 4. Violin
p4 <- VlnPlot(
  data,
  features = tf_list,
  assay = "SCENIC",
  group.by = "celltype"
)
ggsave(
  filename = paste0("4.SCENIC_Violin_", time, ".pdf"),
  plot = p4,
  width = 8,
  height = 6
)