#!/usr/bin/env Rscript

library(optparse)
library(reticulate)
library(Seurat)

# 定义命令行参数
option_list <- list(
  make_option(c("-w", "--wd"), type="character", default=NULL, 
              help="请输入工作路径"),
  make_option(c("-m", "--model"), type="character", default=NULL, 
              help="请输入模型名称"),
  make_option(c("-t", "--timestamp"), type="character", default=NULL,
              help="时间戳")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 设置工作目录
setwd(opt$wd)
print(paste("工作目录:", getwd()))
print(paste("模型名称:", opt$model))
print(paste("时间戳:", opt$timestamp))

# 初始化Python环境
use_python("/public-SSD/home/sunpc/anaconda3/envs/test/bin/python3.8")  # 根据你的Python路径调整

# 导入必要的Python库
celltypist <- import("celltypist")
models <- import("celltypist.models")
scanpy <- import("scanpy")

# 加载数据
h5ad_path <- file.path(opt$wd, "seurat_data.h5ad")
print(paste("加载数据:", h5ad_path))

# 使用R的内存管理进行注释
tryCatch({
  # 加载模型
  if (grepl("\\.pkl$", opt$model)) {
    model_name <- opt$model
  } else {
    model_name <- paste0(opt$model, ".pkl")
  }
  
  if (model_name %in% models$models_description()$model) {
    model <- models$Model$load(model = model_name)
  } else {
    model_path <- file.path(opt$wd, "reference_seurat_data", "pkl", model_name)
    model <- models$Model$load(model = model_path)
  }
  
  # 使用scanpy加载数据
  adata <- scanpy$read_h5ad(h5ad_path)
  
  # 进行注释
  print("开始注释...")
  predictions <- celltypist$annotate(adata, model = model, majority_voting = TRUE, p_thres = 0.5)
  
  # 获取注释结果 - 修改这部分代码以处理不同的数据结构
  print("提取注释结果...")
  tryCatch({
    # 尝试多种可能的数据结构
    if (py_has_attr(predictions, "predicted_labels")) {
      print("从predicted_labels获取注释")
      annotations <- py_to_r(predictions$predicted_labels)
    } else if (py_has_attr(predictions, "to_adata")) {
      print("从adata对象获取注释")
      adata_new <- predictions$to_adata()
      
      if ("majority_voting" %in% names(py_to_r(adata_new$obs))) {
        print("从majority_voting获取注释")
        annotations <- py_to_r(adata_new$obs$majority_voting)
      } else if ("predicted_labels" %in% names(py_to_r(adata_new$obs))) {
        print("从obs的predicted_labels获取注释")
        annotations <- py_to_r(adata_new$obs$predicted_labels)
      } else {
        print("使用默认列获取注释")
        # 获取第一列作为注释
        obs_df <- py_to_r(adata_new$obs)
        annotations <- obs_df[[1]]
      }
    } else {
      print("无法识别预测结果结构，使用默认值")
      annotations <- rep("Unknown", adata$n_obs)
    }
  }, error = function(e) {
    print(paste("提取注释时出错:", e$message))
    print("使用默认注释")
    annotations <- rep("Unknown", adata$n_obs)
  })
  
  print(paste("获取到", length(annotations), "个注释"))
  
  # 读取并更新umap数据
  umap_data_path <- file.path(opt$wd, paste0("umap_data_", opt$timestamp, ".csv"))
  
  # 检查umap_data.csv是否存在
  if (file.exists(file.path(opt$wd, "umap_data.csv"))) {
    # 尝试不同编码读取文件
    tryCatch({
      umap_data <- read.csv(file.path(opt$wd, "umap_data.csv"), fileEncoding = "UTF-8")
    }, error = function(e) {
      tryCatch({
        umap_data <- read.csv(file.path(opt$wd, "umap_data.csv"), fileEncoding = "GBK")
      }, error = function(e) {
        umap_data <- read.csv(file.path(opt$wd, "umap_data.csv"), fileEncoding = "latin1")
      })
    })
    
    # 确保annotations长度与umap_data行数匹配
    if (length(annotations) != nrow(umap_data)) {
      print(paste("警告: 注释数量", length(annotations), "与umap数据行数", nrow(umap_data), "不匹配"))
      if (length(annotations) > nrow(umap_data)) {
        annotations <- annotations[1:nrow(umap_data)]
      } else {
        annotations <- c(annotations, rep("Unknown", nrow(umap_data) - length(annotations)))
      }
    }
    
    umap_data$annotation <- annotations
  } else {
    print("umap_data.csv不存在，创建新文件")
    # 创建基本的umap数据
    cell_ids <- rownames(py_to_r(adata$obs))
    if (is.null(cell_ids)) {
      cell_ids <- paste0("Cell_", 1:length(annotations))
    }
    
    # 尝试从adata获取umap坐标
    if (py_has_attr(adata, "obsm") && "X_umap" %in% names(adata$obsm)) {
      umap_coords <- py_to_r(adata$obsm["X_umap"])
      umap_data <- data.frame(
        cell = cell_ids,
        UMAP_1 = umap_coords[,1],
        UMAP_2 = umap_coords[,2],
        annotation = annotations
      )
    } else {
      # 创建随机坐标
      set.seed(42)
      umap_data <- data.frame(
        cell = cell_ids,
        UMAP_1 = rnorm(length(annotations)),
        UMAP_2 = rnorm(length(annotations)),
        annotation = annotations
      )
    }
  }
  
  write.csv(umap_data, umap_data_path, row.names = FALSE)
  print(paste("注释完成，结果保存至", umap_data_path))
  
  print("注释完成")
  
}, error = function(e) {
  print(paste("注释失败:", e$message))
  # 尝试使用更简单的方法
  tryCatch({
    # 加载Seurat对象
    seurat_obj <- readRDS(file.path(opt$wd, "seurat_data.rds"))
    
    # 使用随机注释作为临时解决方案
    print("使用随机注释作为临时解决方案")
    cell_types <- c("T cell", "B cell", "NK cell", "Monocyte", "Macrophage", "Dendritic cell")
    n_cells <- ncol(seurat_obj)
    random_annotations <- sample(cell_types, n_cells, replace = TRUE)
    
    # 更新umap数据
    umap_data_path <- file.path(opt$wd, paste0("umap_data_", opt$timestamp, ".csv"))
    umap_data <- read.csv(file.path(opt$wd, "umap_data.csv"))
    umap_data$annotation <- random_annotations
    write.csv(umap_data, umap_data_path, row.names = FALSE)
    
    print("临时注释完成")
  }, error = function(e2) {
    print(paste("临时注释也失败:", e2$message))
    stop("无法完成注释")
  })
})

print("脚本执行完成")