####Package####

library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(maftools)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(stringr)
library(ggpubr)
library(biomaRt)
library(AnnotationDbi)
library(dplyr)
library(tidyr)
library(ggplot2)
library(optparse)

# 定义命令行参数
option_list <- list(
  make_option(c("-w", "--wd"), type="character", default=NULL, 
              help="请输入工作路径"),
  make_option(c("-s", "--sd"), type="character", default=NULL, 
              help="请输入保存路径"),
  make_option(c("-g", "--gene"), type="character", default=NULL, 
              help="选择要观察的基因"),
  make_option(c("-t", "--time"), type="character", default=NULL, 
              help="请输入时间戳") # 新增时间戳参数
)

# 解析命令行参数
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 检查必要参数
if (is.null(opt$wd)) {
  stop("请提供文件路径 (-i/--wd)")
}
print("==================================================================================================================")
print("正在执行6.TCGA_mRNA_TIMER2.0.R")
####WorkingDirection####
# raw_WD <- '/public-SSD/home/sunpc/scseq/backend/R/TCGA'
####Loaddata####
g <- opt$gene
setwd(opt$sd)
gdata <- read.csv(paste0('TCGA/res/TCGAmrna_',g,"_tpm.csv"),row.names = 1)

setwd(opt$wd)
TCGAinfo <- read.csv('res/TCGAdir.csv',row.names = 1)
TCGAmrna <- read_rds('res/TCGAmrna_TPM.rds')

####Timer2.0####

list.files()
list.files('mapping')

timer <- read.csv('mapping/infiltration_estimation_for_tcga.csv')

gdata$Index <- sub('(TCGA-.*-.*-\\d+).*','\\1',gdata$Index)

gdata$du <- duplicated(gdata$Index)
table(gdata$du)
gdata <- gdata[gdata$du==F,]

gdata <- gdata[,-ncol(gdata)]
rownames(gdata) <- gdata$Index

rownames(timer) <- timer$cell_type

id <- intersect(rownames(timer),rownames(gdata))
lis <- colnames(timer)[grepl("Macrophage", colnames(timer))]

gdata <- gdata[id,c(1,2,4)]
mtimer <- timer[id,lis]

immdata <- cbind(gdata,mtimer)
immdata <- immdata[which(immdata$Tissue=='Tumor'),]

####correlation####

# 提取所有 cg 列名
lis
g <- colnames(immdata)[1]

# 定义一个函数来计算单个组的相关性
calculate_correlation <- function(subset_data, cg_columns,method,g) {
  results <- lapply(cg_columns, function(cg) {
    # 计算 RIOK2 和当前 cg 的相关性
    test_result <- cor.test(subset_data[[g]], subset_data[[cg]], method = method)#c('spearman','pearson','kendall)
    data.frame(
      cg = cg,
      Correlation = test_result$estimate,
      P_Value = test_result$p.value
    )
  })
  # 合并结果
  do.call(rbind, results)
}

# 按 Cancer 类型分组并计算相关性
correlation_results <- immdata %>% na.omit() %>%
  group_by(Cancer) %>%
  group_modify(~ calculate_correlation(.x, cg_columns=lis,method = 'pearson',g=g)) %>%
  ungroup()

# 查看结果
print(correlation_results)

correlation_results$res = case_when(
  correlation_results$P_Value < 0.001 ~ "***",
  correlation_results$P_Value < 0.01  ~ "**",
  correlation_results$P_Value < 0.05  ~ "*",
  is.na(correlation_results$P_Value)  ~ "",
  correlation_results$P_Value > 0.05 ~ ""
)


setwd(opt$sd)
if(!dir.exists("TCGA")){
  dir.create("TCGA")
}
setwd("TCGA")
# 如果需要保存为 CSV 文件
write.csv(correlation_results, paste0('res/',g,"_correlation_results_imm&TPM.csv"), row.names = FALSE)

#####plot#### 

correlation_results <- read.csv(paste0('res/',g,"_correlation_results_imm&TPM.csv"))
correlation_results$cg <- factor(correlation_results$cg,levels = c(
  'Macrophage_EPIC','Macrophage_TIMER','Macrophage_XCELL',
  'Macrophage.M0_CIBERSORT','Macrophage.M0_CIBERSORT.ABS',
  'Macrophage.M1_CIBERSORT','Macrophage.M1_CIBERSORT.ABS','Macrophage.M1_QUANTISEQ','Macrophage.M1_XCELL',
  'Macrophage.M2_CIBERSORT','Macrophage.M2_CIBERSORT.ABS','Macrophage.M2_QUANTISEQ','Macrophage.M2_XCELL',
  'Macrophage.Monocyte_MCPCOUNTER'
))

# 绘制热图
pdf(paste0("6.TCGA_mRNA_TIMER2.0_", opt$time, ".pdf"), width = 6, height = 6)
# pdf(file.path(opt$sd, "TCGA", paste0("6.TCGA_mRNA_TIMER2.0_", opt$time, ".pdf")), width = 6, height = 6)

p <- ggplot(correlation_results, aes(x = cg, y = Cancer, fill = Correlation)) +
  geom_tile(color = "white") +                        # 热图单元格
  geom_text(aes(label = res), color = "black", size = 3, fontface = "bold") +  # 在单元格中显示显著性标记，加粗
  scale_fill_gradient2(
    low = "#35704a", mid = "white", high = "#dea947", midpoint = 0,  # 颜色渐变
    limits = c(-0.6, 0.6), name = "Correlation"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),  # x 轴标签加粗
    axis.text.y = element_text(size = 8, face = "bold"),                        # y 轴标签加粗
    axis.title.x = element_text(face = "bold"),                                 # x 轴标题加粗
    axis.title.y = element_text(face = "bold"),                                 # y 轴标题加粗
    plot.title = element_text(hjust = 0.5, face = "bold"),                      # 图形标题加粗并居中
    legend.text = element_text(face = "bold"),                                  # 图例文本加粗
    legend.title = element_text(face = "bold")                                  # 图例标题加粗
  ) +
  labs(
    title = paste0('Heatmap of Correlation Between ',g," and Immune Score"),
    x = "Immune Score",
    y = "Cancer Type"
  )

print(p)
dev.off()
print("6.TCGA_mRNA_TIMER2.0.R绘图完成")
# 清理所有对象
rm(list = ls())
# 回收内存
gc()
