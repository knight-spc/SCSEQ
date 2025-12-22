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
library(optparse)

# 定义命令行参数
option_list <- list(

  make_option(c("-w", "--wd"), type="character", default=NULL, 
              help="请输入保存路径"),
  make_option(c("-g", "--gene"), type="character", default=NULL, 
              help="选择要观察的基因"),
  make_option(c("-i", "--file"), type="character", default=NULL, 
              help="请输入文件路径"),
  make_option(c("-o", "--outlier_size"), type="character", default=NULL, 
              help="请输入保存路径"),
  make_option(c("-s", "--frame_size"), type="character", default=NULL, 
              help="选择要观察的基因"),
  make_option(c("-b", "--box_width"), type="character", default=NULL, 
              help="请输入箱线宽度"),
  make_option(c("-a", "--axis_textsize"), type="character", default=NULL, 
              help="请输入标签字体大小"),
  make_option(c("-z", "--base_size"), type="character", default=NULL, 
              help="请输入基础字体大小"),
  make_option(c("-y", "--ylimits"), type="character", default=NULL, 
              help="请输入y轴范围")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 检查必要参数
if (is.null(opt$wd)) {
  stop("请提供文件路径 (-w/--wd)")
}
setwd(opt$wd)

print("==================================================================================================================")
print("正在执行resize.R")

g <- opt$gene#'HSPA1B'#

boxdata <- read.csv(paste0('TCGA/res/TCGAmrna_',g,"_tpm.csv"))
boxdata <- boxdata[, -1]  # 删除第一列（假设它是行名列）
# 计算 p 值
p_values <- boxdata %>%
  group_by(Cancer, Tissue) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Tissue, values_from = n, values_fill = 0) 

for (i in 1:dim(p_values)[1]) {
  
  if(p_values$Tumor[i]>=2&p_values$Normal[i]>=2){
    p_values$p[i] <- round(t.test(boxdata[which(boxdata$Cancer%in%p_values$Cancer[i]),1] ~ boxdata$Tissue[which(boxdata$Cancer%in%p_values$Cancer[i])],
                            data = boxdata[which(boxdata$Cancer%in%p_values$Cancer[i]),]
                            )$p.value,2)
  }else{
    p_values$p[i] <- NA
  }
}


p_values$res = case_when(
  p_values$p < 0.001 ~ "***",
  p_values$p < 0.01  ~ "**",
  p_values$p < 0.05  ~ "*",
  is.na(p_values$p)  ~ "",
  p_values$p > 0.05 ~ ""
)

# 画图
pdf(opt$file, width = 6, height = 6)


p <- ggplot(boxdata, aes(x = Cancer, y = .data[[colnames(boxdata)[1]]], color = Tissue)) +
  geom_boxplot(
    fill = NA, outlier.shape = 16, outlier.size = as.numeric(opt$outlier_size), size = as.numeric(opt$frame_size), width = as.numeric(opt$box_width) ) +#width箱线宽度
  scale_color_manual(values = c("Tumor" = "#C94052", "Normal" = "#4F7A5E")) +
  labs(title = "", y = paste0(g," expression"), x = "") +
  theme_minimal(base_size = as.numeric(opt$base_size)) +#字体基础大小
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", face = "bold", size = as.numeric(opt$axis_textsize)),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(color = "black", face = "bold", size = 14),
    axis.line = element_line(color = "black", size = 1.2),
    plot.title = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(face = "bold", size = 12),
    legend.position = "top"
  ) +
  scale_y_continuous(limits = c(0, as.numeric(opt$ylimits))) +#y轴范围
  geom_text(
    data = p_values, 
    aes(x = Cancer, y = 11.5, label = res),  # 在 y=11.5 处添加显著性标记
    inherit.aes = FALSE,
    size = 6, fontface = "bold"
  )
print(p)
dev.off()
print("2.TCGA_mRNA_boxplot.R绘图完成")



# 清理所有对象
rm(list = ls())
# 回收内存
gc()

