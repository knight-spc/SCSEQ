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
print("正在执行4.TCGA_mRNA_survival.R")

####WorkingDirection####
setwd(opt$sd)
g <- opt$gene
gdata <- read.csv(paste0('TCGA/res/TCGAmrna_',g,"_tpm.csv"),row.names = 1)
# raw_WD <- '/public-SSD/home/sunpc/scseq/backend/R/TCGA'
setwd(opt$wd)
####Loaddata####
TCGAinfo <- read.csv('res/TCGAdir.csv',row.names = 1)
TCGAmrna <- read_rds('res/TCGAmrna_TPM.rds')


TCGAsurv <- list()

for (i in 1:dim(TCGAinfo)[1]) {
  print(paste0("Accessing: ",TCGAinfo$Cancer[i]))
  if(TCGAinfo$Surv[i]!='none'){
    TCGAsurv[[i]] <- read.table(TCGAinfo$Surv[i], 
                                   header = TRUE, sep = "\t", 
                                   row.names = 1, check.names = FALSE)
  }
  else{
    TCGAsurv[[i]] <- data.frame()
    print('!!no Surv data')
  }
}

names(TCGAsurv) <- TCGAinfo$Cancer

saveRDS(TCGAsurv,file = 'res/TCGAsurv.rds')
####survival####

TCGAsurv <- readRDS('res/TCGAsurv.rds')

gdata$Patient <- sub('(TCGA-.*-.*-\\d+).*','\\1',gdata$Index)
table(duplicated(gdata$Patient))
table(sub('(TCGA-.*-.*-\\d+)(.*)','\\2',gdata$Index))
id <- duplicated(gdata$Patient)
gdata <- gdata[!id,]
rownames(gdata) <- sub('(TCGA-.*-.*-\\d+).*','\\1',gdata$Index)

gdata <- gdata[which(sub('TCGA-.*-.*-(\\d+).*','\\1',gdata$Index)!=11),]
table(sub('TCGA-.*-.*-(\\d+).*','\\1',gdata$Index))

survdata <- function(time,gdata){
  i <- intersect(rownames(time),rownames(gdata))
  h <- cbind(gdata[i,c(1,2)],time[i,])
  return(h)
}

survinfo <- list()

#clinical <- read.table(TCGAinfo$Clinical[3],header = TRUE, sep = "\t",row.names = 1, check.names = FALSE)
#clinical <- clinical[rownames(ex),]

#ex <- cbind(ex,clinical[,c('pathologic_stage','pathologic_T','pathologic_N','pathologic_M')])

#ex <- ex[which(!ex$pathologic_stage%in%c('','[Discrepancy]')),]

#table(ex$pathologic_stage)
#table(ex$pathologic_T)

for (i in 1:length(TCGAsurv)) {
  survinfo[[i]] <- survdata(time = TCGAsurv[[i]],gdata = gdata)
}

names(survinfo) <- names(TCGAsurv)

resultOS <- list()
resultDSS <- list()
resultPFI <- list()
resultDFI <- list()

for (i in 1:length(TCGAsurv)){
  print(names(TCGAsurv)[i])
  ex <- survinfo[[i]]

 colnames(ex)

 table(sub('TCGA-.*-.*-(\\d+)','\\1',rownames(ex)))

exOS <- na.omit(ex[,c(1:2,4:5)])
exDSS <- na.omit(ex[,c(1:2,6:7)])
exDFI <- na.omit(ex[,c(1:2,8:9)])
exPFI <- na.omit(ex[,c(1:2,10:11)])

# 使用数据框的第一列名称，而不是硬编码的 'HSPA1B'
gene_col_name <- colnames(ex)[1]

if(dim(exOS)[1]>1){
  cox_formula <- as.formula(paste("Surv(OS.time, OS) ~", gene_col_name))
  cox_model <- coxph(cox_formula, data = exOS)
  summary(cox_model)
  
  cox_summary <- summary(cox_model)
  result <- data.frame(
    Variable = rownames(cox_summary$coefficients),
    HR = round(cox_summary$coefficients[,2], 2),  # 风险比 (Hazard Ratio)
    Lower_CI = round(cox_summary$conf.int[,3], 2), # 95% 置信区间下界
    Upper_CI = round(cox_summary$conf.int[,4], 2), # 95% 置信区间上界
    p_value = round(cox_summary$coefficients[,5], 4)  # p 值
  )
  result$Cancer <- names(TCGAsurv)[i]
  resultOS[[i]] <- result}else{
    resultOS[[i]] <- data.frame(Variable = NA,HR = NA,Lower_CI=NA,Upper_CI=NA,p_value=NA,Cancer=NA)
  }

if(dim(exDSS)[1]>1){
  cox_formula <- as.formula(paste("Surv(DSS.time, DSS) ~", gene_col_name))
  cox_model <- coxph(cox_formula, data = exDSS)
  summary(cox_model)
  
  cox_summary <- summary(cox_model)
  result <- data.frame(
    Variable = rownames(cox_summary$coefficients),
    HR = round(cox_summary$coefficients[,2], 2),  # 风险比 (Hazard Ratio)
    Lower_CI = round(cox_summary$conf.int[,3], 2), # 95% 置信区间下界
    Upper_CI = round(cox_summary$conf.int[,4], 2), # 95% 置信区间上界
    p_value = round(cox_summary$coefficients[,5], 4)  # p 值
  )
  result$Cancer <- names(TCGAsurv)[i]
  resultDSS[[i]] <- result}else{
    resultDSS[[i]] <- data.frame(Variable = NA,HR = NA,Lower_CI=NA,Upper_CI=NA,p_value=NA,Cancer=NA)
  }

if(dim(exDFI)[1]>1){
  cox_formula <- as.formula(paste("Surv(DFI.time, DFI) ~", gene_col_name))
  cox_model <- coxph(cox_formula, data = exDFI)
  summary(cox_model)
  
  cox_summary <- summary(cox_model)
  result <- data.frame(
    Variable = rownames(cox_summary$coefficients),
    HR = round(cox_summary$coefficients[,2], 2),  # 风险比 (Hazard Ratio)
    Lower_CI = round(cox_summary$conf.int[,3], 2), # 95% 置信区间下界
    Upper_CI = round(cox_summary$conf.int[,4], 2), # 95% 置信区间上界
    p_value = round(cox_summary$coefficients[,5], 4)  # p 值
  )
  result$Cancer <- names(TCGAsurv)[i]
  resultDFI[[i]] <- result}else{
    resultDFI[[i]] <- data.frame(Variable = NA,HR = NA,Lower_CI=NA,Upper_CI=NA,p_value=NA,Cancer=NA)
  }

if(dim(exPFI)[1]>1){
  cox_formula <- as.formula(paste("Surv(PFI.time, PFI) ~", gene_col_name))
  cox_model <- coxph(cox_formula, data = exPFI)
  summary(cox_model)
  
  cox_summary <- summary(cox_model)
  result <- data.frame(
    Variable = rownames(cox_summary$coefficients),
    HR = round(cox_summary$coefficients[,2], 2),  # 风险比 (Hazard Ratio)
    Lower_CI = round(cox_summary$conf.int[,3], 2), # 95% 置信区间下界
    Upper_CI = round(cox_summary$conf.int[,4], 2), # 95% 置信区间上界
    p_value = round(cox_summary$coefficients[,5], 4)  # p 值
  )
  result$Cancer <- names(TCGAsurv)[i]
  resultPFI[[i]] <- result}else{
    resultPFI[[i]] <- data.frame(Variable = NA,HR = NA,Lower_CI=NA,Upper_CI=NA,p_value=NA,Cancer=NA)
  }
}

names(resultOS) <- TCGAinfo$Cancer
names(resultDSS) <- TCGAinfo$Cancer
names(resultDFI) <- TCGAinfo$Cancer
names(resultPFI) <- TCGAinfo$Cancer

# 处理四种生存分析结果
process_survival_data <- function(result_list, type) {
  a <- as.data.frame(result_list[[1]])
  for (i in 2:length(result_list)) {
    b <- as.data.frame(result_list[[i]])
    a <- rbind(a,b)
  }
  a <- na.omit(a)
  a$Type <- type
  return(a)
}

# 合并所有生存分析结果
all_results <- rbind(
  process_survival_data(resultOS, "OS"),
  process_survival_data(resultDSS, "DSS"),
  process_survival_data(resultDFI, "DFI"),
  process_survival_data(resultPFI, "PFI")
)


# 如果需要保存为 CSV 文件
write.csv(all_results, paste0('res/',g,"_TCGA_survival_results.csv"), row.names = FALSE)
# 为每种类型创建森林图
plot_forest <- function(data, type, ylim_max = 4) {
  # 为不同类型设置不同的颜色
  type_colors <- list(
    "OS" = "#FF7F00",    # 橙色
    "DSS" = "#006400",   # 深绿色
    "DFI" = "#800080",   # 紫色
    "PFI" = "#2364c2"    # 蓝色
  )
  
  data$Cancer <- factor(data$Cancer, levels = rev(unique(data$Cancer)))
  
  ggplot(data, aes(x = Cancer, y = HR, ymin = Lower_CI, ymax = Upper_CI)) +
    geom_pointrange(aes(color = Color), size = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "#C94052") +
    geom_text(aes(label = Significance, y = Upper_CI + 0.1, color = Color), size = 5) +
    scale_color_manual(values = c("grey" = "grey", "significant" = type_colors[[type]])) +
    coord_flip() +
    labs(
      title = paste0('Forest Plot for ', g, " Across TCGA (", type, ")"),
      x = "Cancer Type",
      y = paste0(type, " Hazard Ratio (HR)")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      legend.position = "none"
    ) +
    ylim(c(0, ylim_max))
}

# 保存所有森林图
setwd(opt$sd)
setwd("TCGA")
pdf(paste0("4.TCGA_mRNA_survival_", opt$time, ".pdf"), width = 24, height = 8)  # 修改宽高比例以适应1行4列

# 使用gridExtra来排列图形
library(gridExtra)
plot_list <- list()

# 绘制四种生存分析的森林图
for (type in c("OS", "DSS", "DFI", "PFI")) {
  data_subset <- subset(all_results, Type == type)
  data_subset <- data_subset[-11,] # 移除特定行，如果需要的话
  data_subset <- data_subset %>%
    mutate(
      Significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ ""
      ),
      Color = case_when(
        p_value < 0.05 ~ "significant",
        TRUE ~ "grey"
      )
    )
  plot_list[[type]] <- plot_forest(data_subset, type)
}

# 使用grid.arrange将图形排列为1行4列
grid.arrange(grobs = plot_list, ncol = 4)

dev.off()

print("4.TCGA_mRNA_survival.R生存分析图像绘制完成")
# 清理所有对象
rm(list = ls())
# 回收内存
gc()
