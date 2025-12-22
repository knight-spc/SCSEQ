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
print("正在执行5.TCGA_mRNA_meth.R")
####WorkingDirection####
setwd(opt$sd)
g <- opt$gene
gdata <- read.csv(paste0('TCGA/res/TCGAmrna_',g,"_tpm.csv"),row.names = 1)
# raw_WD <- '/public-SSD/home/sunpc/scseq/backend/R/TCGA'
setwd(opt$wd)
####Loaddata####
TCGAinfo <- read.csv('res/TCGAdir.csv',row.names = 1)
TCGAmrna <- read_rds('res/TCGAmrna_TPM.rds')


####Methdata####

list.files()
list.files('mapping')

methinfo <- read.table('mapping/probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy.txt', sep = "\t", header = TRUE, stringsAsFactors = FALSE)

TCGAmeth <- list()

for (i in 1:dim(TCGAinfo)[1]) {
  print(paste0(i,".Accessing: ",TCGAinfo$Cancer[i]))
  if(TCGAinfo$Methylation[i]!='none'){
    data <- read.table(TCGAinfo$Methylation[i], sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    colnames(data) <- gsub('\\.','-',colnames(data))
    TCGAmeth[[i]] <- data
  }
  else{
    TCGAmeth[[i]] <- data.frame()
    print("!!No mrna Methylation data")
  }
}

names(TCGAmeth) <- TCGAinfo$Cancer

####intersting gene####
g <- opt$gene

id <- methinfo$id[grepl(g, methinfo$gene)]

saveRDS(TCGAmeth,file = 'res/TCGAmeth.rds')

gmeth <- list()

for (i in 1:dim(TCGAinfo)[1]) {
  ex <- TCGAmeth[[i]]
  data <- ex[which(ex[,1]%in%id),]
  gmeth[[i]] <- data
}

names(gmeth) <- TCGAinfo$Cancer

#extract and merge

gData <- gmeth[[1]]
rownames(gData) <- gData$sample 
gData <- gData[,-1]
gData <- as.data.frame(t(gData))
gData$Cancer <- rep(names(gmeth)[1],dim(gData)[1])
gData$Index <- rownames(gData)

for (i in 2:length(gmeth)) {
  a <- gmeth[[i]]
  rownames(a) <- a$sample 
  a <- a[,-1]
  a <- as.data.frame(t(a))
  a$Cancer <- rep(names(gmeth)[i],dim(a)[1])
  a$Index <- rownames(a)
  gData <- rbind(gData,a)
}

gData$Tissue <- ifelse(sub('TCGA-.*-.*-(\\d+).*','\\1',gData$Index)=='11','Normal','Tumor')
gData$Tissue <- factor(gData$Tissue,levels=c("Tumor","Normal"))

table(gData$Cancer,gData$Tissue)

summary(gData[,1])


setwd(opt$sd)
setwd("TCGA")
if(!dir.exists('plot')){
  dir.create('plot')
}

write.csv(gData,file = paste0('res/','TCGAmeth_',g,'.csv'),row.names = T)

####merge####

gdata$Index <- sub('(TCGA-.*-.*-\\d+).*','\\1',gdata$Index)

gdata$du <- duplicated(gdata$Index)
table(gdata$du)
gdata <- gdata[gdata$du==F,]

gdata <- gdata[,-ncol(gdata)]
rownames(gdata) <- gdata$Index

id <- intersect(rownames(gdata),rownames(gData))
# 只选择甲基化位点列
meth_cols <- grep("^cg", colnames(gData), value = TRUE)

if(length(meth_cols) >= 11) {
    meth_cols <- meth_cols[1:11]
}
print(meth_cols)

gData_subset <- gData[id, meth_cols]
mdata <- cbind(gdata[id,], gData_subset)

mdata <- mdata[which(mdata$Tissue=="Tumor"),]
write.csv(mdata,file = paste0('res/','TCGAmethTPM_',g,'.csv'),row.names = T)

#####correlation####

# 加载必要的包
library(dplyr)

data <- mdata
# 数据准备
# 假设您的数据框名为 data
# data 包含列: HSPA1B, Cancer, cg12346467, cg12676081, ...

# 提取所有 cg 列名
cg_columns <- colnames(data)[grepl("^cg", colnames(data))]

# 定义一个函数来计算单个组的相关性
calculate_correlation <- function(subset_data, cg_columns) {
  results <- lapply(cg_columns, function(cg) {
    # 计算 HSPA1B 和当前 cg 的相关性
    test_result <- cor.test(subset_data[[g]], subset_data[[cg]], method = "pearson")  # 使用g变量作为列名
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
correlation_results <- data %>%
  group_by(Cancer) %>%
  group_modify(~ calculate_correlation(.x, cg_columns)) %>%
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



# 如果需要保存为 CSV 文件
write.csv(correlation_results, paste0('res/',g,"_correlation_results_imm&TPM.csv"), row.names = FALSE)

#####plot#### 

# correlation_results <- read.csv("res/HSPA1B_correlation_results_meth&TPM.csv")

# 绘制热图
pdf(paste0("5.TCGA_mRNA_meth_", opt$time, ".pdf"), width = 10, height = 10)

p <- ggplot(correlation_results, aes(x = Cancer, y = cg, fill = Correlation)) +
  geom_tile(color = "white") +                        # 热图单元格
  geom_text(aes(label = res), color = "black", size = 3, fontface = "bold") +  # 在单元格中显示显著性标记，加粗
  scale_fill_gradient2(
    low = "#2364c2", mid = "white", high = "#C94052", midpoint = 0,  # 颜色渐变
    limits = c(-0.7, 0.5), name = "Correlation"
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
    title = paste0('Heatmap of Correlation Between ',g," and CpG Sites"),
    x = "Cancer Type",
    y = "CpG Sites"
  )




print(p)
dev.off()
print("5.TCGA_mRNA_meth.R绘图完成")
# 清理所有对象
rm(list = ls())
# 回收内存
gc()


