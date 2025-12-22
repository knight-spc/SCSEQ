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
print("正在执行2.TCGA_mRNA_boxplot.R")
####WorkingDirection####
# raw_WD <- '/public-SSD/home/sunpc/scseq/backend/R/TCGA'
setwd(opt$wd)
####Loaddata####

TCGAinfo <- read.csv('res/TCGAdir.csv',row.names = 1)
genelength <- read.table('mapping/gene_length_Table.txt')

TCGAmrna <- list()

for (i in 1:dim(TCGAinfo)[1]) {
  print(paste0("Accessing: ",TCGAinfo$Cancer[i]))
  if(TCGAinfo$mRNA[i]!='none'){
    TCGAmrna[[i]] <- read.table(TCGAinfo$mRNA[i], 
                                      header = TRUE, sep = "\t", 
                                      row.names = 1, check.names = FALSE)
  }
  else{
    TCGAmrna[[i]] <- data.frame()
    print("!!No mrna data")
  }
  
}

names(TCGAmrna) <- TCGAinfo$Cancer

saveRDS(TCGAmrna,file = 'res/TCGAmrna.rds')

####Tpm Gene####


TCGAtpm <- list()

for (i in 1:dim(TCGAinfo)[1]) {
  print(paste0("Accessing: ",TCGAinfo$Cancer[i]))
  if(TCGAinfo$mRNAtpm[i]!='none'){
    
    data <- read.table(TCGAinfo$mRNAtpm[i], sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    data$Ensembl_ID <- substr(data$Ensembl_ID,1,15)
    data <- data[which(data$Ensembl_ID%in%genelength$simple),]
    colnames(data) <- gsub('\\.','-',colnames(data))
    du <- duplicated(data$Ensembl_ID)
    data <- data[!(du),]
    rownames(data) <- data$Ensembl_ID
    
    
    TCGAtpm[[i]] <- data
  }
  else{
    TCGAtpm[[i]] <- data.frame()
    print("!!No mrna TPM data")
  }
}

names(TCGAtpm) <- TCGAinfo$Cancer

du <- genelength$simple

for (i in 1:dim(TCGAinfo)[1]) {
  du <- intersect(du,rownames(TCGAtpm[[i]]))
}

genelength <- genelength[du,]

TCGAtpm0 <- TCGAtpm

for (i in 1:dim(TCGAinfo)[1]) {
  data <- TCGAtpm0[[i]][du,]
  data$Ensembl_ID <- genelength$genename 
  id <- duplicated(data$Ensembl_ID)
  data <- data[!id,]
  rownames(data) <- data$Ensembl_ID
  TCGAtpm[[i]] <- data[,-1]
}

saveRDS(TCGAtpm,file = 'res/TCGAmrna_TPM.rds')



####Counts Gene####


TCGAcount <- list()

for (i in 1:dim(TCGAinfo)[1]) {
  print(paste0("Accessing: ",TCGAinfo$Cancer[i]))
  if(TCGAinfo$mRNAcount[i]!='none'){
    
    data <- read.table(TCGAinfo$mRNAcount[i], sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    data$Ensembl_ID <- substr(data$Ensembl_ID,1,15)
    data <- data[which(data$Ensembl_ID%in%genelength$simple),]
    colnames(data) <- gsub('\\.','-',colnames(data))
    du <- duplicated(data$Ensembl_ID)
    data <- data[!(du),]
    rownames(data) <- data$Ensembl_ID
    
    
    TCGAcount[[i]] <- data
  }
  else{
    TCGAcount[[i]] <- data.frame()
    print("!!No mrna count data")
  }
}

names(TCGAcount) <- TCGAinfo$Cancer

du <- genelength$simple

for (i in 1:dim(TCGAinfo)[1]) {
  du <- intersect(du,rownames(TCGAcount[[i]]))
}

genelength <- genelength[du,]

TCGAcount0 <- TCGAcount

for (i in 1:dim(TCGAinfo)[1]) {
  data <- TCGAcount0[[i]][du,]
  data$Ensembl_ID <- genelength$genename 
  id <- duplicated(data$Ensembl_ID)
  data <- data[!id,]
  rownames(data) <- data$Ensembl_ID
  TCGAcount[[i]] <- data[,-1]
}

saveRDS(TCGAcount,file = 'res/TCGAmrna_count.rds')
####Interesting Gene####
#####Boxplot####

TCGAmrna <- readRDS('res/TCGAmrna_TPM.rds')

#set gene name
# g <- "HSPA1B"#'HSPA1B'#
g <- opt$gene#'HSPA1B'#

#check its position
for (i in 1:dim(TCGAinfo)[1]) {
  print(paste0("Accessing: ",TCGAinfo$Cancer[i]))
  print(dim(TCGAmrna[[i]]))
  print(paste0('Target gene (',g,') Found ',(g%in%rownames(TCGAmrna[[i]]))))
}

#biomarker function
extractgene <- function(g,data,name){
  i <- as.data.frame(t(data[g,]))
  i[,2] <- rep(name,dim(i)[1])
  i[,3] <- rownames(i)
  colnames(i)[c(2,3)] <- c("Cancer",'Index') 
  rownames(i) <- 1:dim(i)[1]
  return(i)
}

#extract and merge

gData <- extractgene(g,TCGAmrna[[1]],names(TCGAmrna)[1])

for (i in 2:length(TCGAmrna)) {
  gData <- rbind(gData,extractgene(g,TCGAmrna[[i]],names(TCGAmrna)[i]))
}

gData$Tissue <- ifelse(sub('TCGA-.*-.*-(\\d+).*','\\1',gData$Index)=='11','Normal','Tumor')
gData$Tissue <- factor(gData$Tissue,levels=c("Tumor","Normal"))

table(gData$Cancer,gData$Tissue)

summary(gData[,1])
if(!dir.exists('TCGA')){
  dir.create('TCGA')
}

setwd(opt$sd)
if(!dir.exists('TCGA')){
  dir.create('TCGA')
}
setwd("TCGA")

if(!dir.exists('plot')){
  dir.create('plot')
}
if(!dir.exists('res')){
  dir.create('res')
}

write.csv(gData,file = paste0('res/','TCGAmrna_',g,'_tpm','.csv'),row.names = T)

#####Plotting####

#pdf(file = paste0('plot/boxplot_',g,'.pdf'),width = 7,height = 4)

boxdata <- gData

#####T test####

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
pdf(paste0("2.TCGA_mRNA_boxplot_", opt$time, ".pdf"), width = 6, height = 6)

p <- ggplot(boxdata, aes(x = Cancer, y = .data[[colnames(boxdata)[1]]], color = Tissue)) +
  geom_boxplot(
    fill = NA, outlier.shape = 16, outlier.size = 3, size = 1.2, width = 0.6
  ) +
  scale_color_manual(values = c("Tumor" = "#C94052", "Normal" = "#4F7A5E")) +
  labs(title = "", y = paste0(g," expression"), x = "") +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(color = "black", face = "bold", size = 14),
    axis.line = element_line(color = "black", size = 1.2),
    plot.title = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(face = "bold", size = 12),
    legend.position = "top"
  ) +
  scale_y_continuous(limits = c(0, 7.5)) +#c(5, 12)
  geom_text(
    data = p_values, 
    aes(x = Cancer, y = 11.5, label = res),  # 在 y=11.5 处添加显著性标记
    inherit.aes = FALSE,
    size = 6, fontface = "bold"
  )



# # 动态计算显著性标记的Y轴位置
# p_values <- p_values %>%
#   group_by(Cancer) %>%
#   mutate(y_pos = max(boxdata[[colnames(boxdata)[1]]][boxdata$Cancer == Cancer], na.rm = TRUE) + 0.5)

# p <- ggplot(boxdata, aes(x = Cancer, y = .data[[colnames(boxdata)[1]]], color = Tissue)) +
#   geom_boxplot(
#     fill = NA, 
#     outlier.shape = 16, 
#     outlier.size = 1.5,  # 减小离群点大小
#     size = 0.8,          # 减小边框粗细
#     width = 0.5          # 减小箱线图宽度
#   ) +
#   scale_color_manual(values = c("Tumor" = "#C94052", "Normal" = "#4F7A5E")) +
#   labs(title = "", y = paste0(g, " expression"), x = "") +
#   theme_minimal(base_size = 15) +
#   theme(
#     panel.grid = element_blank(),
#     axis.text = element_text(color = "black", face = "bold", size = 14),
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10), # 调整标签角度和大小
#     axis.title = element_text(color = "black", face = "bold", size = 14),
#     axis.line = element_line(color = "black", size = 1.0), # 减小轴线粗细
#     plot.title = element_text(face = "bold", size = 16),
#     legend.title = element_text(face = "bold", size = 14),
#     legend.text = element_text(face = "bold", size = 12),
#     legend.position = "top"
#   ) +
#   scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + # 动态缓冲
#   geom_text(
#     data = p_values, 
#     aes(x = Cancer, y = y_pos, label = res), # 动态Y轴位置
#     inherit.aes = FALSE,
#     size = 4, #
#     fontface = "bold"
#   )


print(p)
dev.off()
print("2.TCGA_mRNA_boxplot.R绘图完成")

# 清理所有对象
rm(list = ls())
# 回收内存
gc()

