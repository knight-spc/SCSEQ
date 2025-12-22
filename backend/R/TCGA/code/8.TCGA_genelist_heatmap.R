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
library(reshape2)
library(optparse)

# 设置命令行参数选项
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

####Loaddata####
setwd(opt$wd)

TCGAinfo <- read.csv('res/TCGAdir.csv',row.names = 1)
TCGAmrna <- readRDS('res/TCGAmrna_TPM.rds')
#####surv####  
TCGAsurv <- readRDS('res/TCGAsurv.rds')  
####Interesting Gene list####
if (!is.null(opt$gene)) {
    # 将输入的字符串按逗号分割并去除空格
    index <- strsplit(opt$gene, ",")[[1]]
    index <- trimws(index)  # 去除空格
}
# index <- c('AKT1', 'CCND1', 'EGFR','FOS','HIF1A', 'MAPK3', 'MDM2', 'MMP9',  'NFKB1', 'STAT3')
print(index)
length(index)

#####heatmap####

#biomarker function
extractgene <- function(g,data,name){
  i <- as.data.frame(t(data[g,]))  # 转置数据
  i[,2] <- rep(name,dim(i)[1])     # 添加癌症类型标签
  i[,3] <- rownames(i)             # 添加样本ID
  colnames(i)[c(2,3)] <- c("Cancer",'Index') 
  rownames(i) <- 1:dim(i)[1]
  return(i)
}

# 检查每个癌症类型中目标基因的存在情况
for (i in 1:dim(TCGAinfo)[1]) {
  print(paste0("Accessing: ",TCGAinfo$Cancer[i]))
  print(dim(TCGAmrna[[i]]))
  print(paste0('Target gene (',index,') Found ',(index%in%rownames(TCGAmrna[[i]]))))
}

# 创建列表存储基因表达和统计检验结果
indexlist <- list()   # 存储基因表达数据
pvaluelist <- list()  # 存储统计检验结果

# 对每个基因进行分析
for (j in 1:length(index)){
  #set gene name
  g <- index[j]#'HSPA1B'#

  #extract and merge

  gData <- extractgene(g,TCGAmrna[[1]],names(TCGAmrna)[1])
  for (i in 2:length(TCGAmrna)) {
    gData <- rbind(gData,extractgene(g,TCGAmrna[[i]],names(TCGAmrna)[i]))
  }
  
  # 区分肿瘤和正常样本
  gData$Tissue <- ifelse(sub('TCGA-.*-.*-(\\d+).*','\\1',gData$Index)=='11','Normal','Tumor')
  gData$Tissue <- factor(gData$Tissue,levels=c("Tumor","Normal"))

  table(gData$Cancer,gData$Tissue)

  summary(gData[,1])
  colnames(gData)[1] <- 'Gene'
  indexlist[[j]] <- gData
  names(indexlist)[j] <- index[j]

#pdf(file = paste0('plot/boxplot_',g,'.pdf'),width = 7,height = 4)

  boxdata <- gData


#####T test####

# 计算 p 值
  p_values <- boxdata %>%
    group_by(Cancer, Tissue) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = Tissue, values_from = n, values_fill = 0) %>%
    mutate(p = NA)  # 初始化p值列

for (i in 1:dim(p_values)[1]) {
  if(p_values$Tumor[i]>=2 & p_values$Normal[i]>=2){
    p_values$p[i] <- round(t.test(Gene ~ Tissue,
                                 data = boxdata[which(boxdata$Cancer%in%p_values$Cancer[i]),]
    )$p.value, 2)
  }
}

p_values$res = case_when(
  p_values$p < 0.001 ~ "***",
  p_values$p < 0.01  ~ "**",
  p_values$p < 0.05  ~ "*",
  is.na(p_values$p)  ~ "",
  p_values$p > 0.05 ~ ""
)

pvaluelist[[j]] <- p_values
names(pvaluelist)[j] <- index[j]
}

####merge####

a <- indexlist[[1]][which(indexlist[[1]]$Tissue=='Tumor'),c('Gene',"Cancer")]
colnames(a)[1] <- 'Expression'
# 按 Cancer 和 Genename 分组，计算 Expression 的均值
mean<- a %>%
  group_by(Cancer) %>% # 按 Cancer 和 Genename 分组
  summarise(Expression = mean(Expression, na.rm = TRUE)) %>% # 计算均值
  ungroup() # 取消分组
mean$Genename <- index[1]

for (i in 2:length(indexlist)) {
  b <- indexlist[[i]][which(indexlist[[i]]$Tissue=='Tumor'),c('Gene',"Cancer")]
  colnames(b)[1] <- 'Expression'
  # 按 Cancer 和 Genename 分组，计算 Expression 的均值
  c<- b %>%
    group_by(Cancer) %>% # 按 Cancer 和 Genename 分组
    summarise(Expression = mean(Expression, na.rm = TRUE)) %>% # 计算均值
    ungroup() # 取消分组
  c$Genename <- index[i]
  mean <- rbind(mean,c)
}

setwd(opt$sd)
setwd("TCGA")
# 创建必要的文件夹
dir.create("plot", showWarnings = FALSE)
dir.create("plot/forest_list", recursive = TRUE, showWarnings = FALSE)
dir.create("res", showWarnings = FALSE)
#####heatmap####
# 将数据重塑为宽格式
mean_wide <- mean %>%
  pivot_wider(names_from = Cancer, 
              values_from = Expression) %>%
  column_to_rownames("Genename")

# 只对表达值进行标准化
mean_scaled <- t(scale(t(mean_wide)))

# 转回长格式用于绘图
mean_long <- mean_scaled %>%
  as.data.frame() %>%
  rownames_to_column("Genename") %>%
  pivot_longer(-Genename, 
               names_to = "Cancer", 
               values_to = "Expression")

# 绘制热图
p <- ggplot(mean_long, aes(x = Cancer, y = Genename, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#35704a", mid = "white", high = "#dea947", 
    midpoint = 0,  # 标准化后的中点为0
    name = "Scaled\nExpression"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  labs(
    title = "Pan-cancer Expression of Target Gene",
    x = "Cancer Type",
    y = "Gene"
  )

# 显示热图
pdf(paste0('08.heatmap_', opt$time, '.pdf'), width = 8, height = 8)
print(p)
dev.off()


#####boxplot####


pdf(file = 'plot/boxplot_list.pdf',width = 9,height = 6)

for (i in 1:length(indexlist)){
  
  g <- index[i]
  boxdata <- indexlist[[i]]
  
  b <- ggplot(boxdata, aes(x = Cancer, y = Gene, color = Tissue)) +#y = RIOK2
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
  scale_y_continuous(limits = c(0, max(boxdata$Gene)+1)) +#c(5, 12)
  geom_text(
    data = p_values, 
    aes(x = Cancer, y = max(boxdata$Gene)+0.5, label = res),  # 在 y=11.5 处添加显著性标记
    inherit.aes = FALSE,
    size = 6, fontface = "bold"
  )
  print(b)
}

dev.off()

#####pvalue####

resv <- pvaluelist[[1]]
resv$genename <- index[1]
  
for (i in 2:length(pvaluelist)) {
  a <- pvaluelist[[i]]
  a$genename <- index[i]
  resv <- rbind(resv,a)
}  
  
write.csv(resv,file = 'res/pvalue_list.csv',row.names = T)  
  



for (sh in 1:length(index)) {
  
  shn <- index[sh]
  print(shn)
  gdata <- indexlist[[sh]]
  
  rownames(gdata) <- gdata$Index
  
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
  
  for (i in 1:length(TCGAsurv)) {
    survinfo[[i]] <- survdata(time = TCGAsurv[[i]],gdata = gdata)
  }
  
  names(survinfo) <- names(TCGAsurv)
  
  resultOS  <- list()
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
    
    if(dim(exOS)[1]>1){
      cox_model <- coxph(Surv(OS.time, OS) ~ Gene 
                         , data = exOS)
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
      cox_model <- coxph(Surv(DSS.time, DSS) ~ Gene 
                         , data = exDSS)
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
      cox_model <- coxph(Surv(DFI.time, DFI) ~ Gene 
                         , data = exDFI)
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
      cox_model <- coxph(Surv(PFI.time, PFI) ~ Gene 
                         , data = exPFI)
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
  
  names(resultOS) <-  TCGAinfo$Cancer
  names(resultDSS) <- TCGAinfo$Cancer
  names(resultDFI) <- TCGAinfo$Cancer
  names(resultPFI) <- TCGAinfo$Cancer
  
  #OS
  a <- as.data.frame(resultOS[[1]])
  for (i in 2:length(resultOS)) {
  b <- as.data.frame(resultOS[[i]])
  a <- rbind(a,b)
  }
  a <- na.omit(a)
  dataOS <- a
  
  #DSS
  a <- as.data.frame(resultDSS[[1]])
  for (i in 2:length(resultOS)) {
    b <- as.data.frame(resultDSS[[i]])
    a <- rbind(a,b)
  }
  a <- na.omit(a)
  dataDSS <- a
  
  #DFI
  a <- as.data.frame(resultDFI[[1]])
  for (i in 2:length(resultOS)) {
    b <- as.data.frame(resultDFI[[i]])
    a <- rbind(a,b)
  }
  a <- na.omit(a)
  dataDFI <- a
  
  #PFI
  a <- as.data.frame(resultPFI[[1]])
  for (i in 2:length(resultOS)) {
    b <- as.data.frame(resultPFI[[i]])
    a <- rbind(a,b)
  }
  a <- na.omit(a)
  dataPFI <- a
  datalist <- list(dataOS,dataDSS,dataDFI,dataPFI)
  names(datalist) <- c('OS','DSS','DFI','PFI')
  
  # 修改PDF输出部分
  pdf(height = 8, width = 24, paste0('plot/forest_list_',opt$time, '/', shn, '.pdf'))  # 修改尺寸以适应并排显示
  
  # 创建图形列表
  plot_list <- list()
  
  # 设置每种类型的颜色
  type_colors <- c(
    'OS' = '#ff7e44',
    'DSS' = '#2d534b',
    'DFI' = '#86005a',
    'PFI' = '#2364c2'
  )
  
  for (k in 1:length(datalist)) {
    data <- datalist[[k]] %>%
      mutate(
        Significance = case_when(
          p_value < 0.001 ~ "***",
          p_value < 0.01 ~ "**",
          p_value < 0.05 ~ "*",
          TRUE ~ ""
        ),
        Color = case_when(
          p_value < 0.05 ~ type_colors[k],
          TRUE ~ "grey"
        )
      )
    
    # 反转 Cancer 的显示顺序
    data$Cancer <- factor(data$Cancer, levels = rev(unique(data$Cancer)))
    
    # 绘制森林图
    p <- ggplot(data, aes(x = Cancer, y = HR, ymin = Lower_CI, ymax = Upper_CI)) +
      geom_pointrange(aes(color = Color), size = 0.8) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "#C94052") +
      geom_text(aes(label = Significance, y = Upper_CI + 0.1, color = Color), size = 5) +
      scale_color_identity() +
      coord_flip() +
      labs(
        title = paste0(names(datalist)[k], " Forest Plot for ", shn),
        x = "Cancer Type",
        y = paste0(names(datalist)[k], " Hazard Ratio (HR)")
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
      ylim(c(0, max(data$Upper_CI[which(data$p_value<0.05)]) + 0.5))
    
    plot_list[[k]] <- p
  }
  
  # 使用 gridExtra 包将四张图并排显示
  library(gridExtra)
  grid.arrange(grobs = plot_list, ncol = 4)
  
  dev.off()
  
}

print("完成分析")









