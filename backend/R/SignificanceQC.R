print(.libPaths())
library(Seurat) # 4.2.0
library(RColorBrewer) # 1.1-3
library(dplyr) # 1.0.10
library(harmony) # 1.2.0
library(pals) # 1.7
library(ggplot2) # 3.3.6
library(glue)


library(tidyverse)
library(viridis)
library(gridExtra)
library(SeuratData)
library(ggpubr)



get_significance_QC <- function(wd_path, cellid, gene) {
    setwd(wd_path)
    data <- readRDS("seurat_data.rds")
    print('111111111')
    cellid <- list(unlist(cellid))
    print(str(cellid))
    data.markers = read.csv( 'all_marker.csv')
    top_markers <-data.markers %>%
            group_by(cluster) %>%
            top_n(n = 20, wt = avg_log2FC)
    # celltype <- unique(as.character(data$celltype))
    # my_comparisons <- lapply(celltype[-cellid], function(x)
    #                 c(x,celltype[cellid]))

    pdf("QC_comparisons.pdf", width = 6, height = 6)
    p <- VlnPlot(data, features = gene ,add.noise = F,pt.size = 0.1)&
            theme_bw()& theme(axis.title.x = element_blank(),
                axis.text.x = element_text(color ='black',face ="bold", size = 12,angle = 45),
                axis.text.y = element_text(color ='black', face ="bold"),
                axis.title.y = element_text(color ='black', face ="bold", size = 15),  
                panel.grid.major = element_blank(),  
                panel.grid.minor = element_blank(),   
                panel.border = element_rect(color="black",size = 1.2, linetype="solid"),  
                panel.spacing = unit(0.12,"cm"),   
                plot.title = element_text(hjust = 0.5, face ="bold.italic"),  
                legend.position ='none')&  
            stat_compare_means(method="t.test",
                hide.ns = F, 
                comparisons = cellid,
                label="p.signif",
                bracket.size=0.8,
                tip.length=0,
                size=6)& 
            scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
    print(p)
    dev.off()
    # 清理所有对象
    rm(list = ls())
    # 回收内存
    gc()
}
