###cellchar pipeline

library(CellChat)
library(patchwork)
library(Seurat)
library(glue)
library(NMF)
options(stringsAsFactors = FALSE)


library(optparse)

# 定义命令行参数
option_list <- list(
  make_option(c("-w", "--wd"), type="character", default=NULL, 
              help="请输入工作路径"),
  make_option(c("-a", "--anno"), type="character", default=NULL, 
              help="输入的annotation路径"),
  make_option(c("-t", "--time"), type="character", default=NULL, 
              help="请输入时间戳")
)
# 解析命令行参数
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 检查必要参数
if (is.null(opt$wd)) {
  stop("请提供数据路径 (-i/--wd)")
}


caculate_CII<-function(work_path, annotation_path, timestamp){
	annotation_data = read.csv(annotation_path)  # 这个文件存储的是绘制umap图的数据，其annotation列存放的是最近一次的注释结果

	setwd(work_path)
	data <- readRDS("seurat_data.rds")
	
	print(data)	
	# 下面代码是直接从Seurat对象转换为CellChat对象
	# expression_matrix <- GetAssayData(data, layer  = "scale.data", assay = "RNA")
	# expression_matrix_sparse <- as(expression_matrix, "dgCMatrix")
	data$label_cc <- annotation_data$annotation # 将 Seurat 对象的细胞注释添加到 meta_data 中
	print("成功Seurat对象转换为CellChat对象")
	meta_data <- data@meta.data
	meta_data <- meta_data[,c("seurat_clusters","orig.ident","label_cc")]
	# data.MMM <- subset(data, layer  = "scale.data")
	# data.MMM <- JoinLayers(data.MMM[['RNA']])
	data.input <- data

	cellchat <- createCellChat(
		object = data.input,
		meta = meta_data,                                                                                     
		group.by = "label_cc"
	)
	print('成功创建CellChat对象')
	groupSize <- as.numeric(table(cellchat@idents)) # 每个亚群细胞数

	## 设置数据库
	CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
	# CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
	showDatabaseCategory(CellChatDB)
	CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # 使用数据库中的自分泌/旁分泌信号相互作用部分进行后续分析，可根据试验方案、目的选择
	# 使用数据库所有内容进行分析
	# CellChatDB.use <- CellChatDB
	# 在cellchat对象中设置使用的数据库
	cellchat@DB <- CellChatDB.use 

	# subset表达数据，提取仅在互作数据库中的基因，减少下游分析数据量
	cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
	# future::plan("multiprocess", workers = 4)
	# 识别在单个细胞类型中过表达配/受体
	cellchat <- identifyOverExpressedGenes(cellchat)
	print('成功识别在单个细胞类型中过表达配/受体')

	# 识别过表达互作对
	cellchat <- identifyOverExpressedInteractions(cellchat)
	print('成功识别过表达互作对')

	# 平滑表达值（目的是消除dropout影响，可选不用）
	# We also provide a function to project gene expression data onto protein-protein interaction (PPI) network. Specifically, a diffusion process is used to smooth genes’ expression values based on their neighbors’ defined in a high-confidence experimentally validated protein-protein network. This function is useful when analyzing single-cell data with shallow sequencing depth because the projection reduces the dropout effects of signaling genes, in particular for possible zero expression of subunits of ligands/receptors. 
	# cellchat <- projectData(cellchat, PPI.human)
	print('成功平滑表达值')

	# 互作可能性计算
	cellchat <- computeCommunProb(cellchat, raw.use = TRUE)    
	print('成功互作可能性计算')

	# 过滤表达细胞比例低的互作对
	cellchat <- filterCommunication(cellchat, min.cells = 10)
	print('成功过滤表达细胞比例低的互作对')

	# 提取关注得细胞间通讯关系
	df.net <- subsetCommunication(cellchat)
	print('成功提取关注得细胞间通讯关系')

	cellchat <- computeCommunProbPathway(cellchat)  # 推断通路水平的互作网络
	print('成功推断通路水平的互作网络')



	# 整合通讯网络结果
	# USER can also calculate the aggregated network among a subset of cell groups by setting sources.use and targets.use.
	cellchat <- aggregateNet(cellchat)
    # 确保当前工作目录正确
    print(paste("当前工作目录为:", getwd()))
    
    # 保存结果文件
    write.table(df.net, file = '00.CCI.net.res.tsv', quote = FALSE, sep = '\t')
    saveRDS(cellchat, file = '01.CCI.net.res.rds')
	# return  (cellchat)

	# 生成JSON文件
	annotations <- levels(cellchat@idents)  # 获取注释名
    groupSize <- as.numeric(table(cellchat@idents))  # 每个注释对应的数量
    
    # 对groupSize进行对数归一化
    normalizedGroupSize <- log1p(groupSize)  # 使用log1p避免log(0)问题
    
    # 合成nodes格式
    nodes <- lapply(seq_along(annotations), function(i) {
        list(
            id = annotations[i],
            name = annotations[i],
            category = annotations[i],
            symbolSize = normalizedGroupSize[i],  # 使用归一化后的值
            value = groupSize[i]  # 存入原始数据
        )
    })

    # 合成categories格式
    categories <- lapply(seq_along(annotations), function(i) {
        list(
            name = annotations[i]
        )
    })

    # 合成links格式
    net_count <- cellchat@net$count  # 获取矩阵
    links <- list()
    for (i in seq_len(nrow(net_count))) {
        for (j in seq_len(ncol(net_count))) {
            if (net_count[i, j] > 0) {  # 仅保留非零值
                links <- append(links, list(list(
                    source = annotations[i],
                    target = annotations[j],
                    value = net_count[i, j]
                )))
            }
        }
    }

    # 转换为JSON格式并保存
    library(jsonlite)
    json_data <- toJSON(list(nodes = nodes, categories = categories, links = links), pretty = TRUE, auto_unbox = TRUE)
    write(json_data, file = glue("count_{timestamp}.json"))
	# 合成links格式
    net_weight <- cellchat@net$weight  # 获取矩阵
    links <- list()
    for (i in seq_len(nrow(net_weight))) {
        for (j in seq_len(ncol(net_weight))) {
            if (net_weight[i, j] > 0) {  # 仅保留非零值
                links <- append(links, list(list(
                    source = annotations[i],
                    target = annotations[j],
                    value = net_weight[i, j] * 10
                )))
            }
        }
    }

    # 转换为JSON格式并保存
    json_data <- toJSON(list(nodes = nodes, categories = categories, links = links), pretty = TRUE, auto_unbox = TRUE)
    write(json_data, file = glue("weight_{timestamp}.json"))

	# 先检查可用的信号通路
	pathways <- cellchat@netP$pathways
	print("Available signaling pathways:")
	print(pathways)

	pathway.use<-"GALECTIN"
	print(paste0("分析通路: ", pathway.use))
		
		# 查看该通路的通信概率矩阵
		print("通路权重矩阵:")
		print(head(cellchat@net$weight))
		
		# 查看该通路的配体-受体对
		print("配体-受体对:")
		lr.pathway <- extractEnrichedLR(cellchat, signaling = pathway.use, geneLR.return = TRUE)
		print(lr.pathway)
		
		# 查看该通路的细胞通信数量
		print("细胞通信计数:")
		print(head(cellchat@net$count))
		
		# 查看通路的统计信息
		print("通路统计信息:")
		pathway.stats <- cellchat@netP$stats[pathway.use,]
		print(pathway.stats)
		
		# 绘图
		pdf_path <- glue("cellchat_heatmap_weight_{timestamp}.pdf")
		pdf(pdf_path, width = 8, height = 6)
		p1 <- netVisual_heatmap(cellchat, signaling = pathway.use, color.heatmap = "Reds", 
							measure = "weight", 
							title.name = paste0("MIF signaling ", "(weight)"))
		print(p1)
		dev.off()

		pdf_path <- glue("cellchat_heatmap_count_{timestamp}.pdf")
		pdf(pdf_path, width = 8, height = 6)
		p2 <- netVisual_heatmap(cellchat, signaling = pathway.use, color.heatmap = "Reds", 
							measure = "count",
							title.name = paste0("MIF signaling ", "(count)"))
		print(p2)
    dev.off()
	

	# 可以添加以下代码来查看数值差异
	# print("Weight matrix:")
	# print(cellchat@net$weight)
	# print("Count matrix:")
	# print(cellchat@net$count)
	
	# 清理所有对象
	rm(list = ls())
	# 回收内存
	gc()

	print("完成cellchat分析")
}

work_path = opt$wd
timestamp = opt$time
annotation_path = opt$anno
caculate_CII(work_path, annotation_path, timestamp)
