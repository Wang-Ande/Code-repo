本项目以adjusted_gene_count_matrix_merge数据集进行DESeq2差异分析过程
实验分组：MOLM13_WT（3例），MOLM13_D（5例）
R版本：
R包：tidyverse，lance，DESeq2	
# 设置工作路径 ----
setwd('./03_result') # 设置工作路径
if(!dir.exists('./DE')){
  dir.create('./DE')
} # 判断该工作路径下是否存在名为DE的文件夹，如果不存在则创建，如果存在则pass
setwd('./DE/') # 设置路径到刚才新建的DE下
# 加载包 ----
library(tidyverse)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
}
library(DESeq2)
library(openxlsx)

# count数据导入----
dat <- read.xlsx("./aml数据/MV4__11_genecounts.xlsx")
# 检验数据格式
head(dat,n=8) 
# 提取geneid列作为行名
colnames(dat)[1]<-"gene_id"
rownames(dat)<-dat$gene_id
dat<-dat[,-1]
# group信息导入----
colData <- read.xlsx('./aml数据/MV4__11_group.xlsx')
# 检验数据格式
head(colData,n=8) 
# 显示分组信息
table(colData$group) 
colData$group <- factor(colData$group, levels = c("MV4_11_WT","MV4_11_D"))
# R会默认将因子向量的第一个水平作为参考（对照）组，然后将其他水平与这个参考组进行比较。
table(colData$group)
# 让表达矩阵的样本与分析信息表保持一致
colData$id<-gsub("-","_",colData$id) # 根据需要改变文件内容，非必须
dat <- dat[, colData$id] 

# 创建DESeqDataSet对象 ----
dds <- DESeqDataSetFromMatrix(countData = dat, colData = colData, design = ~group)

# 过滤一些低表达的genes
dds <- dds[rowSums(counts(dds)) > 0, ]
dds <- dds[ rowSums(counts(dds[,1:6])) > 15, ]
dds <- dds[ rowSums(counts(dds[,7:9])) > 15, ]

# 估计每个样本的大小因子，用于归一化原始的表达量数据
dds <- estimateSizeFactors(dds) 

# 差异表达分析 ----
dds <- DESeq(dds)

# 提取结果，results函数中contrast参数来指定比较时
# contrast <- c("condition", "treated", "control")
# R会默认将因子向量的第一个水平作为参考（对照）组，
# 然后将其他水平与这个参考组进行比较。
# 所以要将因子水平对调一下 rev(levels(group)))
## result差异表达分析 ----
res <- results(dds, contrast = c("group",rev(levels(colData$group))))
DE <- as.data.frame(res) 
## 删除NA值 ----
DE<-na.omit(DE)

# volcano plot ------------------------------------------------------------
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)
# group 1为实验组
# group 2为对照组
group_1 <- "MV4_11_D"
group_2 <- "MV4_11_WT"
res_data <- DE
data <- res_data[res_data[,"padj"] <= 1,]
# 颜色划分padj <0.05，且log2FC绝对值大于sd(tempOutput$logFC)*3
# 设定cutoff 
# 计算阈值
cutoff <- sd(res_data$log2FoldChange) * 1

data$sig[data$padj >= 0.05 | abs(data$log2FoldChange < cutoff)] <- "Not"

data$sig[data$padj < 0.05 & data$log2FoldChange >= cutoff] <- "Up"

data$sig[data$padj < 0.05 & data$log2FoldChange <= -cutoff] <- "Down"
# table统计上下调基因数量
table(data$sig)
# 得到DEseq2结果
input <- data
# DEseq2结果导出 ----
# 从初始数据中提取差异基因
output<-dat[rownames(dat)%in%rownames(input),]
# 合并DEseq2数据
DEseq2_output<-cbind(output,input)
# 将结果写入相应文件名导出
write.csv(DEseq2_output,file = "./aml数据/MOLM13_DEseq2.csv")

library(ggrepel)
library(ggplot2)
# 提取基因名
y <- rownames(data)
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),"\\|")[[1]][2]))
id<-data.frame(gene)
data<-cbind(id,data)
# 筛选上调中显著性top10的基因
top10_up <- filter(data, sig == 'Up') %>%  # 筛选出上调的基因
  distinct(gene, .keep_all = TRUE) %>%   # 去除重复的基因，保留第一个出现的行
  top_n(10, -log10(padj))               # 选择-P.Value值最大的前10个基因

# 筛选下调中显著性top10的基因
top10_down <- filter(data, sig == 'Down') %>%  # 筛选出上调的基因
  distinct(gene, .keep_all = TRUE) %>%   # 去除重复的基因，保留第一个出现的行
  top_n(10, -log10(padj)) 
# 绘图
volc <- ggplot(data = data, aes(x = log2FoldChange,
                                y = -log10(padj),
                                color = sig)) +
  geom_point(alpha = 0.9) +  theme_classic() +
  theme(panel.grid = element_blank(),strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#4DBBD5","grey80","#E64B35")) +
  geom_vline(xintercept = c(-cutoff,cutoff),lty = 4,lwd = 0.6,alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05),lty = 4,lwd = 0.6,alpha = 0.8) +
  labs(title = paste0("Volcano Plot: ",group_1,"-",group_2))+
  xlim(c(-8, 8))+
  geom_label_repel(data = top10_up, 
                   aes(x = log2FoldChange, y = -log10(padj), label = gene),
                   show.legend = FALSE,size = 2,max.overlaps = Inf)+  # 添加上调基因的标签
  geom_label_repel(data = top10_down, 
                   aes(x = log2FoldChange, y = -log10(padj), label = gene),
                   show.legend = FALSE,size = 2,max.overlaps = Inf)  # 添加下调基因的标签

volc
# 检查log2FoldChange值高的基因 ----
# 按照 log2FC 值从大到小排序
DE_results_sorted <- DE[order(-abs(DE$log2FoldChange)), ]
# 提取 log2FC 值最大的前10个基因
top10_genes <- DE_results_sorted[1:10, ]
top10_genes_id<-rownames(top10_genes)
top10_genes_data<-dat[top10_genes_id,]
write.csv(top10_genes_data,file = "./aml数据/OCI_AML2_top10_genes_data.csv")
# 检查结果发现部分gene表达只在较少基因当中表达，需要进行进一步筛选
