# 提取差异基因结果
# 设置工作路径 ----
setwd("../AML_project/aml_analysis")

# 加载包
library(dplyr)
library(openxlsx)

# 导入DE数据
data <- read.csv("./03_result/DE_limma/MOLM13_WT_vs_MOLM13_D/DE.csv")

# 加change列
logFC = 1
P.Value = 0.05
k1 <- (data$logFC > logFC)&(data$P.Value < P.Value) 
k2 <- (data$logFC < -logFC)&(data$P.Value < P.Value)
data <- mutate(data,change = ifelse(k1,"up",ifelse(k2,"down","stable")))
table(data$change)

# 提取差异基因，cutoff <- 1，padj <0.05
data_1 <- data[data$change!="stable",]
y <- data_1$Row.names
gene <- unlist(lapply(y, function(y)strsplit(as.character(y),"\\|")[[1]][2]))
gene <- as.data.frame(gene)

# 导出差异基因列表
write.xlsx(gene,file = "03_result/DE_limma/MOLM13_WT_vs_MOLM13_D/Genes.csv")
