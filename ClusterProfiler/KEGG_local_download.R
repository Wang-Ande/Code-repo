# Title: KEGG分析本地部署 (Human, hsa)
# Author: [wangzihao]
# Date: [2025-07-10]
# Description:
#   此脚本展示如何下载适用于Clusterprofiler包的KEGG数据格式(via `createKEGGdb`)。
#   方便使用者本地进行KEGG富集分析，而不依赖KEGG API 路径。
#   优点：适用于网络卡顿或者需要保证结果统一的分析。
#   缺点：需要保证本地KEGG数据的更新。

# pkgs ----       
# pak("YuLab-SMU/createKEGGdb") #创建KEGG数据库的包的包
# library(createKEGGdb)
# 创建名称为KEGG.db_1.0.tar,gz的包。
# createKEGGdb::create_kegg_db("hsa")   # 需要保持更新！！！！！！！！！！！
# 安装这个包(默认的包的路径在当前工作目录，根据实际情况修改路径)
# install.packages("KEGG.db_1.0.tar.gz",repos=NULL,type="source")
library(KEGG.db)
library(clusterProfiler)
library(pathview)
library(org.Hs.eg.db)
library(readxl)
library(readr)
library(dplyr)
library(openxlsx)
library(ggplot2)

# DE_res input ----
DP_result <- read.csv("your DE_res") 

# set P.Value ----
GeneSymbol <- subset(DP_result, P.Value < 0.05)

# set logFC ----
cutoff <- 0.263                 # 对应fc约为1.2

# 转换基因名 
y <- GeneSymbol$Genes
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
GeneSymbol$gene <- gene

# gene list ----
down_genes <- subset(GeneSymbol, logFC < -cutoff)

# gene ID转换 
gene <- clusterProfiler::bitr(down_genes$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)
library(KEGG.db)
kk <- enrichKEGG(gene= gene$ENTREZI,
                 keyType = "kegg",
                 organism= 'hsa',
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 1,
                 readable = T,
                 use_internal_data =T) # 当使用本地的KEGG.db数据库时，use_internal_data设置为T

