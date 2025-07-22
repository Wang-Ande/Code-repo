# Title: KEGG分析本地部署（手动下载 KEGG 映射关系）
# Author: [wangzihao]
# Date: [2025-07-10]
# Description:
#   此脚本展示如何手动下载 KEGG Pathway 数据，并格式化为 
#   `clusterProfiler::enricher()` 可识别的 `TERM2GENE` 与 `TERM2NAME` 格式。
#   适用于不借助 `createKEGGdb` 包、也无需在线访问 KEGG API 的情况下，
#   依然能稳定完成 KEGG 富集分析。
#
# 步骤说明:
#   1. 从网站上下载KEGG Pathway ID 列表（pathway + description）
#   链接：https://rest.kegg.jp/list/pathway/hsa ，右键另存为pathway_list.txt
#   2. 从网站上下载基因与 KEGG Pathway 的映射关系
#   链接：https://rest.kegg.jp/link/pathway/hsa ，右键另存为gene2pathway.txt
#
# 优点:
#   - 脱离 createKEGGdb 和 KEGG API，可灵活控制版本与更新
#   - 可在无网络或 KEGG 访问不稳定时使用


# pkgs ----
library(clusterProfiler)

# 读取你刚才下载的文件
gene2pathway <- read.table("01_Data/gene2pathway.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(gene2pathway) <- c("gene", "pathway")
gene2pathway$gene <- gsub("hsa:", "", gene2pathway$gene)
gene2pathway$pathway <- gsub("path:", "", gene2pathway$pathway)

# 自己定义 TERM2GENE 数据框
TERM2GENE <- gene2pathway[, c("pathway", "gene")]

# 可选：读入 pathway 名称，作为 TERM2NAME
pathway_list <- read.table("01_Data/pathway_list.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(pathway_list) <- c("pathway", "name")
pathway_list$pathway <- gsub("path:", "", pathway_list$pathway)
TERM2NAME <- pathway_list

# 你的基因列表（ENTREZ ID）
# DE_res input ----
DP_result <- read.csv('./03_Result/DEP/OCI_AML2/Low_vs_Ctrl/result_DE.csv')

# set P.Value ----
GeneSymbol <- subset(DP_result, P.Value < 0.05)

# set logFC ----
cutoff <- 0.263                 # 对应fc约为1.2

# 转换基因名 
y <- GeneSymbol$Genes
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
GeneSymbol$gene <- gene

# 设置数据库 
GO_database <- 'org.Hs.eg.db'  # GO是org.Hs.eg.db数据库
KEGG_database <- 'hsa'         # KEGG是hsa数据库

# gene list ----
down_genes <- subset(GeneSymbol, logFC < -cutoff)

# gene ID转换 
gene <- clusterProfiler::bitr(down_genes$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)

# 运行富集分析
kegg_result <- enricher(gene = gene$ENTREZID,
                        TERM2GENE = TERM2GENE,
                        TERM2NAME = TERM2NAME,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        qvalueCutoff = 1)
