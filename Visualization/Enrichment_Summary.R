#' 此脚本为富集分析可视化总结，仅包含富集分析可视化部分, 需要提供各种可视化方法
#' 使用的数据结构，如用dotplot（）就必须提供ClusterProfiler包富集的结果
#' 

# 1. Clusterprofiler
library(enrichplot)
library(stringr)
library(ggplot2)
library(dplyr)
# data input
GO <- go
# 选择ontology （可选）
GO_filtered <- GO@result %>%
  filter(ONTOLOGY %in% c("BP")) 

# 重新构建 enrichResult 对象，只保留 BP 和 CC 项
GO_subset <- GO
GO_subset@result <- GO_filtered

# plot GO
p1 <- dotplot(GO_subset, x = "GeneRatio", color = "p.adjust", 
              showCategory = 20, split = NULL, label_format = 30) + 
  # facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free') +
  theme(axis.text.y = element_text(angle = 0, hjust = 1)) + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40))  # 控制每行最多显示40个字符
print(p1)
ggsave("./03_Result/GO&KEGG/OCI_AML2_single_fill/Low_vs_Con/GO_down_BP_dotplot.pdf", plot = p1, width = 6.5, height = 8)

#plot KEGG
