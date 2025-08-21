#' 此脚本为富集分析可视化总结，仅包含富集分析可视化部分, 需要提供各种可视化方法
#' 使用的数据结构，如用dotplot（）就必须提供ClusterProfiler包富集的结果
#' 

# 1. Clusterprofiler ----
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

# 2. ggplot2
## 2.1 dotplot
# ggplot2画气泡图，scale_color_gradient设置蓝红配色
library(openxlsx)
library(ggplot2)
library(stringr)
library(scales)

# Enrich_res input 
Enrich_res <- read.xlsx("./03_result/03_Enrichment/OE_A175fs_vs_Control/KEGG_up.xlsx")

# Kegg pathway blacklist screening
kegg_blacklist <- read.csv("D:/R/RStudio/Code_Box/ClusterProfiler/kegg_blacklist.csv")
Enrich_res$ID <- gsub("hsa0", "", Enrich_res$ID)
Enrich_res <- Enrich_res[!Enrich_res$ID%in%kegg_blacklist$id,]

# pvalue threshold
Enrich_res <- Enrich_res[Enrich_res$pvalue<0.05,]

# GeneRatio compute
Enrich_res$GeneRatio <- round(
  as.numeric(sub("/.*", "", Enrich_res$GeneRatio)) / 
    as.numeric(sub(".*/", "", Enrich_res$GeneRatio)),
  3
)

# 设置排列顺序，按X轴大小
Enrich_res$Description <- factor(
  Enrich_res$Description, 
  levels = Enrich_res$Description[order(Enrich_res$GeneRatio, decreasing = FALSE)]
)

# 绘图
p2 <- ggplot(Enrich_res, aes(x = GeneRatio, y = Description)) +
  geom_point(
    aes(size = Count, fill = pvalue),  # 用 fill 映射颜色
    shape = 21,                        # 有边框的点
    color = "black",                   # 边框颜色
    stroke = 0.6
  ) +
  theme_bw() +
  labs(y = "", x = "GeneRatio") +
  scale_size(name = "Count") +         # 自定义图例标题
  scale_fill_gradient(name = "pvalue", low = "#e06663", high = "#327eba") +  # 用 fill 控制颜色
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +  # x轴保留两位
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  scale_size_continuous(range = c(2, 8))+  # 调整气泡的大小范围
  guides(
    fill = guide_colorbar(order = 1),  # pvalue 图例在上面
    size = guide_legend(
      order = 2,
      override.aes = list(
        shape = 21,           # 必须指定，不然没边框
        color = "black",      # 边框颜色
        fill = "transparent"  # 填充透明
      )
    )
  ) +
  theme(
    text = element_text(family = "sans", colour = "black"),
    axis.title.x = element_text(size = 12), 
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(angle = 0, hjust = 1 ,size = 12, colour = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 12))
x_max <- max(Enrich_res$GeneRatio, na.rm = TRUE)
p2 <- p2 + 
  scale_x_continuous(
    labels = scales::number_format(accuracy = 0.01),
    limits = c(min(Enrich_res$GeneRatio, na.rm = TRUE), x_max * 1.05)  # 右边界比最大值大 5%
  )
print(p2)
# 用 ggsave 保存
ggsave(
  filename = "your_output_pathway.pdf",
  plot = p2,
  width = 5,
  height = 6,
  units = "in"   # 单位可选 "in", "cm", "mm"
)
# scale_size_continuous(range = c(3, 12))+  # 调整气泡的大小范围
