# VolcanoPlot_Summary
# 注释火山图 ----
library(openxlsx)
library(ggplot2)
library(ggrepel)  # 避免标签重叠

DE_res <- read.xlsx("./03_result/02_DE/OED_vs_SCD/DE_Protein_Names_cut-pvalue.xlsx")
DE_res$Sig <- factor(DE_res$Sig, levels = c("up", "stable"))
max_abs_logfc <- max(abs(DE_res$logFC), na.rm = TRUE)
x_limit_right <- max_abs_logfc * 1.05
# 如果要求对称则不需要设置左界限，直接设置为 -x_limit_right 即可
x_limit_left <- min(DE_res$logFC) * 1.2

# 输出目录
output_dir <- "./03_result/02_DE/OED_vs_SCD/"

# 标记分组
group_1 <- "OED"    # 对照组
group_2 <- "SCD"        # 实验组

# 标记阈值
logfc_threshold <- 1
pvalue_threshold <- 0.00002	


# 提取显著差异基因在火山图上标记
# 比如 top 10 up 
top_up <- DE_res[DE_res$Sig == "up", ]
top_up <- top_up[order(top_up$logFC,decreasing = TRUE), ][1:10, ]
# top_down <- DE_res[DE_res$Sig == "down", ][order(DE_res$P.Value), ][1:10, ]
label_data <- top_up

## plot ----
p1 <- ggplot(data = DE_res, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha = 0.5, size = 3, 
             aes(color = Sig)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(
    name = "Change",
    values = c("up" = "#B30000", "stable" = "grey", "down" = "#003366"),
    labels = c("up" = paste0("Up ：", sum(DE_res$Sig == "up")),
               "stable" = paste0("Stable ：", sum(DE_res$Sig == "stable")),
               "down" = paste0("Down ：", sum(DE_res$Sig == "down"))))+
  geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), lty = 4, 
             col = "black", lwd = 0.8, alpha = 0.4) +
  geom_hline(yintercept = -log10(pvalue_threshold), lty = 4, 
             col = "black", lwd = 0.8, alpha = 0.4) +
  geom_text_repel(data = label_data,   # 标签
                  aes(label = gene), 
                  size = 3.5,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  arrow = arrow(length = unit(0.008,"npc")), # 箭头指向
                  segment.size = 0.5, min.segment.length = 0.5,            # 确保箭头显示
                  force =15, # 标签排斥力
                  max.overlaps = 20) +
  labs(title = paste0(group_2,"-",group_1)) +
  xlim(x_limit_left, x_limit_right)+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(), # x轴网格线
        panel.grid.minor = element_blank(), # y轴网格线
        aspect.ratio = 1)
print(p1)
ggsave(filename = paste0(output_dir,"anno_volc_cut-pvalue.pdf"),
       plot = p1, device = "pdf", 
       width = 6, height = 5)
