# 样本箱线图+散点图+小提琴图 ----
library(ggplot2)
library(ggpubr)
library(dplyr)

# 定义变量
gene <- "CDK1"                    # 目标基因

# output category
dir <- "./Beat-AML/03_Result/Venetoclax/Gene Expression/"

# expr input
expr <- read.csv("./Beat-AML/01_Data/Venetoclax/Expr_corrected.csv",row.names = 1)

# group input
group <- read.xlsx("./Beat-AML/01_Data/Venetoclax Response Top.xlsx")
group$lab_id <- paste0("X",group$lab_id)
group$lab_id <- gsub("-",".",group$lab_id)

# select targeted gene
targeted_gene <- expr[rownames(expr) == gene, ] # 使用 == 精准匹配

# check targeted gene
print(rownames(targeted_gene)) 

if(T){
# 数据整理
expr_targeted <- targeted_gene
expr_targeted <- t(expr_targeted)
colnames(expr_targeted)[1] <- "expr"             # 所有基因统一修改为expr列名，方便后续绘图
expr_targeted <- as.data.frame(expr_targeted)
expr_targeted$group <- group$group


# 查看样本数（假设 group1 是分组列）
dt <- expr_targeted
num_Sensitivity <- sum(dt$group == "Sensitivity")
num_Resistance <- sum(dt$group == "Resistance")

# 确保分组列是因子，并指定顺序
dt$group <- factor(dt$group, levels = c("Sensitivity", "Resistance"))
rownames(dt)


# 定义颜色
mycol <- c("#6388B4","#EF6F6A")   # 配色

# 绘图
p1 <- ggplot(data = dt,
             aes(x = group, 
                 y = expr,        # 列名为 expr
                 color = group, 
                 fill = group)) +
  geom_violin(alpha = 0.6,        # 先画小提琴图（底层）
              width = 0.9,        # 调节小提琴的宽度
              adjust = 0.9,         # 曲线平滑度,越小越陡峭
              trim = TRUE,       # 是否修剪边缘
              scale = "width",    # 宽度标准化
              linewidth = 1) +  # 轮廓线粗细
  geom_boxplot(alpha = 0.1,       # 再画箱线图（中层）
               width = 0.4,       # 调节箱体宽度
               color = "white",
               outlier.shape = NA,
               coef = 1.5,            # 须的长度（异常值范围）
               linewidth = 1.2) + # 箱线轮廓粗细
  geom_jitter(alpha = 0.7,        # 最后叠加散点（顶层）
              width = 0.2,        # 点的抖动距离 
              shape = 21, 
              size = 2.8,
              aes(fill = group), 
              stroke = 0.8) +
  theme_classic() +  
  scale_color_manual(name = "Group", values = mycol) +
  scale_fill_manual(name = "Group", values = mycol) +
  labs(title = paste("The expression of",gene),
       x = NULL, 
       y = "Normalized CPM") + 
  #scale_y_continuous(breaks = seq(1, 6, by = 0.5)) +  # 网格线间距5单位
  theme(plot.title = element_text(hjust = 0.5, family = "sans" , face = "bold"),  # 标题居中
        aspect.ratio = 1.2,                      # 设置纵横比，调整为更高
        axis.text.x = element_text(size = 12, color = "#333333", face = "bold"),  # x轴文本
        axis.title.y = element_text(size = 12, family = "sans", color = "#333333", face = "bold"), 
        #panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),  # 背景水平网格线
        # 图例标题 
        legend.title = element_text(
          family = "serif",        # 字体类型（如serif衬线体）
          size = 12,               # 字体大小
          face = "bold",           # 加粗
          color = "#333333",       # 深灰色
          margin = margin(b = 5)   # 标题与标签的间距
        ),
        
        # 图例标签 
        legend.text = element_text(
          family = "sans",         # 无衬线字体（如Arial）
          size = 11,               
          color = "#666666",       # 中灰色
          margin = margin(r = 10)  # 标签右侧间距
        ),
        
        # 整体布局 
        legend.position = "none",  # 不显示图例
        legend.spacing.y = unit(0.3, "cm"),  # 图例项垂直间距
        legend.key.size = unit(0.8, "cm"),   # 图例符号大小
        legend.background = element_rect(fill = "white", color = NA))  # 背景优化 
# 添加统计学检验（可选）
p1 <- p1 + stat_compare_means(comparisons = list(c("Sensitivity", "Resistance")),
                              method = "t.test",  # 非参数检验 wilcox.test, 参数检验 t.test 
                              label = "p.signif")  # 显示显著性标记（***, ​**, *）
print(p1)
}

ggsave(filename = paste0(dir, gene, "_expr_boxplot.pdf"),
       plot = p1, device = cairo_pdf, 
       width = 6, height = 5)        
