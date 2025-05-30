# 加载所需的库
library(ggplot2)

# 普通箱线图 ----
# 示例数据集：使用内置的 iris 数据集
data(iris)
iris <- as.data.frame(iris)
# 绘制箱线图
ggplot(iris, aes(x = Species, y = Sepal.Length)) + 
  geom_boxplot(
    aes(fill = Species),          # 箱子的填充颜色根据分类变量 (Species)
    outlier.colour = "red",        # 异常值的颜色
    outlier.size = 3,              # 异常值的大小
    outlier.shape = 16,            # 异常值的形状（16是点的形状）
    varwidth = TRUE,              # 是否根据样本量调整箱体宽度
    width = 0.6,                  # 箱子的宽度
    alpha = 0.6                   # 箱子的透明度
  ) +
  # 设置主题
  theme_minimal() +
  # 设置标题、坐标轴标签等
  labs(
    title = "Sepal Length by Species",   # 图表标题
    x = "Species",                      # x轴标签
    y = "Sepal Length (cm)"             # y轴标签
  ) +
  # 自定义坐标轴
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),  # x轴标题字体设置
    axis.title.y = element_text(size = 14, face = "bold"),  # y轴标题字体设置
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1), # x轴标签字体及旋转角度
    axis.text.y = element_text(size = 12)  # y轴标签字体设置
  )

# 样本配对箱线图 ----
# load data
expr <- read.csv("./Beat-AML/01_Data/Gene Counts CPM Venetoclax.csv",row.names = 1)
targeted_gene <- expr[grep("PHGDH",rownames(expr)),]
expr_targeted <- targeted_gene[1,grep("MV4_11",colnames(targeted_gene))]
expr_targeted <- t(log2(expr_targeted))
colnames(expr_targeted)[1] <- "expr"
expr_targeted <- as.data.frame(expr_targeted)

group <- read.xlsx("./Beat-AML/01_Data/Venetoclax Response Top.xlsx")
group_targeted <- group[grep("MV4_11",group$id),c(1,2)]
expr_targeted$group <- group_targeted$group

# 定义变量
gene <- "TP53"          # 目标基因
mycol <- c("#6388B4","#FFAE34","#EF6F6A")  # 配色

# 计算样本数（假设 group1 是分组列）
dt <- expr_targeted
num_normal <- sum(dt$group == "Ctrl")
num_low_resistance <- sum(dt$group == "Low")
num_high_resistance <- sum(dt$group == "High")

# 确保分组列是因子，并指定顺序
dt$group <- factor(dt$group, levels = c("Ctrl", "Low", "High"))
rownames(dt)
dt$pair <- c(1,1,1,2,3,2,3,2,3)
table(dt$pair)

dir <- "./03_Result/DEP/MV4_11/"

# 绘图
p1 <- ggplot(data = dt,
             aes(x = group, 
                 y = expr,  # 假设表达量列名为 expression
                 color = group, 
                 fill = group)) +
  geom_line(aes(group = pair), color = "gray50", linewidth = 0.5, alpha = 0.6) +  # 加入配对线
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +  # 隐藏箱线图异常值
  geom_jitter(width = 0.2, size = 2.8, shape = 21, 
              stroke = 0.8,aes(fill = group), 
              alpha = 0.7) +
  theme_classic() +  scale_color_manual(name = "Group", values = mycol) +
  scale_fill_manual(name = "Group" ,values = mycol) +
  labs(
    title = paste("Expression of", gene, "across three states"),
    x = "MV4_11 Group",
    y = "Log2(expr + 1)"
  ) + 
  theme(plot.title = element_text(hjust = 0.5),  # 标题居中
        aspect.ratio = 1.2,                      # 设置纵横比，调整为更高
        axis.title = element_text(size = 12, face = "plain"),  # 轴标题加粗，字号14
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
        legend.position = "right",
        legend.spacing.y = unit(0.3, "cm"),  # 图例项垂直间距
        legend.key.size = unit(0.8, "cm"),   # 图例符号大小
        legend.background = element_rect(fill = "white", color = NA))  # 背景优化  
# 添加统计学检验（可选）
p1 <- p1 + stat_compare_means(comparisons = list(c("Ctrl", "Low"),
                                                 c("Low", "High"),
                                                 c("Ctrl", "High")),
                              method = "t.test",  # 非参数检验 wilcox.test, 参数检验 t.test 
                              label = "p.signif")  # 显示显著性标记（***, ​**, *）

ggsave(filename = paste0(dir,"/TP53_expr_boxplot.pdf"),
       plot = p1, device = cairo_pdf, 
       width = 6, height = 5)
