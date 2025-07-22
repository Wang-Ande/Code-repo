library(DESeq2)
library(ggplot2)
library(factoextra)
# PCA绘图总结 ----
## 获取PCA数据 ----
# 导入数据矩阵
data <- read.csv("./01_data/aml数据/data_merge_adjusted.csv")
rownames(data) <- data$X
data <- data[,-1]

# 清洗数据，去除标准差为零的行
data <- data[apply(data,1,sd)!=0,]
data <- t(data)
# 导入分组信息
colData <- read.xlsx("./01_data/aml数据/aml_group_merge.xlsx")

# 执行 PCA
pca_result <- prcomp(data,scale. = T) 

# PCA运算之前需要将数据进行转置t（），使行变为样本
# 且数据要经过标准化处理，如运用scale（）

# 分组信息设置
group_list <- colData$group
# 批次信息设置

## factoextra绘图方法 ----
id_list<-sapply(as.character(colData$id), 
                function(x)substr(x,nchar(x)-3,nchar(x)))
pdf(
  "./PCA.pdf")

a <- fviz_pca_ind(pca_result,
                  axes = c(1, 2),          # 显示 PC1 和 PC2
                  geom.ind = c("point"),   # 显示点和文本
                  col.ind = group_list,  # 根据分组着色
                  #shape.ind = batch_list,  # 设定批次信息
                  palette = c("#e74645", "#56B4E9","#F2A200"), # 自定义颜色
                  #shape = c(16,17) , # 不同的批次会被绘制为实心圆形、实心三角形
                  addEllipses = TRUE,  # 添加置信区间椭圆
                  ellipse.type = "t", # 置信区间椭圆,参数confidence(大样本)/t
                  legend.title = "group",
                  title = "PCA of  Dataset") +
  geom_text(aes(label = id_list,color = group_list),nudge_x = -1.5)+
  theme_classic()
dev.off()

## ggplot绘图方法 ----
# 创建数据框用于 ggplot2
pca_df <- data.frame(PC1 = pca_result$x[,1], 
                     PC2 = pca_result$x[,2], 
                     group = group_list,
                     batch = batch_list)

# 绘制 PCA 图
ggplot(pca_df, aes(x = PC1, y = PC2, color = group, shape = batch)) +
  geom_point(size = 3) +
  labs(title = "PCA",
       x = "PC1",
       y = "PC2") +
  theme_minimal()
