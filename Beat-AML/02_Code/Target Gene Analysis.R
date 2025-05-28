# Pipline ----
#' @description
#' 表达矩阵 → SVA矫正（surrogate variables）→ limma-trend差异分析 → FDR筛选 → 
#' removeBatchEffect校正表达矩阵 → Mclust聚类 → clusplot可视化

# Library packages ----
library(openxlsx)
library(readr)
library(dplyr)
library(limma)
library(edgeR)
library(sva)
library(mclust)
library(cluster)

# Data input ----
Drug_Responses <- read.xlsx("./Beat-AML/01_Data/Drug Responses.xlsx")
Ven_Response <- Drug_Responses[grep("Venetoclax",Drug_Responses$inhibitor),]
Gene_Cpm <- read.xlsx("./Beat-AML/01_Data/Gene Counts CPM.xlsx")
rownames(Gene_Cpm) <- Gene_Cpm$Symbol
Gene_Cpm <- Gene_Cpm[,-c(1,2)]

# Group input ----
# 根据文献中的介绍，按照分别取AUC的前后20%作为resistance和sensitivity
# 病人样本仅取在sensitive or resistant中的样本
df <- Ven_Response
df <- df %>%
  mutate(group = case_when(
    auc <= quantile(auc, 0.2, na.rm = TRUE) ~ "Sensitivity",    
    auc >= quantile(auc, 0.8, na.rm = TRUE) ~ "Resistance",
    TRUE ~ "Moderate"
  ))                                    # quantile按升序排列c(1,2,3,4)

# 查看 n
table(df$group)
Ven_Response <- df

# 筛选只在sensitive和resistant中的转录组测序样本
filtered_id <- Ven_Response[Ven_Response$group!= "Moderate", "lab_id"]
Gene_Cpm <- Gene_Cpm[,colnames(Gene_Cpm)%in%filtered_id]

Ven_Response <- Ven_Response[Ven_Response$lab_id %in% colnames(Gene_Cpm),]
table(Ven_Response$group)

# res output
write.csv(Gene_Cpm, file = "./Beat-AML/01_Data/Gene Counts CPM Venetoclax.csv")
write.xlsx(Ven_Response, file = "./Beat-AML/01_Data/Venetoclax Response Top.xlsx")

# PCA analysis ----
library(factoextra)
# PCA计算
data_pca <- t(expr_corrected)
pca_result <- prcomp(data_pca,scale. = T) 

# 设置分组
group_list <- Ven_Response$group

# 绘制 PCA 并保存为 PDF
pdf("./Beat-AML/03_Result/PCA_corrected.pdf", height = 6.5, width = 6.5)
pca1 <- fviz_pca_ind(pca_result,
                  axes = c(1, 2),             # 显示 PC1 和 PC2
                  geom.ind = c("point"),      # 显示点和文本
                  col.ind = group_list,       # 根据分组着色
                  #shape.ind = batch_list,    # 设定批次信息
                  palette = c("#e74645", "#56B4E9","#F2A200"), # 自定义颜色
                  #shape = c(16,17) ,         # 不同的批次会被绘制为实心圆形、实心三角形
                  addEllipses = TRUE,         # 添加置信区间椭圆
                  ellipse.type = "t",         # 置信区间椭圆,参数confidence(大样本)/t
                  legend.title = "Group", 
                  title = "PCA of  Dataset",
                  label = "all",              # 显示所有样本的 ID
                  repel = TRUE)+              # 防止标签重叠（推荐）
      theme_classic()
print(pca1)
dev.off()

# Expr mat sva correct ----
# 表达矩阵: gene x sample
# 假设表达矩阵已经是log2 CPM或RPKM (normalized counts)
expr_mat <- log2(Gene_Cpm + 1)

# 分组信息
# group should be a factor with levels "sensitive" and "resistant"
group <- factor(Ven_Response$group, levels = c("Sensitivity","Resistance"))           # 主要关注耐药状态的影响
design <- model.matrix(~group)

# 计算 SVs 以矫正批次效应
# 构建基础模型矩阵（仅包含已知生物变量）
mod <- model.matrix(~group)  

# 构建空模型矩阵（仅含截距项）
mod0 <- model.matrix(~1, data=data.frame(group=group))

# 方差过滤,选择高变基因
library(matrixStats)
top_var_genes <- order(rowVars(as.matrix(expr_mat)), decreasing = TRUE)[1:2000]
expr_filtered <- expr_mat[top_var_genes, ]

# 使用svaseq推断隐变量（SVs）
n.sv <- min(num.sv(expr_filtered, mod, method="leek"),15)   # 限制最大数量sv为15 ，“leek”为推荐方法，
svobj <- sva(as.matrix(expr_mat), mod, mod0, n.sv = 15)     # 选择n.sv = 15
fsvaobj <- fsva(dbdat = as.matrix(expr_mat),mod = mod,sv = svobj, 
             newdat = matrix(nrow = nrow(expr_mat), ncol = 0)) # 函数源代码有误，需要提供一个数据集作为newdata
expr_corrected <- fsvaobj$db                                # 与removeBatchEffect矫正结果一样，二选一即可

# 查看不同SVs与PC1 explained vraiance 的关系
source("./Beat-AML/02_Code/test_sv_range.R")
result <- test_sv_range(expr_mat, mod, group, sv_range = 0:30)

# Check over correct/fit ----
source("./Beat-AML/02_Code/check_sva_overfitting.R")
check_sva_overfitting(expr_mat = expr_mat, 
                      group = group, 
                      svobj = svobj, 
                      mod = mod)

# 获取矫正后的表达矩阵（如果用fsvaobj得到矫正后的矩阵，即可省略）
expr_corrected <- removeBatchEffect(expr_mat, covariates = svobj$sv, design = mod)
expr_corrected <- as.data.frame(expr_corrected)
write.csv(expr_corrected, file = "./Beat-AML/01_Data/Venetoclax/Expr_corrected.csv")

# DE analysis ----
# 将SVs作为协变量加入设计矩阵
design_sva <- cbind(design, svobj$sv)

# 使用 limma-trend 模式建模
fit <- lmFit(expr_mat, design_sva) 
fit <- eBayes(fit, trend=TRUE)                 

# 提取差异表达结果
res <- topTable(fit, coef="groupResistance", number=Inf, adjust="BH")
sig_genes <- rownames(res[res$adj.P.Val < 0.05, ])

# 合并表达矩阵
expr_merge <- merge(expr_corrected, res, by = "row.names")
write.csv(expr_merge, file = "./Beat-AML/03_Result/Venetoclax/Expr_DE_merge.csv")

# Cluster analysis ----
# 仅用显著差异基因
expr_sig <- expr_corrected[sig_genes, ]

# PCA
pca <- prcomp(t(expr_sig), scale.=TRUE)
pc_data <- pca$x[, 1:2]  # 前两个主成分

# Mclust 聚类
mc <- Mclust(t(expr_sig))  # 注意：Mclust默认用于样本聚类
clusters <- mc$classification

# 可视化聚类
clusplot(pc_data, clusters, color=TRUE, shade=TRUE, labels=2, lines=0, main="Clustering via Mclust")



corrected_expr_1 <- fsva(
  dbdat = as.matrix(expr_mat),  # 需要校正的表达矩阵
  mod = mod,              # 原模型矩阵（与SVA阶段一致）
  sv = svobj,          # SVA结果
  newdat = as.matrix(expr_mat)           # 如果没有新数据，设为NULL
)
corrected_matrix_1 <- corrected_expr_1$db
