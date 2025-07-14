# limma


## Voom with sample quality weights ----
# pkgs
library(GEOquery)
library(limma)
library(edgeR)

# 导入GSE64099 的 series matrix 文件
example_data <- read.table("Limma/01_data/GSE55713_series_matrix.txt.gz",sep="\t",quote = "",
                fill=T,comment.char = "!",header=T)
rownames(example_data)<-example_data[,1] #把第一列的值变为行名
example_data<-example_data[,-1] 
example_group <- colnames(example_data)
example_group <- as.data.frame(example_group)
# 正确的分组向量（不需要是 data.frame）
example_group <- factor(c("Ctrl", "Ctrl", "Ctrl", "TGIF1", "TGIF1", "TGIF1"))
x <- DGEList(counts = example_data, group = example_group)
col <- c("Ctrl" = "blue", "TGIF1" = "red")[example_group]

# PCA （非经典主成分分析算法） 
plotMDS(x,
        gene.selection = "common",  # ✅ 强制使用 PCA 模式
        top = 1000,                 # 用标准差最大的前 1000 个基因
        col = col,                  # 每个样本颜色
        labels = c(1:6),              # 不显示样本名（可设为 colnames(x)）
        main = "PCA plot")
legend("topleft", legend=c("Ctrl", "TGIF1"), col=1:2, pch=15)

# PCoA （主坐标分析）
plotMDS(x, 
        col = col, 
        top = 1000, 
        labels = c(1:6), 
        main = "MDS plot")
legend("topleft", legend=c("Ctrl", "TGIF1"), col=1:2, pch=15)

# 创建设计矩阵：默认因子第一个水平
design <- model.matrix(~-1 + example_group)
# Analysis with voom only
design[1:6,]
colnames(design) <- c("Ctrl", "TGIF1")  # 重命名，方便后续分析
# 创建对比矩阵
contrast <- makeContrasts(TGIF1vsCtrl = TGIF1 - Ctrl, levels = design)
v <- voom(x, design = design)
fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
options(digits = 3)  # 只显示三位有效数字 不改变实际值
topTable(fit2, coef = "TGIF1vsCtrl", sort.by="P")
alltable <- topTable(fit2, coef = "TGIF1vsCtrl", number = Inf, sort.by="P")
sum(alltable$adj.P.Val < 0.05)

# Analysis with combined voom and sample quality weights
vwts <- voomWithQualityWeights(x, design=design, normalize.method ="none", plot=TRUE) 
vfit <- lmFit(vwts) 
vfit2 <- contrasts.fit(vfit, contrast)
vfit2 <- eBayes(vfit2)  
topTable(vfit2, coef="TGIF1vsCtrl", sort.by="P")
alltable2 <- topTable(vfit2, coef = "TGIF1vsCtrl", number = Inf, sort.by="P")
sum(alltable2$adj.P.Val < 0.05)
