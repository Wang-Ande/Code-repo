# WGCNA
rm(ls())     # 清空工作环境
# 加载包
if(!require(biomaRt ,quietly = TRUE)){BiocManager::install("biomaRt")}
library(biomaRt)
library(readr)
library(openxlsx)
library(readxl)
library(dplyr)
library(scales)
library(tidyr)
library(WGCNA)
library(forcats)
library(clusterProfiler)
library(pathview)
library(org.Hs.eg.db)
library(ggplot2)

# 0. set output category ----
folder_path <- "your pathway"
if (!dir.exists(folder_path)) {
  dir.create(folder_path)
  print(paste("Folder", folder_path, "created."))
} else {
  print(paste("Folder", folder_path, "already exists."))
}

dir_WGCNA <- "your pathway"

# 1.data input ----
## 1.1 Expr_input ----
Expr_raw <- read.csv("your pathway")
rownames(Expr_raw) <- Expr_raw$X
Expr_raw <- Expr_raw[,-1]
datExpr <- Expr_raw

## 1.2 标准化 log2+1 ----
datExpr <- log2(datExpr)

## 1.3 filter ---- 
# WGCNA可以不过滤，WGCNA包中有自己的过滤方案
# var计算每个基因方差，筛选基因变化较大的基因，此次选取前75%的基因
if(T){
  vars_res <- apply(datExpr, 1, var)
  
  # 计算百分位数截止值
  per_res <- quantile(vars_res, probs = seq(0, 1, 0.25)) # quantile生成分位数
  per_res
  
  upperGene <- datExpr[ which(vars_res > per_res[4]),]  # 选取方差位于前75%的基因
  dim(upperGene)
  datExpr <- data.matrix(upperGene)
}

# 2 t转置 ----
# 为了执行 WGCNA，通常需要转置数据，使其符合 WGCNA 的标准格式
datExpr <- t(datExpr)                       # 让每一列代表一个基因

# 3. 检查数据好坏 ----
gsg <- goodSamplesGenes(datExpr, verbose = 3)

# 查看标记结果
gsg$allOK
#[1] TRUE为好

# 如果为false则运行下段！
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", 
                     paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", 
                     paste(names(datExpr)[!gsg$goodSamples], collapse = ", ")))
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}

# 对样本聚类查看有无异常样本
sampleTree = hclust(dist(datExpr), method = "average")
png("your pathway", 
    width = 1500, 
    height = 1000,
    res = 300)
par(cex = 0.6)
par(mar = c(0,4,2,0)) 
plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub="", 
     xlab="", 
     cex.lab = 1.5, 
     cex.axis = 1.5, 
     cex.main = 2)
#abline(h = 150, col = 'red')
dev.off()

# 剔除离群样本，无则跳过
clust = cutreeStatic(sampleTree, cutHeight = 150, minSize = 10) 
table(clust)               #显示剔除多少样本， 0代表剔除的个数，1代表保留的个数
datExpr = datExpr[clust == 1, ]

# 4. 选择软阈值 ----
# 定义1:30的power值，用于构建不同的基因共表达网络
powers <- c(c(1:10), 
            seq(from = 12, 
                to = 30,
                by = 2))

sft <-  pickSoftThreshold(datExpr,
                          powerVector = powers,
                          verbose = 5,
                          networkType = 'unsigned')
sft$powerEstimate
# 绘制拟合度图来选择合适的软阈值
pdf(file = "your pathway.pdf",width = 10,height = 6.5)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,
     col="steelblue")
abline(h=0.85,col="red")
plot(sft$fitIndices[,1],  sft$fitIndices[,5], 
     type="n", 
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", 
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,  col="steelblue")
dev.off()

# 使用合适的软阈值（没有合适就选6）
power <- 6

# 5. 构建共表达矩阵 ----
net = blockwiseModules(datExpr,
                       power = power,           # 软阀值选为6
                       TOMType = "unsigned",    # 构建无尺度网络
                       minModuleSize = 30,      # 最小模块基因数为30
                       reassignThreshold = 0,   # 重分配阈值          
                       mergeCutHeight = 0.25,   # 模块合并阀值
                       numericLabels = F,       # 返回字符标签（如模块颜色名称）
                       pamRespectsDendro = FALSE,         
                       saveTOMs = FALSE,         # 不保存TOM矩阵
                       verbose = 3,
                       maxBlockSize = 20000)    # 可处理的最大模块基因数
# 显示所有模块个数和各个模块中基因数量
table(net$colors)

## 5.1 聚类分析 ----
# 使用层次聚类
geneTree = net$dendrograms[[1]] 
geneTree

png("your pathway.png",
    width = 2000,
    height = 1500,
    res = 300,
    type = "cairo")
moduleColors = net$colors
plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],
                    "Module colors",             # 标题
                    dendroLabels = FALSE,        # 不显示基因标签
                    hang = 0.03,                 # 树枝悬挂高度（相对于根部）
                    addGuide = TRUE,             # 颜色引导条
                    guideHang = 0.05)            # 颜色引导条与树状图悬挂距离
dev.off()


# 6. 模块与表型的关系 ----
### 6.1 性状数据输入 ----
# 假设 `traitData` 是你样本的表型数据，可以是疾病分组、临床特征等
trait <- read.csv("./01_Data/trait.csv")
rownames(trait) <- trait$X
trait <- trait[,-1]
colnames(trait)


if(T){
  # 模块特征
  MEs0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes  # 计算模块特征向量
  MEs = orderMEs(MEs0)                                       # 对模块特征向量排序
  
  # 计算模块特征向量与表型的相关系数矩阵
  moduleTraitCor <- cor(MEs,trait,use = "p")  
  
  # 计算相关系数矩阵的p值
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)  
  
  # 构建绘图时用的文本矩阵
  textMatrix = paste(signif(moduleTraitCor,2),"\n(",
                     signif(moduleTraitPvalue,1),")",sep = "")  
  
  # 修改文本矩阵的维度，与相关系数矩阵相同
  dim(textMatrix)=dim(moduleTraitCor)  
}

### 6.2 相关性热图 ----
if(T){
  pdf(file = "your pathway.pdf",
      width = 12,height = 8)
  
  # mar（）分别代表图形边距的底部、左侧、顶部和右侧的边距
  par(mar = c(7, 7, 2, 2)) 
  labeledHeatmap(Matrix = moduleTraitCor,        # 绘制带标签的热图
                 xLabels = colnames(trait),     # x轴标签
                 yLabels = names(MEs),           # y轴标签
                 ySymbols = names(MEs),          # y轴符号
                 colorLabels = FALSE,            # 不显示颜色标签
                 colors = blueWhiteRed(50),      # 颜色范围
                 textMatrix = textMatrix,        # 显示文本矩阵
                 setStdMargins = FALSE,          # 不设置标准边距
                 cex.text = 0.5,                 # 文本大小
                 cex.lab.x = 0.7,                # X轴文本大小
                 zlim = c(-1,1),                 # 颜色映射范围
                 main = paste("Module-trait relationships"))  # 绘图标题
  dev.off()
}

# 提取gene模块信息备用
df <- data.frame(gene = colnames(datExpr), module = moduleColors)
write.csv(df,file = "your pathway.csv")

### 6.3 基因与性状和重要模块的关系：基因重要性和模块成员 ----
Time = as.data.frame(trait$Time)   # 假设关注性状是Time                           
names(Time) = "Time"
modNames = substring(names(MEs), 3)                              # 提取模块名称
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), 
                                          nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")      # 修改列名
names(MMPvalue) = paste("p.MM", modNames, sep="")                # 修改列名

# 基因显著性
geneTraitSignificance = as.data.frame(cor(datExpr, Time, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 
                                          nSamples))
names(geneTraitSignificance) = paste("GS.", names(Time), sep="")
names(GSPvalue) = paste("p.GS.", names(Time), sep="")

### 6.4 模内分析：鉴定具有高GS和MM的基因 ----
module = "royalblue"                 # 选择模块
column = match(module, modNames)     # 匹配模块名称
moduleGenes = moduleColors==module   # 提取模块内基因
table(moduleGenes)

pdf("your pathway.pdf",
    width = 6,height = 7)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Time",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "red")
dev.off()

### 6.5 提取指定模块的基因做PPI网络分析 ----

adjacency = adjacency(datExpr, power =6)  # 计算邻接矩阵
TOM = TOMsimilarity(adjacency)             # 计算拓扑重叠矩阵（TOM）

module="royalblue"                               # 选择要导出的模块
probes = colnames(datExpr)                  # 获取基因名称
inModule = (moduleColors==module)           # 找到属于当前模块的基因
modProbes=probes[inModule]                  # 提取属于当前模块的基因名称
head(modProbes)                             # 显示基因名称前几行
modTOM = TOM[inModule,inModule]             # 提取属于当前模块的基因之间的TOM值
dimnames(modTOM)=list(modProbes,modProbes)  # 修改维度名称

# 这里只选了top100的基因
nTop=165                                       # 设置要选择的基因数目
IMConn = softConnectivity(datExpr[,modProbes])  # 计算当前模块中基因之间的相似性
top=(rank(-IMConn)<=nTop)                      # 找到相似性排名前nTop的基因
filterTOM=modTOM[top,top]                      # 提取相似性排名前nTop的基因之间的TOM值
# for visANT
vis = exportNetworkToVisANT(filterTOM,
                            file = paste("visANTinput-",module,".txt",sep = ""),
                            weighted = T,threshold = 0)  

# for cytoscape
cyt = exportNetworkToCytoscape(filterTOM,
                               edgeFile = paste("your pathway/CytoscapeInput-edges-", 
                                                paste(module, collapse="-"), 
                                                ".txt", sep=""),
                               nodeFile = paste("your pathway/CytoscapeInput-nodes-", 
                                                paste(module, collapse="-"), 
                                                ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               altNodeNames = 
                                 nodeNames = modProbes[top], 
                               nodeAttr = moduleColors[inModule][top])  
# 7. KEGG GO -------------------------------------------------------------------
## 7.1 Set output catagory----
# 指定文件夹路径
dir.create("your pathway")
dir_WGCNA <- "your pathway/"

## 7.2 WGCNA_res input ----
WGCNA_gene <- read.csv('./')

# 转换基因名 
y <- WGCNA_gene$Genes
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
WGCNA_gene$gene <- gene

## 7.3 Module genes ----
Module_genes <- subset(WGCNA_gene, WGCNA_gene$module == "turquoise")

if(T){
  # 设置数据库 
  GO_database <- 'org.Hs.eg.db'  # GO是org.Hs.eg.db数据库
  KEGG_database <- 'hsa'         # KEGG是hsa数据库
  
  # gene ID转换 
  gene <- bitr(Module_genes$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)
  
  ## 7.4 GO ----
  # GO富集分析
  GO <- enrichGO(gene = gene$ENTREZID,  # 导入基因的ENTREZID编号
                  OrgDb = GO_database,   # 用到的数据库（人类是：org.Hs.eg.db）
                  keyType = "ENTREZID",  # 设定读取的gene ID类型
                  ont = "ALL",           # (ont为ALL因此包括 BP,CC,MF三部分）
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2)    # 设定q值阈值
  
  ## 7.5 KEGG ----
  # KEGG富集分析
  KEGG <- enrichKEGG(gene = gene$ENTREZID,
                   keyType = "kegg",
                   organism = KEGG_database,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2)
  
  ## GO、KEGG结果整合 
  result <- list(enrichGO = GO, enrichKEGG = KEGG)
  GO_res <- result$enrichGO
  KEGG_res <- result$enrichKEGG
  
  ## 7.6  res_output ----
  # 导出enrichGO 
  write.csv(GO_res@result, file = paste0(dir_WGCNA, "/GO_res.csv"), quote = F, row.names = F)
  
  # dotplot
  pdf(file = paste0(dir_WGCNA, "/GO.pdf"), width = 6, height = 7)
  p1 <- dotplot(GO_res, showCategory = 5, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free') 
  print(p1)
  dev.off()
  
  # 导出enrichKEGG
  write.csv(KEGG_res@result, file = paste0(dir_WGCNA, "/KEGG_res.csv"), quote = F, row.names = F)
  
  # dotplot
  pdf(file = paste0(dir_WGCNA, "/KEGG.pdf"), width = 6, height = 5)
  p2 <- dotplot(KEGG_gene,showCategory = 10)
  print(p2)
  dev.off()
}

## 7.7 统计通路的数量 ----
# GO 
table(GO_res@result$p.adjust<0.05)
# KEGG 
table(KEGG_res@result$p.adjust<0.05)
