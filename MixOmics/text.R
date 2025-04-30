# text
library(mixOmics)
# load data ----
data(srbct)
X <- srbct$gene

# Outcome y that will be internally coded as dummy:
Y <- srbct$class 
dim(X); length(Y)

summary(B)

# Initial exploration with PCA ----
pca.srbct <- pca(X, ncomp = 3, scale = TRUE)

plotIndiv(pca.srbct, group = srbct$class, ind.names = FALSE,
          legend = TRUE, 
          title = 'SRBCT, PCA comp 1 - 2')

plsda.srbct <- plsda(X,Y, ncomp = 10)

# Number of components in PLS-DA ----
set.seed(30) # For reproducibility with this handbook, remove otherwise
perf.plsda.srbct <- perf(plsda.srbct, validation = 'Mfold', folds = 3, 
                         progressBar = FALSE,  # Set to TRUE to track progress
                         nrepeat = 10)         # We suggest nrepeat = 50

plot(perf.plsda.srbct, sd = TRUE, legend.position = 'horizontal')


final.plsda.srbct <- plsda(X,Y, ncomp = 3)

plotIndiv(final.plsda.srbct, ind.names = FALSE, legend=TRUE,
          comp=c(1,2), ellipse = TRUE, 
          title = 'PLS-DA on SRBCT comp 1-2',
          X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2')
plotIndiv(final.plsda.srbct, ind.names = FALSE, legend=TRUE,
          comp=c(2,3), ellipse = TRUE, 
          title = 'PLS-DA on SRBCT comp 2-3',
          X.label = 'PLS-DA comp 2', Y.label = 'PLS-DA comp 3')

# Classification Performance ----
set.seed(30) # 设定随机种子，保证结果可复现
perf.final.plsda.srbct <- perf(final.plsda.srbct, validation = 'Mfold', 
                               folds = 3, 
                               progressBar = FALSE, # 设为 TRUE 以显示进度条
                               nrepeat = 10) # 建议使用 50 次重复，提高稳定性
perf.final.plsda.srbct$error.rate$BER[, 'max.dist']
##      comp1      comp2      comp3 
## 0.53850543 0.25986413 0.04481884
perf.final.plsda.srbct$error.rate.class$max.dist
##         comp1      comp2      comp3
## EWS 0.2565217 0.08695652 0.08260870
## BL  0.8125000 0.51250000 0.00000000
## NB  0.3000000 0.37500000 0.04166667
## RMS 0.7850000 0.06500000 0.05500000

# Background Prediction ----
# 1.计算 PLS-DA 预测背景区域
background.max <- background.predict(final.plsda.srbct, 
                                     comp.predicted = 2,  # 计算前两个主成分的背景预测区域
                                     dist = 'mahalanobis.dist')   # "centroids.dist", "mahalanobis.dist" or "max.dist" 
# max.dist 适用于类别区分明显的数据（简单但可能误分类）。
# centroids.dist 适用于类内差异较大的情况（更稳定）。
# mahalanobis.dist 适用于复杂数据（但计算量较大）。

# 2.生成 PLS-DA 样本图，并添加预测背景
plotIndiv(final.plsda.srbct, comp = 1:2, group = srbct$class,
          ind.names = FALSE, title = 'mahalanobis.dist',
          legend = TRUE, background = background.max)
