# limma包去除差异分析----
# 导入分组信息 ----
group<-read.xlsx("./")
## 因子化(不是必要操作) ----
design<-model.matrix(~factor(group$group))
# 导入计数矩阵 ----
exprSet<-read.csv("./")
# 设置批次 ----
batch<-c(rep(1,x),rep(2,y)) #
# 去除批次效应 ----
adj_experSet<-removeBatchEffect(experSet,batch=batch,design<-design)
