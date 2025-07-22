rm(list=ls())
library(matrixTests)
library(tidyr)
library(dplyr)
library(data.table)
library(limma)
library(agricolae)
library(graphics)
library(preprocessCore)
library(ggplot2)
library(ggpubr)
setwd("E:\\work")

#输入的数据为删除基因名重复后的表格
a<-xlsx::read.xlsx("2-tidy.xlsx",sheetIndex = 1,row.names = 1,header = T,check.names = F)
names(a)

#---------------按分组排序样本--------------
SampleType <- c("Ctrl","Allprep")

expr <- list()
for (i in 1:2) {
expr[[i]] <- a[,which(colnames(a) %like% SampleType[i])]
}
expr <- data.frame(do.call(cbind,expr))
names(expr)

#查看数据定量的批次
library(RColorBrewer)
cols <- brewer.pal(8, "Set1")
boxplot(log2(expr),col=cols,las=2,outline=F,ylab="log2(intensity)",frontsize=8)


#---------------median normalization---------------
#将所有样本的中位数与第一个样本的对齐
nor_facter<-apply(expr,2,median,na.rm = T)/apply(expr,2,median,na.rm = T)[1]  
nor_facter
expr<-t(t(expr)/nor_facter)
write.csv(expr,"3-median norm_data.csv")

expr <- log2(expr)
write.csv(expr,"4-log2 trans.csv")
boxplot(expr,col=cols,las=2,outline=F,ylab="log2(intensity)",frontsize=8)

#--------------------筛选可定量蛋白---------------
expr <- as.matrix(expr)
miss_value <- nrow(data.frame(which(is.na(expr))))*100/(nrow(expr)*ncol(expr))
miss_value

#统计每个蛋白在每组中的缺失比例
index <- list()
for (i in 1:nrow(data.frame(SampleType))) {
  num <- which(colnames(expr) %like% SampleType[i])
  index[[i]] <- apply(expr[,num], 1, function(x) sum(is.na(x)))*100/nrow(data.frame(num))
}

index1 <- data.frame(do.call(cbind,index))
colnames(index1) <- SampleType


#过滤至少有一组缺失值不高于33.3%的蛋白
expr1 <- data.frame(expr[which(apply(index1, 1, function(x) sum(x<34))>0),])
nrow(expr);nrow(expr1)
write.csv(expr1,"5-all screen PG.csv")

#---------------------缺失值填补------------------
#在至少一组中缺失值不高于33.3%,另一组缺失值不低于66.7%的蛋白则使用最小值填补
screen1 <- data.frame(expr[which(index1$Ctrl<34 & index1$Allprep>66|
                                 index1$Ctrl>66 & index1$Allprep<34),])
screen1 <- screen1[intersect(row.names(expr1),row.names(screen1)),]
screen1[is.na(screen1)] <- min(expr,na.rm=T)
screen1_tidy <- data.frame(mean_Ctrl=rowMeans(screen1[,1:3],na.rm = T),
                           mean_Allprep=rowMeans(screen1[,4:6],na.rm = T))
screen1_tidy$logFC <- screen1_tidy$mean_Allprep-screen1_tidy$mean_Ctrl
screen1_tidy$group <- ifelse(screen1_tidy$logFC>0,"Allprep","Ctrl")
screen1 <- data.frame(screen1_tidy,screen1)
write.csv(screen1,"6-unique PG.csv")

###对剩余的蛋白进行KNN填补
library(DMwR2)
screen2 <- knnImputation(expr1[setdiff(row.names(expr1),row.names(screen1)),])
write.csv(screen2,"7-shared PG.csv")

#----------------------limma差异分析-----------------
## 构建分组
m1=SampleType[1]
m2=SampleType[2]
indexm1 <- which(colnames(screen2) %like% m1);indexm1
indexm2 <- which(colnames(screen2) %like% m2);indexm2
n1=nrow(as.data.frame(indexm1));n1
n2=nrow(as.data.frame(indexm2));n2
df1=screen2[,c(indexm1,indexm2)]

##构建design
group_list <- factor(c(rep(m1,n1), rep(m2,n2)))
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(df1)
contrastIndex <- paste0(m2,"-",m1)
contrast.matrix<-makeContrasts(contrastIndex,levels = design)
contrast.matrix 

##差异分析
fit <- lmFit(df1, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res_limma = topTable(fit2, n=Inf,adjust.method="BH")
res_limma <- na.omit(res_limma)[,c(1,4,5)]
colnames(res_limma)[2:3] <- c("pvalue","FDR")
means <- data.frame(mean_Ctrl=rowMeans(df1[,c(indexm1)],na.rm = T),
                    mean_Allprep=rowMeans(df1[,c(indexm2)],na.rm = T))
limma_expr=merge(res_limma,means,by="row.names",all=F)
limma_expr2=merge(limma_expr,df1,by.x="Row.names",by.y="row.names",all=F)
limma_sig=limma_expr2[which(limma_expr2$FDR<0.05 & abs(limma_expr2$logFC)>log2(1.5)),]
dim(limma_sig)

colnames(limma_expr2)[1] <- "Gene"
colnames(limma_sig)[1] <- "Gene"
write.csv(limma_expr2,paste0("8-",m2," vs ",m1," shared_expr.csv"),row.names = F)
write.csv(limma_sig,paste0("9-",m2," vs ",m1," shared_DEPs.csv"),row.names = F)


