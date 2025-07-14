library(dplyr)
library(readr)
library(openxlsx)
library(ggplot2)
# QC methods summary

# 直方图 ----
# 展示数据分布特点 可以用于过滤低表达等等
# 读入数据
expr_raw <- read.csv("./01_data/Cpm_data/data_merged_adjusted.csv")
colnames(expr_raw)
expr_raw <- as.data.frame(expr_raw)
rownames(expr_raw) <- expr_raw$X
expr_raw <- expr_raw[,-1]
expr_raw <- as.matrix(expr_raw)
expr_raw <- log2(expr_raw+1)
expr_mean <- apply(expr_raw,1,mean)

# 绘制直方图
pdf("./03_result/QC/CPM/All/Hitogram.pdf",width = 5,height = 6)
hist(expr_mean, 
     main = " Data Distribution with Density", 
     xlab = "log2(Cpm+1)", 
     ylab = "Frequency", 
     col = "dodgerblue3", 
     border = "black", 
     prob = T ) # 让y轴显示密度而非频数
lines(density(expr_mean), col = "firebrick3", lwd = 2)  # 添加红色的密度曲线
dev.off()

# 选择过滤阈值
expr_raw <- as.data.frame(expr_raw)
expr_filter <- expr_raw[apply(expr_raw,1,mean)>0,] # 23238
