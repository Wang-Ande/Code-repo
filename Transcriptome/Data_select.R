# 清空环境变量 ----
rm(list=ls())
# 设置工作路径 ----
setwd("./01_data/")
# 创建输出路径 ----
dir.create("aml数据_time")
dir <- "./01_data/fpkm_data_by_time/"
# 加载包 ----
library(dplyr)
library(openxlsx)
library(readr)
# 数据导入 ----
## 表达矩阵 ----
input_counts <- read.csv("./01_data/Fpkm_data/data_merge_adjusted_fpkm.csv")
## 分组矩阵 ----
input_group <- read.xlsx("./01_data/Counts_data/aml_group_merge.xlsx")
# 分拣数据 ----
## 表达矩阵 ----
counts_1 <- input_counts[,c(2,17:19,24:26)] # WT_2W
counts_2 <- input_counts[,c(2,17:19,21,22,29)] # 2W_6W
counts_3 <- input_counts[,c(2,21:25,29)] # WT_6W
## 分组矩阵 ----
group_1 <- input_group[c(15:17,21:23),] # WT_2W
group_2 <- input_group[c(15:17,19,20,27),] # 2W_6W
group_3 <- input_group[c(19:23,27),] # WT_6W
# 结果导出 ----
#OCI_AML2/MV4_11/MOLM13
## 数据一
write.xlsx(counts_1,file = paste0(dir,"OCI_AML2_fpkm_sva_WT_2W.xlsx"))
write.xlsx(group_1,file = paste0(dir,"OCI_AML2_group_WT_2W.xlsx"))
## 数据二
write.xlsx(counts_2,file = paste0(dir,"OCI_AML2_fpkm_sva_2W_6W.xlsx"))
write.xlsx(group_2,file = paste0(dir,"OCI_AML2_group_2W_6W.xlsx"))
## 数据三
write.xlsx(counts_3,file = paste0(dir,"OCI_AML2_fpkm_sva_WT_6W.xlsx"))
write.xlsx(group_3,file = paste0(dir,"OCI_AML2_group_WT_6W.xlsx"))

