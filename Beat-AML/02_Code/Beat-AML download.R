# Beat-AML download
# packages
library(devtools)
devtools::install_github("radivot/AMLbeatR",subdir="AMLbeatR")
library(AMLbeatR)
library(tidyverse)

# load data
mkBeatAML(beatHome = "./Beat-AML/Beat-AML1.0",
          beatFile = "41586_2018_623_MOESM3_ESM.xlsx",outFile ="beatAML")
load("./Beat-AML/Beat-AML1.0/beatAML.RData")

# clinical
surv = tidyClin(clin)                        # 修改为行为患者
surv = muts(surv, v, av, n=3)                # 整理突变信息
surv = as.data.frame(surv)                   # 变成数据框
rownames(surv) = surv$lid
surv = surv[,c("status","surv","age","sex")] # 提取状态、时间、年龄、性别
surv = na.omit(surv)                         # 去除NA值
colnames(surv) = c("event","time","age","gender")
surv$time = surv$time/30                     # 时间按月为单位
surv = surv[surv$time>1,]                    # 删除随访时间小于一个月的
surv = surv[rownames(surv)%in%colnames(cpm),] # 保留有表达矩阵的患者

# expression
exp = cpm
rownames(exp) = exp$Symbol
exp = exp[,rownames(surv)]                   # 保留有生存信息的患者表达矩阵 
exp = as.data.frame(t(exp))                  # 改为数据框格式

# merge
identical(rownames(surv), rownames(exp))    # 检查行名是否一致
surv = cbind(surv, exp)                     # 合并
BeatAML_surv = surv
dim(BeatAML_surv)

# res output
save(BeatAML_surv, file = "./Beat-AML/Beat-AML1.0/BeatAML_surv.Rdata")
write.csv(BeatAML_surv,file = "./Beat-AML/Beat-AML1.0/BeatAML_surv.csv")
