library(dplyr)
library(openxlsx)
library(readxl)
neg_data1 <- read.xlsx("./01_data/MV4_11_2w.vs.MV4_11_WT_neg_Diff_order.xlsx")
pos_data2 <- read.xlsx("./01_data/MV4_11_2w.vs.MV4_11_WT_pos_Diff_order.xlsx")

neg_data1_up <- neg_data1[which(neg_data1$Up.Down == "up"),]
neg_data1_up <- neg_data1_up[,c(2,19)]
neg_data1_down <- neg_data1[which(neg_data1$Up.Down == "down"),]
neg_data1_down <- neg_data1_down[,c(2,19)]
pos_data2_up <- pos_data2[which(pos_data2$Up.Down == "up"),]
pos_data2_up <- pos_data2_up[,c(2,19)]
pos_data2_down <- pos_data2[which(pos_data2$Up.Down == "down"),]
pos_data2_down <- pos_data2_down[,c(2,19)]

data_up <- rbind(neg_data1_up,pos_data2_up)
data_down <- rbind(neg_data1_down,pos_data2_down)

colnames(data_up)[colnames(data_up) == "group"] <- "MV4_11_2w.vs.MV4_11_WT"
write.xlsx(data_up,file = "./01_data/MV4_11_2w.vs.MV4_11_WT_allup_Diff_order.xlsx")
colnames(data_down)[colnames(data_down) == "Up.Down"] <- "MV4_11_2w.vs.MV4_11_WT"
write.xlsx(data_down,file = "./01_data/MV4_11_2w.vs.MV4_11_WT_alldown_Diff_order.xlsx")
