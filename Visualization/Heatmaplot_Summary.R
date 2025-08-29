# 设置工作路径 
setwd("../AML_project/aml_analysis")
# 载入包 
source("QC_heatmap.R")
library(readxl)
library(openxlsx)
library(ggplot2)
library(readr)
library(dplyr)
# 导入分组信息
data_group <- read_excel("./aml数据/MOLM13_group.xlsx")
table(data_group$group)
data_group <- as.data.frame(data_group)
# 发现列名当中_被转换为-,使用正则表达式转换回来
#data_group$id <- gsub('-', '_', data_group$id)
head(data_group,n=6) #检验是否转换为正确格式
# 配色设置
value_colour <- c("MOLM13_WT" = "#E64B35FF",# Experimental group
                  "MOLM13_D" = "#4DBBD5FA" # other group1
                  ) # "MOLM13" "MV4_11" "OCI_AML2"
rownames(data_group) <- data_group$id
# RNA-Seq matrix input
data_input <- read_csv("./04_DEseq2_nordata_output/MOLM13_svabatch_rld.csv")
head(data_input,n=6) #检验是否读取正确
# 若文件格式有错误则需要调整，若无则直接跳转到 as.data.frame 转换
# 发现列名当中-被转换为.,使用正则表达式转换回来
colnames(data_input) #显示当前colnames
#colnames(data_input) <- gsub('-', '_', colnames(data_input))
#head(data_input,n=6) #检验是否转换为正确格式

# ---------------------------- 1. 常规热图 -------------------------------------
data_input <- as.data.frame(data_input)
names(data_input)[1]<-"gene_id"
rownames(data_input) <- data_input$gene_id
data_input<-data_input[,-1]
# 处理RNA-Seq数据，将低表达的基因删除
# 使用 apply 函数
#data_before<-data_input[apply(data_input,1,sum)>15,]
# 创建输入路径
dir <- "04_DEseq2_nordata_output/"

# 进行heatmap图绘制时，需要去除标准差为零的基因，否则会出现缺失值
# 去除标准差不等于0后的结果赋值为data_before_1
#data_before_1<-data_before[apply(data_before,1,sd)!=0,]
pdf(file = paste0(dir,"MOLM13_heatmap_svabatch.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = data_input,data_group = data_group,
           value_colour = value_colour)
dev.off()

# ------------------------- 2. 基因标注表达热图 --------------------------------
expr <- read.csv("./01_Data/report.pg_matrix_fill_norma.csv",row.names = 1)
anno <- read.xlsx("./01_Data/data_anno.xlsx",rowNames = TRUE)
expr_anno <- merge(expr, anno, by.x = 0, by.y = 0, all.x = TRUE)

# 转换基因名 
y <- expr_anno$Genes
gene1 <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
expr_anno$gene <- gene1

# 选择目标蛋白
targeted_prote <- expr_anno[grep("PHGDH|PSAT1|PSPH|SLC1A4|SLC1A5|TP53|FLT3|CD38|CXCR4|PCNA|AKT|CDK1|ENSA|MKI67|BCAT1", expr_anno$gene),]
rownames(targeted_prote) <- targeted_prote$gene
targeted_prote <- targeted_prote[,grep("MOLM13|MV4_11|OCI", colnames(targeted_prote))]

rownames(targeted_prote)
targeted_prote <- targeted_prote[-grep("TP53I11|TP53BP1|TP53BP2|TP53I3|TP53RK|CDK11B|CDK13|CDK10|CDK19|CDK12|CDK11A|CDK2P1|URB1|ARRB1|ADARB1|GRB1|CCDC25|LILRB1|ZCRB1|RB1CC1|SCARB1",rownames(targeted_prote)),]
targeted_prote <- targeted_prote[,order(colnames(targeted_prote))]
colnames(targeted_prote)
Sample_order <- c("MOLM13_WT_1", "MOLM13_WT_2", "MOLM13_WT_3", 
                  "MOLM13_2W_2", "MOLM13_2W_3", 
                  "MOLM13_6W_1", "MOLM13_6W_2", "MOLM13_6W_3",
                  "MV4_11_WT_1", "MV4_11_WT_2", "MV4_11_WT_3",
                  "MV4_11_2W_1", "MV4_11_2W_2", "MV4_11_2W_3", 
                  "MV4_11_6W_1", "MV4_11_6W_2", "MV4_11_6W_3",
                  "OCI_WT_1", "OCI_WT_2", "OCI_WT_3",
                  "OCI_2W_1", "OCI_2W_2", "OCI_2W_3", 
                  "OCI_4W_2","OCI_6W_1", "OCI_6W_2", "OCI_6W_3")
# Sample_order <- gsub("OCI","MV4_11",Sample_order)
targeted_prote <- targeted_prote[,Sample_order]

library(circlize)
library(ComplexHeatmap)
expr_matrix <- targeted_prote
expr_matrix <- targeted_prote[,-grep("4W",colnames(targeted_prote))]

# 定义 Cell line 分组
cell_lines <- case_when(
  grepl("MOLM13", colnames(expr_matrix)) ~ "MOLM13",
  grepl("MV4_11", colnames(expr_matrix)) ~ "MV4_11",
  grepl("OCI", colnames(expr_matrix)) ~ "OCI_AML2"
)
cell_lines <- factor(cell_lines, levels = c("MOLM13", "MV4_11", "OCI_AML2"))

# 为 Cell line 设置颜色
cell_line_colors <- c(
  MOLM13 = "#66C2A5",
  MV4_11 = "#8DA0CB",
  OCI_AML2 = "#FC8D62"
)

# 创建分组信息
sample_groups <- case_when(
  grepl("WT", colnames(expr_matrix)) ~ "WT",
  grepl("2W|4W", colnames(expr_matrix)) ~ "Low",
  grepl("6W", colnames(expr_matrix)) ~ "High"
)

# 转换为因子，设置顺序
sample_groups <- factor(sample_groups, levels = c("WT", "Low", "High"))

# 创建颜色映射
group_colors <- c(WT = "#6388B4", Low = "#FFAE34", High = "#EF6F6A")

# 创建包含两个分组注释的列注释对象（Cell line 在上面）
ha_col <- HeatmapAnnotation(
  Cell_line = cell_lines,
  Group = sample_groups,
  col = list(
    Cell_line = cell_line_colors,
    Group = group_colors
  ),
  annotation_name_side = "right",
  annotation_legend_param = list(
    Cell_line = list(title = "Cell line", title_gp = gpar(fontface = "bold", fontsize = 10), labels_gp = gpar(fontsize = 8)),
    Group = list(title = "Group", title_gp = gpar(fontface = "bold", fontsize = 10), labels_gp = gpar(fontsize = 8))
  ),
  annotation_height = unit(c(6, 6), "mm")  # 控制两个注释行的高度
)

# 创建 group 的列注释对象
ha_col <- HeatmapAnnotation(
  Group = sample_groups,
  col = list(Group = group_colors),
  annotation_name_side = "right",
  annotation_legend_param = list(title = "Group", 
                                 title_gp = gpar(fontface = "bold",fontsize = 10),
                                 labels_gp = gpar(fontsize = 8))
)

# scale
log_expr <- log2(expr_matrix + 1)
scaled_expr <- t(scale(t(log_expr)))  # 每行（每个基因）标准化

# plot
ht <- Heatmap(
  scaled_expr,
  name = "Z-score",
  top_annotation = ha_col,           # 加上分组注释
  col = circlize::colorRamp2(c(-2, 0, 2), colors = c("#1E90FF", "white", "#FF4500")),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = TRUE,
  show_row_names = TRUE,
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title = "Z-score", 
                              title_gp = gpar(fontface = "bold",fontsize = 10),
                              legend_height = unit(3, "cm"),   # 控制图例高度
                              grid_width = unit(0.35, "cm"),    # 色块宽度
                              labels_gp = gpar(fontsize = 8)   # 标签字体
                              ))
ht <- draw(ht,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE,
     padding = unit(c(2, 12, 2, 2), "mm"),  # 调整边距
     legend_gap = unit(6, "mm")) # 图例边距
# res output
cairo_pdf("./03_Result/DEP/Targeted_proteins_expr_heatmap.pdf", width = 7, height = 5)
ht
dev.off()

# ------------------------- 3. 聚类分析标注热图 --------------------------------
rm(list=ls())
## 加载R包
library(ComplexHeatmap)
library(shinipsum)
library(ggplot2)
library(ggstatsplot)
library(patchwork)
library(reshape2)
library(stringr)
library(limma)
library(tidyverse)
getOption('timeout')
options(timeout=10000)

# 导入数据
data <- readxl::read_excel("./01_data/BLOOD_example.xlsx",skip = 3)
dim(data)
data[1:6,1:6]

# 标记的基因
marker <- list(cluster1=c("CTSA","SMPD1"), cluster2=c("CD14","CD86","FGR","TLR8","ITGAM"),
               cluster3=c("MEIS1","PBX3","MEF2C"),
               cluster4=c("MYC","LRP5","RUNX3"))
marker

## 修饰每一个文字注释条的参数
tcol <- c("black","black","#34b3e6","#f46951")
bg<- c("#e8f8f8","#fdefe9","#e9f6fb","#fde8e9")

for(i in 1:length(marker)) {
  #i <- 1
  marker[[i]] <- data.frame(text=marker[[i]],fontsize=12,col=tcol[i],
                            fontface="italic")
}
marker
str(marker)

# 注释通路
annotation <- list(cluster1=c("Lysosome"),
                   cluster2=c("CD14 Positive Monocyte","CD33 Positive Myeloid","Regulation of Inflammatory Response",
                              "Hematopoietic cell lineage","Increased Monocyte Cell number"),
                   cluster3=c("KMT2A Fusion Target Genes","MLL-AF6-Spreading Target Genes"),
                   cluster4=c("Ribosome Biogenesis","Myc Targets V2"))

annotation
for(i in 1:length(annotation)) {
  #i <- 1
  annotation[[i]] <- data.frame(text=annotation[[i]],fontsize=11.3,col="black")
}
annotation

# 注释pvalue
adj_p_values <- list(cluster1=c(1.862e-9),
                     cluster2=c(6.094e-55, 9.956e-47, 2.308e-12, 3.595e-8, 5.908e-5),
                     cluster3=c(4.557e-11, 1.715e-7),
                     cluster4=c(1.030e-16, 2.745e-14))

adj_p_values
for(i in 1:length(adj_p_values )) {
  #i <- 1
  adj_p_values [[i]] <- data.frame(text=as.character(adj_p_values[[i]]),fontsize=11.3,col="black")
}
adj_p_values

## 绘图
## 热图与文本框链接
mat <- as.matrix(data[,-c(1:3)])
rownames(mat) <- data$gene_id
dim(mat)
head(mat)
# 表达值归一化
mat <- t(scale(t(mat)))
colnames(mat) <- str_split(colnames(mat), "_rep",simplify = T)[,1]


# 行分割
split <- data$cluster
split

ha <- rowAnnotation(textbox1 = anno_textbox(split, marker,          # split 与 text的名字对应
                                            add_new_line = TRUE,    # 每一个word一行
                                            word_wrap = TRUE,       # 控制单词换行和新行：
                                            line_space = unit(3,"mm"),
                                            background_gp = gpar( fill=bg,col = NA), #：设置文本框的背景和边框颜色
                                            by = "anno_block"
))


#
ha1 <- rowAnnotation(textbox2 = anno_textbox(split, annotation,      # split 与 text的名字对应
                                             add_new_line = TRUE,    # 每一个word一行
                                             word_wrap = TRUE,       # 控制单词换行和新行：
                                             line_space = unit(3,"mm"),
                                             background_gp = gpar( fill=bg,col = NA), #：设置文本框的背景和边框颜色
                                             by = "anno_block"
))

#
ha2 <- rowAnnotation(textbox3 = anno_textbox(split, adj_p_values,    # split 与 text的名字对应
                                             add_new_line = TRUE,    # 每一个word一行
                                             word_wrap = TRUE,       # 控制单词换行和新行：
                                             line_space = unit(3, "mm"),
                                             background_gp = gpar( fill=bg,col = NA), #：设置文本框的背景和边框颜色
                                             by = "anno_block",
))

# 文本框注释条在右边
pdf(file = "Fig4a.pdf",height = 7,width = 11)
ht_list = Heatmap(mat, name = "mat", cluster_rows = FALSE, row_split = split,
                  show_row_names = F,  # 去掉行名
                  cluster_columns = F, # 列不聚类
                  row_title_rot = 0,   # 右边的标题水平放着
                  row_title_gp = gpar(fontsize = 12),
                  heatmap_legend_param = list(direction = "horizontal",nrow = 1), # 图例水平放置
                  column_names_side = c("top"),
                  right_annotation = ha
) +
  ha1 + ha2

draw(ht_list,heatmap_legend_side = "bottom")# 图例放在底部

# 给三个文本注释框添加标题
decorate_annotation("textbox1",{
  grid.text("Genes", y = unit(1,"npc") + unit(2,"mm"), just = "bottom")
})
decorate_annotation("textbox2", {
  grid.text("Annotation", y = unit(1,"npc") + unit(2,"mm"), just =c("right","bottom"))
})
decorate_annotation("textbox3", {
  grid.text("Adj.Pvalue", y = unit(1,"npc") + unit(2,"mm"), just ="bottom")
})
dev.off()
