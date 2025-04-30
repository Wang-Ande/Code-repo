# 设置工作路径 ----
setwd("../AML_project/aml_analysis")
# 载入包 ----
source("QC_heatmap.R")
library(readxl)
library(openxlsx)
library(ggplot2)
library(readr)
library(dplyr)
# 导入分组信息 ----
data_group <- read_excel("./aml数据/MOLM13_group.xlsx")
table(data_group$group)
data_group <- as.data.frame(data_group)
# 发现列名当中_被转换为-,使用正则表达式转换回来
#data_group$id <- gsub('-', '_', data_group$id)
head(data_group,n=6) #检验是否转换为正确格式
# 配色设置 ----
value_colour <- c("MOLM13_WT" = "#E64B35FF",# Experimental group
                  "MOLM13_D" = "#4DBBD5FA" # other group1
                  ) # "MOLM13" "MV4_11" "OCI_AML2"
rownames(data_group) <- data_group$id
# RNA-Seq matrix input ---- 
data_input <- read_csv("./04_DEseq2_nordata_output/MOLM13_svabatch_rld.csv")
head(data_input,n=6) #检验是否读取正确
# 若文件格式有错误则需要调整，若无则直接跳转到 as.data.frame 转换
# 发现列名当中-被转换为.,使用正则表达式转换回来
colnames(data_input) #显示当前colnames
#colnames(data_input) <- gsub('-', '_', colnames(data_input))
#head(data_input,n=6) #检验是否转换为正确格式
#分隔线-----------------------------------------------------
data_input <- as.data.frame(data_input)
names(data_input)[1]<-"gene_id"
rownames(data_input) <- data_input$gene_id
data_input<-data_input[,-1]
# 处理RNA-Seq数据，将低表达的基因删除
# 使用 apply 函数
#data_before<-data_input[apply(data_input,1,sum)>15,]
# 创建输入路径 ----
dir <- "04_DEseq2_nordata_output/"
# heatmap -----------------------------------------------------------------
# 进行heatmap图绘制时，需要去除标准差为零的基因，否则会出现缺失值
# 去除标准差不等于0后的结果赋值为data_before_1
#data_before_1<-data_before[apply(data_before,1,sd)!=0,]

pdf(file = paste0(dir,"MOLM13_heatmap_svabatch.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = data_input,data_group = data_group,
           value_colour = value_colour)
dev.off()

# 多蛋白表达热图 ----
expr <- read.csv("./01_Data/report.pg_matrix_fill_norma.csv",row.names = 1)
anno <- read.xlsx("./01_Data/data_anno.xlsx",rowNames = TRUE)
expr_anno <- merge(expr, anno, by.x = 0, by.y = 0, all.x = TRUE)

# 转换基因名 
y <- expr_anno$Genes
gene1 <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
expr_anno$gene <- gene1

# 选择目标蛋白
targeted_prote <- expr_anno[grep("PHGDH|PSAT1|PSPH|SLC1A4|SLC1A5|TP53|FLT3|CD38|CDK1|CDK2|CDK4|CDK6|CD69|CDKN1A|CDKN1B|WEE1|PKMYT1|E2F1|CDC25|RB1|MYC|PLK1|PCNA|MCM2|ENSA|MKI67|BCAT1", expr_anno$gene),]
rownames(targeted_prote) <- targeted_prote$gene
targeted_prote <- targeted_prote[,grep("OCI", colnames(targeted_prote))]

rownames(targeted_prote)
targeted_prote <- targeted_prote[-grep("TP53I11|TP53BP1|TP53BP2|TP53I3|TP53RK|CDK11B|CDK13|CDK10|CDK19|CDK12|CDK11A|CDK2P1|URB1|ARRB1|ADARB1|GRB1|CCDC25|LILRB1|ZCRB1|RB1CC1|SCARB1",rownames(targeted_prote)),]
targeted_prote <- targeted_prote[,order(colnames(targeted_prote))]
colnames(targeted_prote)
Sample_order <- c("OCI_WT_1", "OCI_WT_2", "OCI_WT_3","OCI_2W_1", "OCI_2W_2", "OCI_2W_3", "OCI_4W_2","OCI_6W_1", "OCI_6W_2", "OCI_6W_3")
Sample_order <- gsub("OCI","MV4_11",Sample_order)
targeted_prote <- targeted_prote[,Sample_order]

library(circlize)
library(ComplexHeatmap)
expr_matrix <- targeted_prote
expr_matrix <- targeted_prote[,-grep("4W",colnames(targeted_prote))]
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

# 创建 ComplexHeatmap 的列注释对象
ha_col <- HeatmapAnnotation(
  Group = sample_groups,
  col = list(Group = group_colors),
  annotation_name_side = "right",
  annotation_legend_param = list(title = "Group", 
                                 title_gp = gpar(fontface = "bold",fontsize = 10),
                                 labels_gp = gpar(fontsize = 8))
)

log_expr <- log2(expr_matrix + 1)
scaled_expr <- t(scale(t(log_expr)))  # 每行（每个基因）标准化

ht <- Heatmap(
  scaled_expr,
  name = "Z-score",
  top_annotation = ha_col,          # 加上分组注释
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
cairo_pdf("./03_Result/DEP/OCI_AML2/Proliferation_proteins_expr_heatmap.pdf", width = 6, height = 5)
ht
dev.off()
