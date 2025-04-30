library(sva)
# input ----
## gene expression input ----
library(sva)
library(readr)
data_batch1 <- read_csv("01_rawdata/gene_count_matrix_batch1.csv")
data_batch2 <- read_csv("01_rawdata/gene_count_matrix_batch2.csv")
colnames(data_batch2) <- gsub("-","_",colnames(data_batch2))
data_merge <- merge(data_batch1,data_batch2,by = "gene_id")
data_merge <- as.data.frame(data_merge)
rownames(data_merge) <- data_merge$gene_id
data_merge <- subset(data_merge,select = -c(gene_id))
## group input ----
library(readxl)
data_group_all <- read_excel("01_rawdata/aml_group_merge.xlsx")
colnames(data_group_all)
data_group <- subset(data_group_all,select = c(id,batch))
colnames(data_group)[2] <- "group"
value_colour <- value_colour <- c("experiment" = "#E64B35FF",# Experimental group
                                  "control" = "#4DBBD5FF"# other group1
)
value_colour <- value_colour <- c("batch1" = "#E64B35FF",# Experimental group
                                  "batch2" = "#4DBBD5FF"# other group1
)
value_colour <- value_colour <- c("MOLM13" = "#E64B35FF",# Experimental group
                                  "MV4_11" = "#4DBBD5FF",# other group1
                                  "OCI_AML2" = "#4DBAFF"
)
data_merge <- data_merge[apply(data_merge,1,sd) != 0,]
# boxplot -----------------------------------------------------------------
pdf(file = paste0(dir,"QC_boxplot_before.pdf"),
    width = 6,
    height = 4)

QC_boxplot(data_merge,data_group = data_group,
           value_colour = value_colour,
           title = "raw data")
dev.off()

# heatmap -----------------------------------------------------------------

pdf(file = paste0(dir,"QC_heatmap_before.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = data_merge,data_group = data_group,
           value_colour = value_colour)
dev.off()

# PCA ---------------------------------------------------------------------
pdf(file = paste0(dir,"QC_pca_before.pdf"),
    width = 5,
    height = 5)
QC_PCA(data = data_merge,
       data_group = data_group,
       value_colour = value_colour)
dev.off()

data_group_all$biological <- paste0(data_group_all$group,data_group_all$cell)
mod <- model.matrix(~as.factor(biological), data=data_group_all)
rownames(mod) <- data_group_all$id
data_merge_adjusted <- sva::ComBat_seq(counts = as.matrix(data_merge), batch=data_group_all$batch,
                                       group=data_group_all$biological)
data_merge_adjusted <- data_merge_adjusted[apply(data_merge_adjusted,1,sd) != 0,]
write.csv(data_merge_adjusted,file = "01_rawdata/data_merge_adjusted.csv")

QC_PCA(data = data_merge_adjusted,
       data_group = data_group,
       value_colour = value_colour)
QC_heatmap(data = data_merge_adjusted,data_group = data_group,
           value_colour = value_colour)

# limma ----
source("02_code/run_DE.R")
# 注意，data和data_anno的行名应一致
# 根据分组选择要进行差异分析的组别
table(data_group_all$cell)
data_group <- subset(data_group_all,data_group_all$cell%in%c("OCI_AML2"))
data_group <- subset(data_group,select = c(id,group))

data_merge_adjusted_select <- data_merge_adjusted[,colnames(data_merge_adjusted)%in%data_group$id]
data_merge_adjusted_select <- data_merge_adjusted_select[apply(data_merge_adjusted_select,1,sd) != 0,]
# group 1为实验组
# group 2为对照组
group_1 <- "experiment"
group_2 <- "control"
y <- rownames(data_merge_adjusted_select)
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),"[|]")[[1]][2]))
data_anno <- as.data.frame(gene)
rownames(data_anno) <- y

result_merge <- run_DE(data = data_merge_adjusted_select,
                       data_group = data_group,
                       log2 = T,
                       data_anno = data_anno,
                       group_1 = group_1,group_2 = group_2,
                       dir = "03_result/DE/OCI_AML2/")

# volcano plot ------------------------------------------------------------
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)

dif_logFC_up <- subset(result_merge,result_merge$logFC > 0&result_merge$P.Value < 0.05)
dif_logFC_up <- dif_logFC_up[order(dif_logFC_up$logFC,decreasing = T),]
dif_logFC_down <- subset(result_merge,result_merge$logFC < 0&result_merge$P.Value < 0.05)
dif_logFC_down <- dif_logFC_down[order(dif_logFC_down$logFC),]

# y <- c(rownames(dif_logFC_up)[1:10],rownames(dif_logFC_down)[1:10])
# y <- na.omit(y)

res_data <- result_merge
data <- res_data[res_data[,"P.Value"] <= 1,]

y <- result_merge$Genes
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
data$gene <- gene
# 颜色划分padj <0.05，且log2FC绝对值大于sd(tempOutput$logFC)*3
data$sig[data$P.Value >= 0.05 | abs(data$logFC) < 0] <- "Not"

data$sig[data$P.Value < 0.05 & data$logFC >= 0] <- "Up"

data$sig[data$P.Value < 0.05 & data$logFC <= -0] <- "Down"
input <- data
library(ggrepel)
library(ggplot2)
volc <- ggplot(data = data, aes(x = logFC,
                                y = -log10(P.Value),
                                color = sig)) +
    geom_point(alpha = 0.9) +  theme_classic() +
    theme(panel.grid = element_blank(),strip.text = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = c("#4DBBD5","grey80","#E64B35")) +
    geom_hline(yintercept = -log10(0.05),lty = 4,lwd = 0.6,alpha = 0.8) +
    geom_vline(xintercept = c(-1,1),lty = 4,lwd = 0.6,alpha = 0.8) +

    labs(title = paste0(group_1,"-",group_2))
volc
ggsave(plot = last_plot(),filename = paste0("03_result/DE/OCI_AML2/","volc.pdf"),
       height = 5,
       width = 6)

