# 3. GO KEGG --------------------------------------------------------------
## 3.1 Output category ----
dir_gokegg <- "./03_result/03_Enrichment/MV4/High_vs_WT/"

## 3.2 DE_res input ----
DE_result <- read.csv('./03_result/02_DE/MV4/High_vs_WT/DE.csv')

## 3.3 Cut P.Value ----
Genesymbol <- subset(DE_result, P.Value < 0.05)

## 3.4 Cut logfc ----
cutoff <- 1

# 转换基因名 
y <- Genesymbol$Row.names
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),"\\|")[[1]][2]))
id <- data.frame(gene)
Genesymbol <- cbind(id,Genesymbol)

# 设置数据库 
GO_database <- 'org.Hs.eg.db'  # GO是org.Hs.eg.db数据库
KEGG_database <- 'hsa'         # KEGG是hsa数据库

if(T){
  ## 3.5 down genes ----
  down_genes <- subset(Genesymbol, logFC < -cutoff)
  
  # gene ID转换 
  gene <- bitr(down_genes$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)
  
  ### 3.5.1 GO ----
  # GO富集分析
  kkd <- enrichGO(gene = gene$ENTREZID, # 导入基因的ENTREZID编号
                  OrgDb = GO_database, # 用到的数据库（人类是：org.Hs.eg.db）
                  keyType = "ENTREZID", # 设定读取的gene ID类型
                  ont = "ALL", 
                  pvalueCutoff = 1,  
                  pAdjustMethod = "BH",
                  qvalueCutoff = 1,
                  readable = T
  )
  
  ### 3.5.2 KEGG ----
  # KEGG富集分析
  kk <- enrichKEGG(gene = gene$ENTREZID,
                   keyType = "kegg",
                   organism = KEGG_database,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 1,
                   use_internal_data = T
  )
  # 设置geneID可读格式
  kk <- setReadable(kk, OrgDb='org.Hs.eg.db', keyType='ENTREZID')
  
  ## GO、KEGG结果整合 
  result <- list(enrichGO = kkd, enrichKEGG = kk)
  
  # 结果标记为下调 
  result_down <- result
  kkd_down <- result_down$enrichGO
  kk_down <- result_down$enrichKEGG
  
  ### 3.5.3 下调res_output ----
  # 导出下调enrichGO 
  write.xlsx(kkd_down@result, file = paste0(dir_gokegg, "GO_down.xlsx"), rowNames = FALSE)
  
  # dotplot
  pdf(file = paste0(dir_gokegg, "GO_down.pdf"), width = 5.5, height = 7)
  p1 <- dotplot(kkd_down, showCategory = 5, color = "pvalue", split = "ONTOLOGY") + 
    facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free') 
  print(p1)
  dev.off()
  
  # 导出下调enrichKEGG
  write.xlsx(kk_down@result, file = paste0(dir_gokegg, "KEGG_down.xlsx"), rowNames = FALSE)
  
  # dotplot
  pdf(file = paste0(dir_gokegg, "KEGG_down.pdf"), width = 6, height = 5)
  p2 <- dotplot(kk_down, showCategory = 10, color = "pvalue")
  print(p2)
  dev.off()
}

if(T){
  ## 3.6 up genes ----
  up_genes <- subset(Genesymbol, logFC > cutoff)
  
  ### 3.6.1 GO-up ----
  # gene ID转换 
  gene <- bitr(up_genes$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)
  
  # GO富集分析
  kkd <- enrichGO(gene = gene$ENTREZID, # 导入基因的ENTREZID编号
                  OrgDb = GO_database, # 用到的数据库（人类是：org.Hs.eg.db）
                  keyType = "ENTREZID", # 设定读取的gene ID类型
                  ont = "ALL", 
                  pvalueCutoff = 1,
                  pAdjustMethod = "BH",
                  qvalueCutoff = 1, # 设定q值阈值
                  readable = T
  )
  ### 3.6.2 KEGG-up ----
  # KEGG富集分析
  kk <- enrichKEGG(gene = gene$ENTREZID,
                   keyType = "kegg",
                   organism = KEGG_database,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 1,
                   use_internal_data = T
  )
  # 设置geneID可读格式
  kk <- setReadable(kk, OrgDb='org.Hs.eg.db', keyType='ENTREZID')
  
  # GO、KEGG结果整合
  result <- list(enrichGO = kkd, enrichKEGG = kk)
  # 结果标记为上调
  result_up <- result
  kkd_up <- result_up$enrichGO
  kk_up <- result_up$enrichKEGG
  
  ### 3.6.3 上调res_output ----
  
  # 导出上调enrichGO
  write.xlsx(kkd_up@result, file = paste0(dir_gokegg, "/GO_up.xlsx"), rowNames = FALSE)
  
  # dotplot
  pdf(file = paste0(dir_gokegg, "/GO_up.pdf"), width = 5.5, height = 7)
  p3 <- dotplot(kkd_up, showCategory = 5, color = "pvalue", split = "ONTOLOGY") + 
    facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free')
  print(p3)
  dev.off()
  
  # 导出上调enrichKEGG
  write.xlsx(kk_up@result, file = paste0(dir_gokegg, "/KEGG_up.xlsx"), rowNames = FALSE)
  
  # dotplot
  pdf(file = paste0(dir_gokegg, "/KEGG_up.pdf"), width = 6, height = 5)
  p4 <- dotplot(kk_up,showCategory = 10, color = "pvalue")
  print(p4)
  dev.off()
}

if(T){
  ## 3.7 All DE genes ----
  All_genes <- subset(Genesymbol, abs(logFC) > cutoff)
  
  ### 3.7.1 GO ----
  # gene ID转换 
  gene <- clusterProfiler::bitr(All_genes$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)
  
  # GO富集分析
  go <- clusterProfiler::enrichGO(gene = gene$ENTREZID, 
                                  OrgDb = GO_database, 
                                  keyType = "ENTREZID", 
                                  ont = "ALL", 
                                  pvalueCutoff = 1,
                                  qvalueCutoff = 1, 
                                  readable = T)
  ### 3.7.2 KEGG ----
  # KEGG富集分析
  kegg <- clusterProfiler::enrichKEGG(gene = gene$ENTREZID,
                                      keyType = "kegg",
                                      organism = KEGG_database,
                                      pAdjustMethod = "BH",
                                      pvalueCutoff = 0.05,
                                      qvalueCutoff = 1,
                                      use_internal_data = T
  )
  # 设置geneID可读格式
  kk <- setReadable(kk, OrgDb='org.Hs.eg.db', keyType='ENTREZID')
  
  # GO、KEGG结果整合
  result <- list(enrichGO = go, enrichKEGG = kegg)
  # 结果标记为上调
  result_up <- result
  GO <- result_up$enrichGO
  KEGG <- result_up$enrichKEGG
  
  ### 3.7.3 res_output ----
  
  # 导出上调enrichGO
  write.xlsx(GO@result, file = paste0(dir_gokegg, "/GO_all_DE.xlsx"), rowNames = FALSE)
  
  # dotplot
  pdf(file = paste0(dir_gokegg, "/GO_all_DE.pdf"), width = 5.5, height = 7) # 如果用重叠 width = 6.5, height = 8.5 
  p5 <- dotplot(GO, showCategory = 5, color = "pvalue", split = "ONTOLOGY") + 
    facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free') 
  # theme(axis.text.y = element_text(angle = 0, hjust = 1)) +
  # scale_y_discrete(labels = function(x) str_wrap(x, width = 35))  # 控制每行最多显示40个字符
  print(p5)
  dev.off()
  
  # 导出上调enrichKEGG
  write.xlsx(KEGG@result, file = paste0(dir_gokegg, "/KEGG_all_DE.xlsx"), rowNames = FALSE)
  
  # dotplot
  pdf(file = paste0(dir_gokegg, "/KEGG_all_DE.pdf"), width = 6, height = 5)
  p6 <- dotplot(KEGG,showCategory = 10, color = "pvalue")
  print(p6)
  dev.off()
}

## 3.8 查看上下调的通路的数量 ----
# 下调GO 
table(kkd_down@result$p.adjust<0.05)
# 下调KEGG 
table(kk_down@result$p.adjust<0.05)
# 上调GO 
table(kkd_up@result$p.adjust<0.05)
# 上调KEGG 
table(kk_up@result$p.adjust<0.05)
