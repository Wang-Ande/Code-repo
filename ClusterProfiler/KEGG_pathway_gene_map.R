# library(pathview)
# library(clusterprofiler)
kegg <- read.csv("./03_Result/GO&KEGG/OCI_AML2/Low_vs_Con/KEGG_down.csv")
gene_ids <- "9538/8795/10572/355/3732/1026/1029/581/50484/64326/1643/51512/27244"
gene_list <- unlist(strsplit(gene_ids, "/"))
gene <- clusterProfiler::bitr(gene_list, fromType = 'ENTREZID', toType = 'SYMBOL', OrgDb = org.Hs.eg.db)


# 输入：一个命名向量，名称是Entrez ID，数值是logFC或任意表达量
enriched_genes  <- gene$ENTREZID
gene_data <- rep(-1, length(enriched_genes))
names(gene_data) <- enriched_genes

# 绘制并高亮基因
pathview(gene.data = gene_data, 
         pathway.id = "04115",       # 通路ID
         species = "hsa",            # 人类
         gene.idtype = "entrez",     # ID类型
         limit = list(gene=1, cpd=1),
         bins = list(gene=10, cpd=10),
         low = list(gene="#66F1A9", cpd="#08519C"),     # 深蓝（下调）
         mid = list(gene="#F7F7F7", cpd="#F0F0F0"),      # 白灰色（中性）
         high = list(gene="#BD0026", cpd="#FC4E2A"))     # 深红（上调）


