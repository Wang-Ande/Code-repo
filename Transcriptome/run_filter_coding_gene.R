filter_coding_genes <- function(expr_matrix, coding_gene_file = "../../Code_Box/Transcriptome/coding_gene_ids.rds") {
  # expr_matrix: 行名为 gene_id，列为样本
  # coding_gene_file: rds 文件路径
  # 使用示例
  # clean_matrix <- filter_coding_genes(expr_matrix)
  coding_gene_ids <- readRDS(coding_gene_file)
  
  # 提取 ENSEMBL ID（点前部分）和 symbol（管道符号第二部分）
  ense <- sapply(rownames(expr_matrix), function(x) strsplit(as.character(x), "[.]")[[1]][1])
  symbol <- sapply(rownames(expr_matrix), function(x) strsplit(as.character(x), "[|]")[[1]][2])
  
  # 判断哪些是编码基因
  keep_idx <- ense %in% coding_gene_ids
  
  # 过滤矩阵
  expr_matrix_clean <- expr_matrix[keep_idx, , drop = FALSE]
  
  # 添加两列信息
  expr_matrix_clean <- cbind(
    ense = ense[keep_idx],
    symbol = symbol[keep_idx],
    expr_matrix_clean
  )
  
  return(expr_matrix_clean)
}

# 使用示例
# clean_matrix <- filter_coding_genes_with_symbol(expr_matrix)
