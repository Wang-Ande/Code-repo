
#' @title check_sva_overfitting
#' @description
#' 检查 SVA 校正是否存在过拟合或信号损失的风险。
#' 包括：差异表达前后 DEG 数量比较、SV 与分组变量的相关性、PCA 主成分解释度分析、PCA 和火山图可视化。
#'
#' @param expr_mat 表达矩阵 (log2(CPM+1))，行为基因，列为样本。
#' @param group 因子向量，表示样本所属组，levels=（，） 默认第一个为对照组。
#' @param svobj sva() 返回的对象，包含 surrogate variables。
#' @param mod 设计矩阵，通常由 group 构建（如 model.matrix(~ group)）。
#' @param coef_name 指定差异分析中感兴趣的变量名称（默认取 mod 的第二列）。
#' @param lfc_threshold logFC 阈值（默认 1）。
#' @param fdr_threshold FDR 阈值（默认 0.05）。
#'
#' @return 控制台输出诊断信息并绘图：PCA前后图 + 火山图前后对比。
#' @export

check_sva_overfitting <- function(expr_mat, group, svobj, mod, coef_name = NULL, lfc_threshold = 1, fdr_threshold = 0.05) {
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(reshape2)

  message(">>> [1] 差异表达基因数量比较")

  fit0 <- lmFit(expr_mat, mod)
  fit0 <- eBayes(fit0, trend = TRUE)
  if (is.null(coef_name)) {
    coef_name <- colnames(mod)[2]
  }
  res0 <- topTable(fit0, coef = coef_name, number = Inf, adjust = "BH")
  deg0 <- sum(res0$adj.P.Val < fdr_threshold & abs(res0$logFC) > lfc_threshold)

  mod_sva <- cbind(mod, svobj$sv)
  fit1 <- lmFit(expr_mat, mod_sva)
  fit1 <- eBayes(fit1, trend = TRUE)
  res1 <- topTable(fit1, coef = coef_name, number = Inf, adjust = "BH")
  deg1 <- sum(res1$adj.P.Val < fdr_threshold & abs(res1$logFC) > lfc_threshold)

  message(sprintf("    → 校正前 DEGs: %d 个", deg0))
  message(sprintf("    → 校正后 DEGs: %d 个", deg1))

  if (deg1 < 0.5 * deg0) {
    message("⚠️  差异表达基因数量下降过多，可能存在过矫正或过拟合。")
  } else {
    message("✅  差异表达基因数量变化合理。")
  }

  message("\n>>> [2] SV 与分组变量的相关性")
  cor_vec <- apply(svobj$sv, 2, function(x) cor(x, as.numeric(group)))
  print(round(cor_vec, 3))

  if (any(abs(cor_vec) > 0.5)) {
    message("⚠️  存在与分组强相关的 surrogate variable（|cor| > 0.5），可能吃掉生物信号。")
  } else {
    message("✅  surrogate variables 与 group 变量相关性低，符合预期。")
  }

  message("\n>>> [3] PCA 主成分解释度分析")
  pca_raw <- prcomp(t(expr_mat), scale. = TRUE)
  expr_corrected <- removeBatchEffect(expr_mat, covariates = svobj$sv, design = mod)
  pca_corr <- prcomp(t(expr_corrected), scale. = TRUE)

  var_raw <- (pca_raw$sdev^2) / sum(pca_raw$sdev^2)
  var_corr <- (pca_corr$sdev^2) / sum(pca_corr$sdev^2)

  pc1_before <- round(var_raw[1] * 100, 2)
  pc1_after  <- round(var_corr[1] * 100, 2)
  message(sprintf("    → PC1 variance before correction: %.2f%%", pc1_before))
  message(sprintf("    → PC1 variance after correction:  %.2f%%", pc1_after))

  if (pc1_after < 0.1 * pc1_before) {
    message("⚠️  PCA 主成分变异急剧降低，可能存在过拟合或数据扁平化。")
  } else {
    message("✅  PCA 主成分解释度变化正常。")
  }

  message("\n>>> [4] 绘制PCA和火山图")

  # PCA plot
  # 计算方差解释度
  pc_var_raw <- (pca_raw$sdev^2 / sum(pca_raw$sdev^2))[1:2] * 100
  pc_var_corr <- (pca_corr$sdev^2 / sum(pca_corr$sdev^2))[1:2] * 100
  
  # 构建 PCA dataframe
  pca_df_raw <- data.frame(PC1 = pca_raw$x[,1], PC2 = pca_raw$x[,2], group = group)
  pca_df_corr <- data.frame(PC1 = pca_corr$x[,1], PC2 = pca_corr$x[,2], group = group)
  
  # 添加主成分百分比到标签
  xlab_raw <- sprintf("PC1 (%.1f%%)", pc_var_raw[1])
  ylab_raw <- sprintf("PC2 (%.1f%%)", pc_var_raw[2])
  
  xlab_corr <- sprintf("PC1 (%.1f%%)", pc_var_corr[1])
  ylab_corr <- sprintf("PC2 (%.1f%%)", pc_var_corr[2])
  
  # 绘图
  p1 <- ggplot(pca_df_raw, aes(PC1, PC2, color = group)) +
    geom_point(size = 3) +
    labs(title = "PCA Before Correction", x = xlab_raw, y = ylab_raw) +
    scale_color_discrete(name = "Group") +
    theme_minimal()
  
  p2 <- ggplot(pca_df_corr, aes(PC1, PC2, color = group)) +
    geom_point(size = 3) +
    labs(title = "PCA After Correction", x = xlab_corr, y = ylab_corr) +
    scale_color_discrete(name = "Group") +
    theme_minimal()
  
  print(p1); print(p2)

  # Volcano plots
  res0$threshold <- as.factor(res0$adj.P.Val < fdr_threshold & abs(res0$logFC) > lfc_threshold)
  res1$threshold <- as.factor(res1$adj.P.Val < fdr_threshold & abs(res1$logFC) > lfc_threshold)

  res0$gene <- rownames(res0)
  res1$gene <- rownames(res1)

  v1 <- ggplot(res0, aes(logFC, -log10(adj.P.Val), color = threshold)) +
    geom_point(alpha = 0.6) + theme_minimal() +
    labs(title = "Volcano Plot (Before Correction)") +
    scale_color_manual(values = c("grey", "red"), name = "Significant")

  v2 <- ggplot(res1, aes(logFC, -log10(adj.P.Val), color = threshold)) +
    geom_point(alpha = 0.6) + theme_minimal() +
    labs(title = "Volcano Plot (After Correction)") +
    scale_color_manual(values = c("grey", "red"), name = "Significant")

  print(v1); print(v2)

  message("\n✅ 检查与图形输出完成。")
}
