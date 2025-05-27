
#' @title test_sv_range
#' @description
#' 通过逐步调整 n.sv 的数量，评估批次效应矫正对表达矩阵主成分结构（特别是 PC1）的影响，
#' 以辅助判断最优 surrogate variable 个数。
#'
#' @param expr_mat 表达矩阵（基因 × 样本），建议为 log2(CPM+1) 数据。
#' @param mod 差异分析设计矩阵（如 model.matrix(~ group)）。
#' @param group 因子型分组变量（如 sensitive/resistant）。
#' @param sv_range 要评估的 surrogate variable 数量范围，例如 0:20。
#'
#' @return data.frame，包含每个 n.sv 对应的 PC1 方差解释百分比。
#' @export

test_sv_range <- function(expr_mat, mod, group, sv_range) {
  library(limma)
  library(sva)
  library(ggplot2)

  expr_mat <- as.matrix(expr_mat)
  mod <- as.matrix(mod)
  pc1_list <- numeric(length(sv_range))

  for (i in seq_along(sv_range)) {
    n <- sv_range[i]

    if (n == 0) {
      expr_corr <- expr_mat  # 不做矫正
    } else {
      mod0 <- model.matrix(~1, data = data.frame(group = group))
      svobj <- sva(expr_mat, mod, mod0, n.sv = n)
      expr_corr <- removeBatchEffect(expr_mat, covariates = svobj$sv, design = mod)
    }

    pca <- prcomp(t(expr_corr), scale. = TRUE)
    var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
    pc1_list[i] <- var_explained[1] * 100
  }

  data.frame(n_sv = sv_range, PC1_var = pc1_list)
  
  # 绘制趋势图
  p1 <- ggplot(result, aes(x = n_sv, y = PC1_var)) +
       geom_point() + geom_line() +
       labs(title = "PC1 Explained Variance vs n.sv",
            x = "Number of Surrogate Variables (n.sv)",
            y = "PC1 Variance Explained (%)") +
       theme_minimal()
  print(p1)
}
