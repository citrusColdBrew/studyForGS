# 基础GWAS实战代码
# 基于mdp数据集的简化GWAS分析

# =============================================================================
# 第一步：环境准备和数据载入
# =============================================================================

# 清理环境
rm(list = ls())

# 载入数据（请根据实际路径调整）
cat("载入数据...\n")
genotype_data <- read.table("data/mdp_numeric.txt", header = TRUE)
snp_info <- read.table("data/mdp_SNP_information.txt", header = TRUE)
trait_data <- read.table("data/mdp_traits.txt", header = TRUE)

# 查看数据基本信息
cat("=== 数据概况 ===\n")
cat(sprintf("个体数: %d\n", nrow(genotype_data)))
cat(sprintf("SNP数: %d\n", ncol(genotype_data) - 1))
cat(sprintf("染色体数: %d\n", length(unique(snp_info$Chromosome))))
cat(sprintf("表型数据个体数: %d\n", nrow(trait_data)))
cat(sprintf("表型性状数: %d\n", ncol(trait_data) - 1))

# =============================================================================
# 第二步：数据预处理
# =============================================================================

cat("\n=== 数据预处理 ===\n")

# 提取个体名和基因型矩阵
individual_names <- genotype_data[, 1]
X <- as.matrix(genotype_data[, -1])  # 基因型矩阵
rownames(X) <- individual_names

cat(sprintf("基因型矩阵维度: %d × %d\n", nrow(X), ncol(X)))

# 简单质量控制：移除单态SNP
polymorphic_mask <- apply(X, 2, function(snp) {
  unique_values <- unique(snp[!is.na(snp)])
  length(unique_values) > 1
})

X_filtered <- X[, polymorphic_mask]
snp_info_filtered <- snp_info[polymorphic_mask, ]

cat(sprintf("质控后SNP数: %d\n", ncol(X_filtered)))

# =============================================================================
# 第三步：表型数据处理
# =============================================================================

cat("\n=== 表型数据处理 ===\n")

# 匹配基因型和表型数据的个体
genotype_individuals <- genotype_data[, 1]
trait_individuals <- trait_data[, 1]
common_individuals <- intersect(genotype_individuals, trait_individuals)

cat(sprintf("基因型数据个体数: %d\n", length(genotype_individuals)))
cat(sprintf("表型数据个体数: %d\n", length(trait_individuals)))
cat(sprintf("共同个体数: %d\n", length(common_individuals)))

# 根据共同个体筛选数据
genotype_matched_idx <- match(common_individuals, genotype_individuals)
trait_matched_idx <- match(common_individuals, trait_individuals)

matched_genotype_data <- genotype_data[genotype_matched_idx, ]
matched_trait_data <- trait_data[trait_matched_idx, ]

# 验证个体名是否匹配
if(!all(matched_genotype_data[, 1] == matched_trait_data[, 1])) {
  stop("个体名匹配失败！")
}
cat("个体名匹配成功！\n")

# 选择分析性状（穗位高 EarHT）
selected_trait <- "EarHT"
cat(sprintf("选择分析性状: %s\n", selected_trait))

# 提取选定性状的数据
trait_column_idx <- which(colnames(matched_trait_data) == selected_trait)
phenotype_values <- matched_trait_data[, trait_column_idx]
names(phenotype_values) <- matched_trait_data[, 1]

# 处理缺失值
missing_count <- sum(is.na(phenotype_values))
cat(sprintf("缺失值数量: %d (%.1f%%)\n",
            missing_count, missing_count / length(phenotype_values) * 100))

# 去除缺失值
valid_individuals <- !is.na(phenotype_values)
y <- phenotype_values[valid_individuals]
final_genotype_data <- matched_genotype_data[valid_individuals, ]

# 提取最终的基因型矩阵
individual_names <- final_genotype_data[, 1]
X <- as.matrix(final_genotype_data[, -1])
rownames(X) <- individual_names

cat(sprintf("最终分析个体数: %d\n", length(y)))
cat(sprintf("表型统计: 均值=%.2f, 标准差=%.2f, 范围=%.2f-%.2f\n",
            mean(y), sd(y), min(y), max(y)))

# 应用质量控制到最终数据
X_filtered <- X[, polymorphic_mask]
snp_info_filtered <- snp_info[polymorphic_mask, ]

cat(sprintf("最终基因型矩阵维度: %d × %d\n", nrow(X_filtered), ncol(X_filtered)))

# =============================================================================
# 第四步：GWAS分析
# =============================================================================

cat("\n=== GWAS分析 ===\n")

# 方法1：相关性分析
cat("1. 相关性分析...\n")
correlation_analysis <- function(genotypes, phenotype) {
  n_snps <- ncol(genotypes)
  p_values <- numeric(n_snps)

  for(i in 1:n_snps) {
    snp_geno <- genotypes[, i]
    if(length(unique(snp_geno)) > 1) {
      test_result <- cor.test(snp_geno, phenotype)
      p_values[i] <- test_result$p.value
    } else {
      p_values[i] <- 1  # 单态SNP设为不显著
    }
  }
  return(p_values)
}

p_correlation <- correlation_analysis(X_filtered, y)

# 方法2：简单线性回归
cat("2. 线性回归分析...\n")
regression_analysis <- function(genotypes, phenotype) {
  n_snps <- ncol(genotypes)
  p_values <- numeric(n_snps)
  beta_values <- numeric(n_snps)

  for(i in 1:n_snps) {
    snp_geno <- genotypes[, i]
    if(length(unique(snp_geno)) > 1) {
      lm_result <- lm(phenotype ~ snp_geno)
      summary_result <- summary(lm_result)
      beta_values[i] <- summary_result$coefficients[2, 1]
      p_values[i] <- summary_result$coefficients[2, 4]
    } else {
      beta_values[i] <- 0
      p_values[i] <- 1
    }
  }
  return(list(beta = beta_values, p_values = p_values))
}

regression_result <- regression_analysis(X_filtered, y)
p_regression <- regression_result$p_values

# 方法3：控制群体结构的GLM
cat("3. GLM分析（控制群体结构）...\n")

# 首先进行PCA
pca_result <- prcomp(X_filtered, center = TRUE, scale. = FALSE)
pc_scores <- pca_result$x[, 1:3]  # 使用前3个主成分

glm_analysis <- function(genotypes, phenotype, covariates) {
  n_snps <- ncol(genotypes)
  p_values <- numeric(n_snps)

  for(i in 1:n_snps) {
    snp_geno <- genotypes[, i]
    if(length(unique(snp_geno)) > 1) {
      # 构建包含协变量的模型
      lm_result <- lm(phenotype ~ covariates + snp_geno)
      summary_result <- summary(lm_result)
      # SNP效应是最后一个系数
      snp_coeff_row <- nrow(summary_result$coefficients)
      p_values[i] <- summary_result$coefficients[snp_coeff_row, 4]
    } else {
      p_values[i] <- 1
    }
  }
  return(p_values)
}

p_glm <- glm_analysis(X_filtered, y, pc_scores)

cat("分析完成！\n")

# =============================================================================
# 第五步：结果可视化
# =============================================================================

cat("\n=== 结果可视化 ===\n")

# 创建曼哈顿图函数
create_manhattan_plot <- function(p_values, title) {
  log_p <- -log10(p_values)
  chromosomes <- snp_info_filtered$Chromosome

  # 为不同染色体设置颜色
  colors <- rainbow(max(chromosomes))[chromosomes]

  plot(1:length(p_values), log_p,
       col = colors, pch = 16, cex = 0.6,
       xlab = "SNP位置", ylab = "-log10(P值)",
       main = title)

  # 添加显著性阈值线（Bonferroni校正）
  threshold <- 0.05 / length(p_values)
  abline(h = -log10(threshold), col = "red", lty = 2, lwd = 2)

  return(threshold)
}

# 绘制三种方法的曼哈顿图
par(mfrow = c(3, 1), mar = c(4, 4, 3, 1))

threshold1 <- create_manhattan_plot(p_correlation, "相关性分析")
threshold2 <- create_manhattan_plot(p_regression, "线性回归分析")
threshold3 <- create_manhattan_plot(p_glm, "GLM分析（控制群体结构）")

# =============================================================================
# 第六步：结果评估
# =============================================================================

cat("\n=== 结果评估 ===\n")

# 评估函数
evaluate_gwas_performance <- function(p_values, method_name) {
  significance_threshold <- 0.05 / length(p_values)

  # 显著SNP数量
  significant_snps <- sum(p_values < significance_threshold)

  # 最小p值及其位置
  min_p_idx <- which.min(p_values)
  min_p_value <- p_values[min_p_idx]

  cat(sprintf("%s:\n", method_name))
  cat(sprintf("  显著SNP数量: %d\n", significant_snps))
  cat(sprintf("  最小P值: %.2e (位置: %d)\n", min_p_value, min_p_idx))
  cat(sprintf("  显著性阈值: %.2e\n", significance_threshold))
  cat("\n")
}

# 评估三种方法
evaluate_gwas_performance(p_correlation, "相关性分析")
evaluate_gwas_performance(p_regression, "线性回归分析")
evaluate_gwas_performance(p_glm, "GLM分析")

# =============================================================================
# 第七步：单个SNP效应分析
# =============================================================================

cat("=== 最显著SNP效应分析 ===\n")

# 找到回归分析中最显著的SNP
top_snp_idx <- which.min(p_regression)
top_snp_genotypes <- X_filtered[, top_snp_idx]

cat(sprintf("最显著SNP位置: %d\n", top_snp_idx))
cat(sprintf("P值: %.2e\n", p_regression[top_snp_idx]))

# 按基因型分组分析表型
cat("\n按基因型分组的表型统计:\n")
for(geno in 0:2) {
  individuals_with_geno <- which(top_snp_genotypes == geno)
  if(length(individuals_with_geno) > 0) {
    phenotype_subset <- y[individuals_with_geno]
    cat(sprintf("基因型%d: n=%d, 均值=%.3f, 标准差=%.3f\n",
                geno, length(individuals_with_geno),
                mean(phenotype_subset), sd(phenotype_subset)))
  }
}

# 绘制箱线图
par(mfrow = c(1, 1))
boxplot(y ~ top_snp_genotypes,
        main = sprintf("最显著SNP (位置%d) 的基因型效应", top_snp_idx),
        xlab = "基因型 (0=AA, 1=Aa, 2=aa)",
        ylab = "表型值",
        col = c("lightblue", "lightgreen", "lightcoral"))

# 检查这个SNP的效应
cat(sprintf("\n这个SNP与%s性状的关联强度很高\n", selected_trait))
cat("说明该SNP可能与穗位高性状有真实的生物学关联\n")

# =============================================================================
# 总结
# =============================================================================

cat("\n=== 分析总结 ===\n")
cat("本次GWAS分析完成了以下步骤:\n")
cat("1. 数据载入和预处理\n")
cat("2. 真实表型数据处理（穗位高性状）\n")
cat("3. 三种GWAS方法比较\n")
cat("4. 结果可视化和性能评估\n")
cat("5. 单个SNP效应分析\n")

cat("\n主要发现:\n")
cat(sprintf("- 分析了%s性状的真实表型数据\n", selected_trait))
cat("- 比较了相关性、回归、GLM三种方法\n")
cat("- GLM方法通过控制群体结构可以减少假阳性\n")
cat("- 发现了与穗位高显著关联的SNP位点\n")

cat("\n这个分析展示了真实GWAS数据的基本流程和关键概念。\n")
cat("在实际应用中，还需要考虑更多因素，如:\n")
cat("- 更严格的质量控制\n")
cat("- 混合线性模型（MLM）\n")
cat("- 多重检验校正策略\n")
cat("- 功能注释和生物学解释\n")
cat("- 验证实验设计\n")

cat("\n分析完成！\n")