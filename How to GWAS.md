# How to GWAS：从零开始的全基因组关联分析教程

## 目录
1. [GWAS是什么？](#1-gwas是什么)
2. [数据准备和理解](#2-数据准备和理解)
3. [第一步：数据载入和检查](#3-第一步数据载入和检查)
4. [第二步：数据预处理](#4-第二步数据预处理)
5. [第三步：表型数据获取或模拟](#5-第三步表型数据获取或模拟)
6. [第四步：基本关联分析](#6-第四步基本关联分析)
7. [第五步：结果可视化](#7-第五步结果可视化)
8. [第六步：统计质量控制](#8-第六步统计质量控制)
9. [第七步：改进的分析方法](#9-第七步改进的分析方法)
10. [完整代码示例](#10-完整代码示例)

---

## 1. GWAS是什么？

### 生物学背景
GWAS（Genome-Wide Association Study，全基因组关联分析）是寻找基因组中与性状相关的DNA变异的方法。

**核心问题**：哪些基因变异影响了我们关心的性状（如身高、疾病易感性、农作物产量）？

**基本原理**：
- 如果某个DNA位点真的影响性状，那么不同基因型的个体在表型上应该有显著差异
- 通过统计方法检测这种关联关系

### 数据需求
- **基因型数据**：每个个体在每个DNA位点的变异情况
- **表型数据**：每个个体的性状测量值
- **可选的协变量**：年龄、性别、群体信息等

---

## 2. 数据准备和理解

### 我们的数据集
1. **mdp_numeric.txt**：玉米多样性面板的基因型数据
2. **mdp_SNP_information.txt**：SNP位置信息
3. **mdp_traits.txt**：玉米农艺性状的表型数据

### 数据格式说明

**mdp_numeric.txt结构**：
```
taxa        SNP1  SNP2  SNP3  ...
Individual1   0     1     2   ...
Individual2   1     2     0   ...
...
```
- 第一列：个体名称
- 其余列：SNP基因型，编码为0/1/2
  - 0 = AA（参考纯合子）
  - 1 = Aa（杂合子）
  - 2 = aa（变异纯合子）

**mdp_SNP_information.txt结构**：
```
SNP          Chromosome  Position
SNP1         1           157104
**mdp_traits.txt结构**：
```
Taxa        EarHT   dpoll   EarDia
Individual1 59.5    NaN     NaN
Individual2 65.5    59.5    32.21933
Individual3 81.13   71.5    32.421
...
```
- 第一列：个体名称（与基因型数据匹配）
- EarHT：穗位高（cm）
- dpoll：开花期（天数）
- EarDia：穗径（mm）
- NaN：缺失值

---

## 3. 第一步：数据载入和检查

### 3.1 设置工作环境

```r
# 清理环境
rm(list = ls())

# 设置工作目录（根据你的实际路径调整）
setwd("/Users/your_path/statistical genomics")

# 检查数据文件是否存在
if(!file.exists("data/mdp_numeric.txt")) {
  stop("基因型数据文件不存在！请检查路径")
}
if(!file.exists("data/mdp_SNP_information.txt")) {
  stop("SNP信息文件不存在！请检查路径")
}
if(!file.exists("data/mdp_traits.txt")) {
  stop("表型数据文件不存在！请检查路径")
}
```

### 3.2 载入数据

```r
# 载入基因型数据
cat("正在载入基因型数据...\n")
genotype_data <- read.table("data/mdp_numeric.txt", header = TRUE, sep = "\t")

# 载入SNP信息
cat("正在载入SNP位置信息...\n")
snp_info <- read.table("data/mdp_SNP_information.txt", header = TRUE, sep = "\t")

# 载入表型数据
cat("正在载入表型数据...\n")
trait_data <- read.table("data/mdp_traits.txt", header = TRUE, sep = "\t")

cat("数据载入完成！\n")
```

### 3.3 数据基本信息检查

```r
# 检查数据维度
cat("=== 数据基本信息 ===\n")
cat(sprintf("个体数量: %d\n", nrow(genotype_data)))
cat(sprintf("SNP数量: %d\n", ncol(genotype_data) - 1))  # 减1是因为第一列是个体名
cat(sprintf("SNP信息条目: %d\n", nrow(snp_info)))
cat(sprintf("表型数据个体数: %d\n", nrow(trait_data)))
cat(sprintf("表型性状数: %d\n", ncol(trait_data) - 1))  # 减1是因为第一列是个体名

# 检查前几行数据
cat("\n=== 基因型数据前5行前5列 ===\n")
print(genotype_data[1:5, 1:6])

cat("\n=== SNP信息前5行 ===\n")
print(snp_info[1:5, ])

cat("\n=== 表型数据前5行 ===\n")
print(trait_data[1:5, ])

# 检查表型数据的基本统计
cat("\n=== 表型数据基本统计 ===\n")
for(i in 2:ncol(trait_data)) {
  trait_name <- colnames(trait_data)[i]
  trait_values <- trait_data[, i]
  # 去除缺失值计算统计量
  valid_values <- trait_values[!is.na(trait_values)]
  if(length(valid_values) > 0) {
    cat(sprintf("%s: n=%d, 均值=%.2f, 标准差=%.2f, 范围=%.2f-%.2f\n",
                trait_name, length(valid_values),
                mean(valid_values), sd(valid_values),
                min(valid_values), max(valid_values)))
  }
}

# 检查染色体分布
cat("\n=== 染色体分布 ===\n")
chr_counts <- table(snp_info$Chromosome)
print(chr_counts)
```

---

## 4. 第二步：数据预处理

### 4.1 提取和整理基因型矩阵

```r
# 提取个体名称
individual_names <- genotype_data[, 1]
cat(sprintf("个体名称示例: %s\n", paste(head(individual_names, 3), collapse = ", ")))

# 提取基因型矩阵（去除第一列个体名）
genotype_matrix <- as.matrix(genotype_data[, -1])

# 设置行名和列名
rownames(genotype_matrix) <- individual_names
colnames(genotype_matrix) <- colnames(genotype_data)[-1]

cat(sprintf("基因型矩阵维度: %d个体 × %d个SNP\n",
            nrow(genotype_matrix), ncol(genotype_matrix)))
```

### 4.2 数据质量检查

```r
# 检查基因型编码是否正确（应该只有0, 1, 2）
unique_values <- unique(as.vector(genotype_matrix))
cat(sprintf("基因型编码: %s\n", paste(sort(unique_values), collapse = ", ")))

# 检查缺失值
missing_count <- sum(is.na(genotype_matrix))
cat(sprintf("缺失值数量: %d (%.2f%%)\n",
            missing_count, missing_count / length(genotype_matrix) * 100))

# 计算每个SNP的等位基因频率
calculate_allele_freq <- function(genotypes) {
  # 计算变异等位基因频率
  freq <- mean(genotypes, na.rm = TRUE) / 2
  return(freq)
}

allele_freqs <- apply(genotype_matrix, 2, calculate_allele_freq)
cat(sprintf("等位基因频率范围: %.3f - %.3f\n",
            min(allele_freqs, na.rm = TRUE), max(allele_freqs, na.rm = TRUE)))
```

### 4.3 质量控制过滤

```r
# 过滤单态SNP（所有个体基因型相同的SNP）
polymorphic_snps <- apply(genotype_matrix, 2, function(x) {
  length(unique(x[!is.na(x)])) > 1
})

cat(sprintf("单态SNP数量: %d\n", sum(!polymorphic_snps)))
cat(sprintf("保留的多态SNP数量: %d\n", sum(polymorphic_snps)))

# 过滤极低频率SNP（MAF < 0.05）
maf_threshold <- 0.05
maf <- pmin(allele_freqs, 1 - allele_freqs)  # 取较小等位基因频率
high_maf_snps <- maf >= maf_threshold

cat(sprintf("低频SNP数量 (MAF < %.2f): %d\n",
            maf_threshold, sum(!high_maf_snps)))

# 应用过滤条件
good_snps <- polymorphic_snps & high_maf_snps
filtered_genotype_matrix <- genotype_matrix[, good_snps]
filtered_snp_info <- snp_info[good_snps, ]

cat(sprintf("质控后SNP数量: %d\n", ncol(filtered_genotype_matrix)))
```

---

## 5. 第三步：表型数据处理

### 5.1 数据匹配和整理

```r
# 匹配基因型和表型数据的个体
# 基因型数据的个体名
genotype_individuals <- genotype_data[, 1]
# 表型数据的个体名
trait_individuals <- trait_data[, 1]

# 找到两个数据集的共同个体
common_individuals <- intersect(genotype_individuals, trait_individuals)
cat(sprintf("基因型数据个体数: %d\n", length(genotype_individuals)))
cat(sprintf("表型数据个体数: %d\n", length(trait_individuals)))
cat(sprintf("共同个体数: %d\n", length(common_individuals)))

# 根据共同个体筛选数据
genotype_matched_idx <- match(common_individuals, genotype_individuals)
trait_matched_idx <- match(common_individuals, trait_individuals)

# 获取匹配的数据
matched_genotype_data <- genotype_data[genotype_matched_idx, ]
matched_trait_data <- trait_data[trait_matched_idx, ]

# 验证个体名是否匹配
if(!all(matched_genotype_data[, 1] == matched_trait_data[, 1])) {
  stop("个体名匹配失败！")
}
cat("个体名匹配成功！\n")
```

### 5.2 选择分析性状

```r
# 查看可用的性状
trait_names <- colnames(matched_trait_data)[-1]  # 去除第一列个体名
cat("可用性状:\n")
for(i in 1:length(trait_names)) {
  trait_values <- matched_trait_data[, i + 1]
  valid_count <- sum(!is.na(trait_values))
  cat(sprintf("%d. %s (有效数据: %d/%d)\n",
              i, trait_names[i], valid_count, length(trait_values)))
}

# 选择一个性状进行分析（这里选择穗位高 EarHT）
selected_trait <- "EarHT"
cat(sprintf("\n选择分析性状: %s\n", selected_trait))

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
final_phenotypes <- phenotype_values[valid_individuals]
final_genotype_data <- matched_genotype_data[valid_individuals, ]

cat(sprintf("最终分析个体数: %d\n", length(final_phenotypes)))
cat(sprintf("表型统计: 均值=%.2f, 标准差=%.2f, 范围=%.2f-%.2f\n",
            mean(final_phenotypes), sd(final_phenotypes),
            min(final_phenotypes), max(final_phenotypes)))
```

### 5.3 表型数据可视化

```r
# 绘制表型分布直方图
par(mfrow = c(1, 2))

# 直方图
hist(final_phenotypes, breaks = 20, col = "lightblue",
     main = paste("表型分布:", selected_trait),
     xlab = "表型值", ylab = "频数")

# 箱线图
boxplot(final_phenotypes, col = "lightgreen",
        main = paste("表型分布:", selected_trait),
        ylab = "表型值")

# 检查数据是否符合正态分布
shapiro_test <- shapiro.test(final_phenotypes)
cat(sprintf("正态性检验 (Shapiro-Wilk): W=%.4f, p=%.4f\n",
            shapiro_test$statistic, shapiro_test$p.value))

if(shapiro_test$p.value > 0.05) {
  cat("表型数据近似正态分布，适合线性模型分析\n")
} else {
  cat("表型数据偏离正态分布，可能需要数据转换\n")
}
```

---

## 6. 第四步：基本关联分析

### 6.1 相关性分析法

```r
cat("=== 开始GWAS分析 ===\n")

# 方法1：简单相关性分析
perform_correlation_gwas <- function(genotypes, phenotype) {
  n_snps <- ncol(genotypes)
  n_individuals <- nrow(genotypes)

  # 存储结果
  correlations <- numeric(n_snps)
  p_values <- numeric(n_snps)

  for(i in 1:n_snps) {
    snp_genotypes <- genotypes[, i]

    # 计算相关性
    if(length(unique(snp_genotypes)) > 1) {  # 确保有变异
      cor_result <- cor.test(snp_genotypes, phenotype)
      correlations[i] <- cor_result$estimate
      p_values[i] <- cor_result$p.value
    } else {
      correlations[i] <- 0
      p_values[i] <- 1
    }
  }

  return(list(
    correlations = correlations,
    p_values = p_values
  ))
}

# 执行分析
cat("正在进行相关性分析...\n")
correlation_results <- perform_correlation_gwas(final_genotype_matrix, final_phenotypes)

cat(sprintf("完成了%d个SNP的分析\n", length(correlation_results$p_values)))
```

### 6.2 线性回归分析法

```r
# 方法2：线性回归分析
perform_regression_gwas <- function(genotypes, phenotype) {
  n_snps <- ncol(genotypes)

  # 存储结果
  beta_values <- numeric(n_snps)
  p_values <- numeric(n_snps)

  for(i in 1:n_snps) {
    snp_genotypes <- genotypes[, i]

    # 线性回归: phenotype ~ genotype
    if(length(unique(snp_genotypes)) > 1) {
      lm_result <- lm(phenotype ~ snp_genotypes)
      summary_result <- summary(lm_result)

      # 提取SNP的效应和p值
      beta_values[i] <- summary_result$coefficients[2, 1]  # 回归系数
      p_values[i] <- summary_result$coefficients[2, 4]     # p值
    } else {
      beta_values[i] <- 0
      p_values[i] <- 1
    }
  }

  return(list(
    beta = beta_values,
    p_values = p_values
  ))
}

# 执行回归分析
cat("正在进行线性回归分析...\n")
regression_results <- perform_regression_gwas(final_genotype_matrix, final_phenotypes)
```

---

## 7. 第五步：结果可视化

### 7.1 曼哈顿图

```r
# 创建曼哈顿图
create_manhattan_plot <- function(p_values, snp_info, title = "Manhattan Plot") {
  # 准备数据
  log_p <- -log10(p_values)
  chromosomes <- snp_info$Chromosome

  # 为不同染色体设置颜色
  colors <- c("darkblue", "darkgreen", "darkred", "darkorange", "purple",
              "brown", "pink", "gray", "olive", "cyan")
  chr_colors <- colors[chromosomes]

  # 绘制曼哈顿图
  plot(1:length(p_values), log_p,
       col = chr_colors,
       pch = 16, cex = 0.8,
       xlab = "SNP Position",
       ylab = "-log10(P-value)",
       main = title)

  # 添加显著性阈值线
  significance_threshold <- 0.05 / length(p_values)  # Bonferroni校正
  abline(h = -log10(significance_threshold), col = "red", lty = 2, lwd = 2)

  # 如果是模拟数据，标记真实QTN位置
  if(exists("qtn_indices")) {
    abline(v = qtn_indices, col = "black", lty = 3, lwd = 1)
    text(qtn_indices, max(log_p) * 0.9, "QTN", cex = 0.8, col = "black")
  }

  cat(sprintf("显著性阈值 (Bonferroni): %.2e\n", significance_threshold))
}

# 绘制相关性分析的曼哈顿图
cat("=== 绘制曼哈顿图 ===\n")
par(mfrow = c(2, 1))
create_manhattan_plot(correlation_results$p_values, snp_info_filtered,
                     "Manhattan Plot - Correlation Analysis")

create_manhattan_plot(regression_results$p_values, snp_info_filtered,
                     "Manhattan Plot - Regression Analysis")
```

### 7.2 QQ图检查

```r
# QQ图：检查p值分布是否正常
create_qq_plot <- function(p_values, title = "QQ Plot") {
  observed <- -log10(sort(p_values))
  expected <- -log10(ppoints(length(p_values)))

  plot(expected, observed,
       xlab = "Expected -log10(P)",
       ylab = "Observed -log10(P)",
       main = title,
       pch = 16, cex = 0.8)

  # 添加对角线
  abline(a = 0, b = 1, col = "red", lwd = 2)

  # 计算膨胀因子
  inflation_factor <- median(qchisq(1 - p_values, 1)) / qchisq(0.5, 1)
  text(max(expected) * 0.1, max(observed) * 0.9,
       sprintf("λ = %.3f", inflation_factor), cex = 1.2)

  cat(sprintf("基因组膨胀因子 λ = %.3f\n", inflation_factor))
}

# 绘制QQ图
par(mfrow = c(1, 2))
create_qq_plot(correlation_results$p_values, "QQ Plot - Correlation")
create_qq_plot(regression_results$p_values, "QQ Plot - Regression")
```

---

## 8. 第六步：统计质量控制

### 8.1 结果评估

```r
cat("=== 结果评估 ===\n")

# 找出最显著的SNP
top_snps_correlation <- order(correlation_results$p_values)[1:10]
top_snps_regression <- order(regression_results$p_values)[1:10]

cat("相关性分析top 10 SNP:\n")
for(i in 1:10) {
  snp_idx <- top_snps_correlation[i]
  cat(sprintf("  SNP %d: P = %.2e\n", snp_idx, correlation_results$p_values[snp_idx]))
}

cat("\n回归分析top 10 SNP:\n")
for(i in 1:10) {
  snp_idx <- top_snps_regression[i]
  cat(sprintf("  SNP %d: P = %.2e\n", snp_idx, regression_results$p_values[snp_idx]))
}

# 如果是模拟数据，检查是否找到了真实QTN
if(exists("qtn_indices")) {
  cat("\n=== 方法验证 ===\n")
  significance_threshold <- 0.05 / length(correlation_results$p_values)

  # 检查真实QTN的检出情况
  detected_qtns_corr <- sum(correlation_results$p_values[qtn_indices] < significance_threshold)
  detected_qtns_reg <- sum(regression_results$p_values[qtn_indices] < significance_threshold)

  cat(sprintf("真实QTN数量: %d\n", length(qtn_indices)))
  cat(sprintf("相关性分析检出: %d/%d\n", detected_qtns_corr, length(qtn_indices)))
  cat(sprintf("回归分析检出: %d/%d\n", detected_qtns_reg, length(qtn_indices)))

  # 假阳性检查
  false_positives_corr <- sum(correlation_results$p_values[-qtn_indices] < significance_threshold)
  false_positives_reg <- sum(regression_results$p_values[-qtn_indices] < significance_threshold)

  cat(sprintf("相关性分析假阳性: %d\n", false_positives_corr))
  cat(sprintf("回归分析假阳性: %d\n", false_positives_reg))
}
```

### 8.2 单个SNP效应分析

```r
# 分析最显著SNP的效应
analyze_top_snp <- function(genotypes, phenotype, snp_index, snp_name) {
  snp_genotypes <- genotypes[, snp_index]

  cat(sprintf("\n=== %s 效应分析 ===\n", snp_name))

  # 按基因型分组统计
  for(genotype in 0:2) {
    individuals <- which(snp_genotypes == genotype)
    if(length(individuals) > 0) {
      phenotype_subset <- phenotype[individuals]
      cat(sprintf("基因型 %d: n=%d, 平均表型=%.3f (±%.3f)\n",
                  genotype, length(individuals),
                  mean(phenotype_subset), sd(phenotype_subset)))
    }
  }

  # 绘制箱线图
  boxplot(phenotype ~ snp_genotypes,
          main = sprintf("%s Genotype Effect", snp_name),
          xlab = "Genotype (0=AA, 1=Aa, 2=aa)",
          ylab = "Phenotype Value")
}

# 分析top SNP
top_snp_idx <- which.min(regression_results$p_values)
top_snp_name <- colnames(final_genotype_matrix)[top_snp_idx]
analyze_top_snp(final_genotype_matrix, final_phenotypes, top_snp_idx, top_snp_name)
```

---

## 9. 第七步：改进的分析方法

### 9.1 主成分分析控制群体结构

```r
cat("=== 主成分分析 ===\n")

# 对基因型矩阵进行PCA
genotype_pca <- prcomp(final_genotype_matrix, center = TRUE, scale. = TRUE)

# 检查前几个主成分的方差解释比例
var_explained <- summary(genotype_pca)$importance[2, 1:10]
cat("前10个主成分的方差解释比例:\n")
for(i in 1:10) {
  cat(sprintf("  PC%d: %.3f\n", i, var_explained[i]))
}

# 绘制主成分散点图
plot(genotype_pca$x[, 1], genotype_pca$x[, 2],
     xlab = sprintf("PC1 (%.1f%%)", var_explained[1] * 100),
     ylab = sprintf("PC2 (%.1f%%)", var_explained[2] * 100),
     main = "Population Structure (PCA)",
     pch = 16, cex = 0.8)
```

### 9.2 控制群体结构的GLM分析

```r
# GLM分析：控制前几个主成分
perform_glm_gwas <- function(genotypes, phenotype, covariates, n_pcs = 3) {
  n_snps <- ncol(genotypes)

  # 存储结果
  beta_values <- numeric(n_snps)
  p_values <- numeric(n_snps)

  for(i in 1:n_snps) {
    snp_genotypes <- genotypes[, i]

    if(length(unique(snp_genotypes)) > 1) {
      # 构建设计矩阵：截距 + 主成分 + SNP
      design_matrix <- cbind(1, covariates[, 1:n_pcs], snp_genotypes)

      # 线性回归
      lm_result <- lm(phenotype ~ design_matrix - 1)  # -1是因为已经包含截距项
      summary_result <- summary(lm_result)

      # SNP效应是最后一个系数
      snp_coeff_idx <- ncol(design_matrix)
      beta_values[i] <- summary_result$coefficients[snp_coeff_idx, 1]
      p_values[i] <- summary_result$coefficients[snp_coeff_idx, 4]
    } else {
      beta_values[i] <- 0
      p_values[i] <- 1
    }
  }

  return(list(
    beta = beta_values,
    p_values = p_values
  ))
}

# 执行GLM分析
cat("正在进行GLM分析（控制群体结构）...\n")
glm_results <- perform_glm_gwas(final_genotype_matrix, final_phenotypes,
                               genotype_pca$x, n_pcs = 3)

# 比较不同方法的结果
par(mfrow = c(3, 1))
create_manhattan_plot(correlation_results$p_values, snp_info_filtered,
                     "Correlation Analysis")
create_manhattan_plot(regression_results$p_values, snp_info_filtered,
                     "Simple Regression")
create_manhattan_plot(glm_results$p_values, snp_info_filtered,
                     "GLM with Population Structure Control")
```

---

## 10. 完整代码示例

### 10.1 整合的分析脚本

```r
# ======================================================================
# 完整的GWAS分析流程
# ======================================================================

# 清理环境并载入数据
rm(list = ls())
setwd("/Users/your_path/statistical genomics")

# 载入数据
genotype_data <- read.table("data/mdp_numeric.txt", header = TRUE, sep = "\t")
snp_info <- read.table("data/mdp_SNP_information.txt", header = TRUE, sep = "\t")
trait_data <- read.table("data/mdp_traits.txt", header = TRUE, sep = "\t")

# 数据预处理
individual_names <- genotype_data[, 1]
genotype_matrix <- as.matrix(genotype_data[, -1])
rownames(genotype_matrix) <- individual_names

# 质量控制
polymorphic_snps <- apply(genotype_matrix, 2, function(x) length(unique(x[!is.na(x)])) > 1)
allele_freqs <- apply(genotype_matrix, 2, function(x) mean(x, na.rm = TRUE) / 2)
maf <- pmin(allele_freqs, 1 - allele_freqs)
good_snps <- polymorphic_snps & (maf >= 0.05)

filtered_genotype_matrix <- genotype_matrix[, good_snps]
filtered_snp_info <- snp_info[good_snps, ]

# 表型数据处理
selected_trait <- "EarHT"
genotype_individuals <- genotype_data[, 1]
trait_individuals <- trait_data[, 1]
common_individuals <- intersect(genotype_individuals, trait_individuals)

genotype_matched_idx <- match(common_individuals, genotype_individuals)
trait_matched_idx <- match(common_individuals, trait_individuals)

matched_genotype_data <- genotype_data[genotype_matched_idx, ]
matched_trait_data <- trait_data[trait_matched_idx, ]

trait_column_idx <- which(colnames(matched_trait_data) == selected_trait)
phenotype_values <- matched_trait_data[, trait_column_idx]
names(phenotype_values) <- matched_trait_data[, 1]

valid_individuals <- !is.na(phenotype_values)
final_phenotypes <- phenotype_values[valid_individuals]
final_genotype_data <- matched_genotype_data[valid_individuals, ]
final_genotype_matrix <- as.matrix(final_genotype_data[, -1])
rownames(final_genotype_matrix) <- final_genotype_data[, 1]

# 应用质量控制到最终基因型矩阵
final_genotype_matrix <- final_genotype_matrix[, good_snps]

# GWAS分析
gwas_results <- list()

# 1. 相关性分析
gwas_results$correlation <- apply(final_genotype_matrix, 2, function(x) {
  if(length(unique(x)) > 1) {
    cor.test(x, final_phenotypes)$p.value
  } else {
    1
  }
})

# 2. 简单回归
gwas_results$regression <- apply(final_genotype_matrix, 2, function(x) {
  if(length(unique(x)) > 1) {
    summary(lm(final_phenotypes ~ x))$coefficients[2, 4]
  } else {
    1
  }
})

# 3. GLM（控制群体结构）
pca_result <- prcomp(final_genotype_matrix, center = TRUE, scale. = TRUE)
gwas_results$glm <- apply(final_genotype_matrix, 2, function(x) {
  if(length(unique(x)) > 1) {
    design_matrix <- cbind(1, pca_result$x[, 1:3], x)
    summary(lm(final_phenotypes ~ design_matrix - 1))$coefficients[5, 4]
  } else {
    1
  }
})

# 结果可视化
par(mfrow = c(3, 1), mar = c(4, 4, 3, 1))

for(method in names(gwas_results)) {
  p_values <- gwas_results[[method]]
  log_p <- -log10(p_values)
  chromosomes <- filtered_snp_info$Chromosome
  colors <- rainbow(max(chromosomes))[chromosomes]

  plot(1:length(p_values), log_p, col = colors, pch = 16, cex = 0.8,
       xlab = "SNP Position", ylab = "-log10(P-value)",
       main = paste("GWAS Results -", method))

  # 显著性阈值
  abline(h = -log10(0.05 / length(p_values)), col = "red", lty = 2)
}

# 结果评估
cat("=== 真实数据GWAS分析结果 ===\n")
significance_threshold <- 0.05 / length(gwas_results$correlation)

for(method in names(gwas_results)) {
  p_values <- gwas_results[[method]]

  # 显著SNP数量
  significant_snps <- sum(p_values < significance_threshold)

  # 最显著SNP
  min_p_idx <- which.min(p_values)
  min_p_value <- p_values[min_p_idx]

  cat(sprintf("%s: 显著SNP数量 %d, 最小P值 %.2e (位置: %d)\n",
              method, significant_snps, min_p_value, min_p_idx))
}

cat("\n分析完成！\n")
```

---

## 总结

这个教程展示了从原始基因型数据和真实表型数据到GWAS结果的完整流程：

### 关键步骤
1. **数据理解**：明确数据格式和含义
2. **质量控制**：过滤低质量SNP
3. **表型处理**：载入和预处理真实表型数据
4. **统计分析**：从简单到复杂的关联分析方法
5. **结果解读**：可视化和结果评估

### 方法进阶
- **基础方法**：相关性分析、简单回归
- **改进方法**：GLM控制群体结构
- **高级方法**：MLM、贝叶斯方法等（后续学习）

### 实际应用建议
1. 理解真实数据的复杂性和处理方法
2. 逐步增加分析复杂度
3. 重视质量控制和结果验证
4. 理解方法的假设和局限性
5. 注意缺失值处理和个体匹配

这个框架为你提供了处理真实GWAS数据的坚实基础，可以根据具体需求进行扩展和改进。