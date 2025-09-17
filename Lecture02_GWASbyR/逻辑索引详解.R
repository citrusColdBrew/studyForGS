# R语言逻辑索引详解：X[, index1to5] 的工作原理
# 解答：如何用TRUE/FALSE向量筛选矩阵的列

# =============================================================================
# 第一部分：逻辑索引的基本概念
# =============================================================================

print("=== R语言逻辑索引机制 ===")
print("逻辑索引是R语言的核心特性之一")
print("用TRUE/FALSE向量来选择数据的子集")

# =============================================================================
# 第二部分：从简单向量开始理解
# =============================================================================

print("\n=== 从向量开始理解 ===")

# 创建一个简单的向量
snp_names <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6")
chromosomes <- c(1, 1, 2, 3, 6, 7)

print("SNP名称:")
print(snp_names)
print("对应染色体:")
print(chromosomes)

# 创建逻辑索引
logical_index <- chromosomes < 6
print("\n逻辑索引 (chromosomes < 6):")
print(logical_index)

# 使用逻辑索引筛选
selected_snps <- snp_names[logical_index]
print("\n筛选出的SNP (只保留TRUE位置的元素):")
print(selected_snps)

print("\n详细对照:")
for(i in 1:length(snp_names)) {
  cat(sprintf("位置%d: %s, 染色体%d, 逻辑值%s, %s\n",
              i, snp_names[i], chromosomes[i], logical_index[i],
              ifelse(logical_index[i], "被选中", "被过滤")))
}

# =============================================================================
# 第三部分：矩阵的列筛选 - 核心概念
# =============================================================================

print("\n=== 矩阵列筛选的工作原理 ===")

# 创建示例基因型矩阵 (行=个体, 列=SNP)
set.seed(123)
genotype_matrix <- matrix(sample(0:2, 6*6, replace=TRUE), nrow=6, ncol=6)
colnames(genotype_matrix) <- paste0("SNP", 1:6)
rownames(genotype_matrix) <- paste0("Individual", 1:6)

print("原始基因型矩阵 X:")
print(genotype_matrix)

# 对应的染色体信息
chromosome_info <- c(1, 1, 2, 3, 6, 7)
print(sprintf("\n每个SNP的染色体: %s", paste(chromosome_info, collapse=", ")))

# 创建逻辑索引
index1to5 <- chromosome_info < 6
print(sprintf("逻辑索引 index1to5: %s", paste(index1to5, collapse=", ")))

print("\n解释逻辑索引的含义:")
for(i in 1:length(index1to5)) {
  cat(sprintf("第%d列 (%s): 染色体%d, 逻辑值%s -> %s\n",
              i, colnames(genotype_matrix)[i], chromosome_info[i],
              index1to5[i], ifelse(index1to5[i], "保留此列", "删除此列")))
}

# =============================================================================
# 第四部分：执行列筛选 X[, index1to5]
# =============================================================================

print("\n=== 执行 X[, index1to5] ===")

# 使用逻辑索引筛选列
X1to5 <- genotype_matrix[, index1to5]

print("筛选后的矩阵 X1to5 (只保留1-5号染色体的SNP):")
print(X1to5)

print("\n对比原始矩阵和筛选后矩阵:")
print("原始矩阵的列名:")
print(colnames(genotype_matrix))
print("筛选后矩阵的列名:")
print(colnames(X1to5))

# =============================================================================
# 第五部分：矩阵索引语法详解
# =============================================================================

print("\n=== 矩阵索引语法 matrix[行索引, 列索引] ===")

print("R语言矩阵索引的通用格式: matrix[row_index, col_index]")
print("- 行索引：选择哪些行")
print("- 列索引：选择哪些列")
print("- 空白表示选择全部")

# 各种索引方式示例
print("\n各种索引方式:")

# 1. 选择特定行和列
print("1. genotype_matrix[1:3, 1:2] - 前3行，前2列:")
print(genotype_matrix[1:3, 1:2])

# 2. 选择所有行，特定列
print("\n2. genotype_matrix[, c(1,3,5)] - 所有行，第1,3,5列:")
print(genotype_matrix[, c(1,3,5)])

# 3. 使用逻辑索引选择列
print("\n3. genotype_matrix[, index1to5] - 所有行，逻辑索引选择的列:")
print(genotype_matrix[, index1to5])

# 4. 使用逻辑非选择相反的列
index6to10 <- !index1to5  # 逻辑非操作
print(sprintf("\n逻辑非 !index1to5: %s", paste(index6to10, collapse=", ")))
print("4. genotype_matrix[, !index1to5] - 6-10号染色体:")
print(genotype_matrix[, !index1to5])

# =============================================================================
# 第六部分：实际GWAS中的应用场景
# =============================================================================

print("\n=== 在GWAS中的实际应用 ===")

# 模拟更大的数据集
set.seed(456)
n_individuals <- 281  # 玉米数据集的个体数
n_snps <- 20         # 简化的SNP数量
large_genotype <- matrix(sample(0:2, n_individuals * n_snps, replace=TRUE),
                        nrow=n_individuals, ncol=n_snps)

# 模拟染色体分布
simulated_chromosomes <- c(rep(1, 3), rep(2, 3), rep(3, 2), rep(4, 2),
                          rep(5, 2), rep(6, 2), rep(7, 2), rep(8, 2), rep(9, 2), rep(10, 2))

print("模拟的大数据集:")
print(sprintf("个体数: %d", nrow(large_genotype)))
print(sprintf("SNP数: %d", ncol(large_genotype)))
print(sprintf("染色体分布: %s", paste(simulated_chromosomes, collapse=", ")))

# 执行筛选
training_index <- simulated_chromosomes < 6
testing_index <- !training_index

training_set <- large_genotype[, training_index]
testing_set <- large_genotype[, testing_index]

print(sprintf("\n训练集 (1-5号染色体): %d个体 × %d个SNP",
              nrow(training_set), ncol(training_set)))
print(sprintf("测试集 (6-10号染色体): %d个体 × %d个SNP",
              nrow(testing_set), ncol(testing_set)))

# =============================================================================
# 第七部分：逻辑索引的高级应用
# =============================================================================

print("\n=== 逻辑索引的高级应用 ===")

# 1. 多条件筛选
print("1. 多条件筛选:")
# 选择1-3号染色体
early_chr <- simulated_chromosomes <= 3
print(sprintf("1-3号染色体: %s", paste(which(early_chr), collapse=", ")))

# 选择奇数染色体
odd_chr <- simulated_chromosomes %% 2 == 1
print(sprintf("奇数染色体: %s", paste(which(odd_chr), collapse=", ")))

# 组合条件：1-3号且奇数染色体
combined <- early_chr & odd_chr
print(sprintf("1-3号且奇数染色体: %s", paste(which(combined), collapse=", ")))

# 2. 基于p值筛选显著SNP
print("\n2. 基于p值筛选显著SNP:")
simulated_pvalues <- runif(n_snps, 0, 1)  # 模拟p值
significant <- simulated_pvalues < 0.05

print(sprintf("模拟的p值: %s",
              paste(round(simulated_pvalues, 3), collapse=", ")))
print(sprintf("显著性筛选 (p < 0.05): %s",
              paste(significant, collapse=", ")))
print(sprintf("显著SNP的位置: %s",
              paste(which(significant), collapse=", ")))

# =============================================================================
# 第八部分：常见错误和注意事项
# =============================================================================

print("\n=== 常见错误和注意事项 ===")

print("1. 维度匹配:")
print("   逻辑向量的长度必须等于要筛选的维度")
print(sprintf("   矩阵列数: %d", ncol(large_genotype)))
print(sprintf("   逻辑向量长度: %d", length(training_index)))

print("\n2. 数据类型:")
wrong_index <- c("TRUE", "FALSE", "TRUE")  # 字符串，不是逻辑值
correct_index <- c(TRUE, FALSE, TRUE)      # 逻辑值

print("   错误: 使用字符串 c('TRUE', 'FALSE')")
print("   正确: 使用逻辑值 c(TRUE, FALSE)")

print("\n3. 空结果处理:")
empty_filter <- simulated_chromosomes > 100  # 没有染色体号>100
empty_result <- large_genotype[, empty_filter]
print(sprintf("   空筛选结果的维度: %d × %d",
              nrow(empty_result), ncol(empty_result)))

# =============================================================================
# 总结
# =============================================================================

print("\n=== 总结 ===")
print("X[, index1to5] 的工作原理:")
print("1. X 是基因型矩阵 (281个体 × 3093个SNP)")
print("2. index1to5 是逻辑向量 (长度=3093，TRUE/FALSE)")
print("3. [, index1to5] 表示选择所有行，但只选择TRUE对应的列")
print("4. 结果是只包含1-5号染色体SNP的子矩阵")

print("\n筛选过程:")
print("原始矩阵: 281 × 3093")
print("逻辑索引: [T,T,T,...,F,F,F] (长度3093)")
print("筛选结果: 281 × (TRUE的数量)")

print("\n生物学意义:")
print("- 将基因组分为训练集(1-5号染色体)和测试集(6-10号染色体)")
print("- 在训练集模拟QTN，在测试集评估方法性能")
print("- 这是GWAS方法验证的标准策略")