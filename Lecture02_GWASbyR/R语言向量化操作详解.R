# R语言向量化比较操作详解
# 解答：myGM[, 2] < 6 的工作原理

# =============================================================================
# 第一部分：理解数据结构
# =============================================================================

# 首先，让我们看看myGM数据的结构
print("=== myGM数据前10行 ===")
# myGM数据示例：
#     SNP           Chromosome  Position
# 1   PZB00859.1    1          157104
# 2   PZA01271.1    1          1947984
# 3   PZA03613.2    1          2914066
# 4   ...           2          ...
# 5   ...           3          ...
# 6   ...           6          ...
# 7   ...           7          ...

# myGM[, 2] 表示提取第2列（Chromosome列）
print("=== 提取染色体列 myGM[, 2] ===")
# 假设myGM[, 2]的内容是这样的：
chromosome_column <- c(1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10)
print("染色体列的内容：")
print(chromosome_column)

# =============================================================================
# 第二部分：向量化比较操作的神奇之处
# =============================================================================

print("\n=== R语言的向量化比较 ===")
print("当我们写 chromosome_column < 6 时：")

# R语言会自动对向量中的每个元素进行比较
result <- chromosome_column < 6
print("结果：")
print(result)

print("\n详细解释每一步：")
for(i in 1:length(chromosome_column)) {
  cat(sprintf("chromosome_column[%d] = %d, %d < 6 = %s\n",
              i, chromosome_column[i], chromosome_column[i],
              chromosome_column[i] < 6))
}

# =============================================================================
# 第三部分：生物学意义解释
# =============================================================================

print("\n=== 生物学意义 ===")
print("在GWAS分析中，我们经常需要按染色体分组分析：")
print("- 训练集：1-5号染色体 (用于模拟QTN)")
print("- 测试集：6-10号染色体 (用于评估方法性能)")

print("\nindex1to5 <- myGM[, 2] < 6 的作用：")
print("- TRUE: 表示该SNP位于1-5号染色体")
print("- FALSE: 表示该SNP位于6-10号染色体")

# =============================================================================
# 第四部分：逻辑索引的实际应用
# =============================================================================

print("\n=== 逻辑索引的使用 ===")

# 模拟一些基因型数据
set.seed(123)
simulated_genotypes <- matrix(sample(0:2, 25*5, replace=TRUE), 25, 5)
colnames(simulated_genotypes) <- paste0("SNP", 1:5)

print("模拟的基因型数据（5个SNP）：")
print(head(simulated_genotypes, 3))

# 假设前3个SNP在1-5号染色体，后2个在6-10号染色体
chromosome_info <- c(1, 2, 3, 7, 8)
index1to5 <- chromosome_info < 6
print(sprintf("\n染色体信息: %s", paste(chromosome_info, collapse=", ")))
print(sprintf("index1to5: %s", paste(index1to5, collapse=", ")))

# 使用逻辑索引提取数据
X1to5 <- simulated_genotypes[, index1to5]
X6to10 <- simulated_genotypes[, !index1to5]  # !表示逻辑非

print("\n提取1-5号染色体的SNP:")
print(head(X1to5, 3))
print("\n提取6-10号染色体的SNP:")
print(head(X6to10, 3))

# =============================================================================
# 第五部分：为什么要这样做？
# =============================================================================

print("\n=== 为什么要分染色体分析？ ===")
print("1. 避免数据泄露：")
print("   - 在1-5号染色体上模拟QTN（已知答案）")
print("   - 在6-10号染色体上测试方法（未知答案）")
print("   - 这样可以客观评估GWAS方法的性能")

print("\n2. 控制实验：")
print("   - 训练集：包含真实信号")
print("   - 测试集：主要是噪音")
print("   - 通过QQ图检查假阳性率")

print("\n3. 模拟真实情况：")
print("   - 真实GWAS中，我们不知道QTN在哪里")
print("   - 通过这种分割，模拟盲测条件")

# =============================================================================
# 第六部分：其他向量化比较操作示例
# =============================================================================

print("\n=== 其他常用的向量化比较 ===")

ages <- c(18, 25, 30, 35, 40, 45, 50)
print(sprintf("年龄向量: %s", paste(ages, collapse=", ")))

print("\n各种比较操作：")
print(sprintf("ages >= 30: %s", paste(ages >= 30, collapse=", ")))
print(sprintf("ages == 35: %s", paste(ages == 35, collapse=", ")))
print(sprintf("ages != 25: %s", paste(ages != 25, collapse=", ")))

# 多条件组合
young_or_old <- (ages < 25) | (ages > 45)
print(sprintf("年轻(<25)或年老(>45): %s", paste(young_or_old, collapse=", ")))

middle_aged <- (ages >= 30) & (ages <= 40)
print(sprintf("中年(30-40): %s", paste(middle_aged, collapse=", ")))

# =============================================================================
# 第七部分：在数据框中的应用
# =============================================================================

print("\n=== 在数据框中筛选数据 ===")

# 创建示例数据框
sample_data <- data.frame(
  SNP = paste0("SNP", 1:8),
  Chromosome = c(1, 1, 2, 3, 5, 6, 8, 10),
  Position = c(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000),
  P_value = c(0.001, 0.05, 0.1, 0.02, 0.3, 0.08, 0.006, 0.4)
)

print("示例数据框：")
print(sample_data)

# 筛选1-5号染色体的SNP
early_chr <- sample_data[sample_data$Chromosome < 6, ]
print("\n1-5号染色体的SNP：")
print(early_chr)

# 筛选显著的SNP (p < 0.05)
significant_snps <- sample_data[sample_data$P_value < 0.05, ]
print("\n显著的SNP (p < 0.05)：")
print(significant_snps)

# 组合条件
early_and_sig <- sample_data[sample_data$Chromosome < 6 & sample_data$P_value < 0.05, ]
print("\n1-5号染色体且显著的SNP：")
print(early_and_sig)

# =============================================================================
# 总结
# =============================================================================

print("\n=== 总结 ===")
print("myGM[, 2] < 6 这个操作：")
print("1. myGM[, 2] 提取第2列（染色体编号）")
print("2. < 6 对每个元素进行比较")
print("3. 返回逻辑向量（TRUE/FALSE）")
print("4. TRUE表示染色体编号小于6（即1-5号染色体）")
print("5. 这个逻辑向量用作索引，选择特定的SNP子集")

print("\n这是R语言向量化操作的强大之处：")
print("- 一行代码完成整个向量的元素级操作")
print("- 比写循环更简洁、更高效")
print("- 是R语言数据分析的核心特性")