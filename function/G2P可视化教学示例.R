# G2P函数可视化教学示例
# 帮助理解从基因型到表型的模拟过程

# =============================================================================
# 简化版G2P演示
# =============================================================================

print("=== G2P函数工作流程可视化演示 ===")

# 设置参数
set.seed(12345)  # 确保结果可重现
n_individuals <- 8    # 8个个体（简化演示）
n_snps <- 20         # 20个SNP
n_qtn <- 3           # 3个QTN
h2 <- 0.75           # 遗传力75%

# 步骤1：创建基因型矩阵
print("\n=== 步骤1：基因型数据 ===")
genotype_matrix <- matrix(sample(0:2, n_individuals * n_snps, replace=TRUE),
                         nrow=n_individuals, ncol=n_snps)
rownames(genotype_matrix) <- paste0("个体", 1:n_individuals)
colnames(genotype_matrix) <- paste0("SNP", 1:n_snps)

print("基因型矩阵（部分显示）：")
print(genotype_matrix[, 1:8])

# 步骤2：随机选择QTN
print("\n=== 步骤2：选择QTN ===")
qtn_positions <- sample(n_snps, n_qtn, replace=FALSE)
qtn_positions <- sort(qtn_positions)  # 排序便于显示
print(sprintf("从%d个SNP中选择的%d个QTN位置：%s",
              n_snps, n_qtn, paste(qtn_positions, collapse=", ")))

# 提取QTN基因型
qtn_genotypes <- genotype_matrix[, qtn_positions]
colnames(qtn_genotypes) <- paste0("QTN", 1:n_qtn)
print("\nQTN基因型数据：")
print(qtn_genotypes)

# 步骤3：生成QTN效应
print("\n=== 步骤3：QTN效应值 ===")
qtn_effects <- rnorm(n_qtn, 0, 1)
names(qtn_effects) <- paste0("QTN", 1:n_qtn)
print("每个QTN的加性效应：")
for(i in 1:n_qtn) {
  cat(sprintf("QTN%d (SNP%d): 效应 = %.3f\n",
              i, qtn_positions[i], qtn_effects[i]))
}

# 步骤4：计算遗传值
print("\n=== 步骤4：计算遗传值 ===")
print("遗传值 = Σ(基因型 × 效应)")

genetic_values <- qtn_genotypes %*% qtn_effects
print("\n详细计算过程：")
for(i in 1:n_individuals) {
  calculation_parts <- paste(sprintf("%d×%.3f",
                                   qtn_genotypes[i,], qtn_effects),
                           collapse=" + ")
  cat(sprintf("个体%d: %s = %.3f\n",
              i, calculation_parts, genetic_values[i]))
}

# 步骤5：计算环境效应
print("\n=== 步骤5：添加环境效应 ===")
genetic_var <- var(genetic_values)
environmental_var <- (genetic_var * (1 - h2)) / h2
environmental_effects <- rnorm(n_individuals, 0, sqrt(environmental_var))

print(sprintf("遗传方差: %.6f", genetic_var))
print(sprintf("环境方差: %.6f", environmental_var))
print(sprintf("理论遗传力: %.3f", h2))

# 步骤6：生成最终表型
print("\n=== 步骤6：最终表型 ===")
final_phenotypes <- genetic_values + environmental_effects
actual_h2 <- genetic_var / (genetic_var + var(environmental_effects))
print(sprintf("实际遗传力: %.3f", actual_h2))

print("\n表型组成分解：")
phenotype_data <- data.frame(
  个体 = paste0("个体", 1:n_individuals),
  遗传值 = round(genetic_values, 3),
  环境效应 = round(environmental_effects, 3),
  最终表型 = round(final_phenotypes, 3)
)
print(phenotype_data)

# =============================================================================
# 可视化效果
# =============================================================================

print("\n=== 可视化分析 ===")

# 遗传值vs表型的关系
correlation <- cor(genetic_values, final_phenotypes)
print(sprintf("遗传值与表型的相关性: %.3f", correlation))
print(sprintf("遗传力的平方根: %.3f", sqrt(h2)))
print("理论上，相关性应接近遗传力的平方根")

# 不同基因型的表型分布示例
print("\n=== QTN效应验证 ===")
# 选择效应最大的QTN进行分析
max_effect_qtn <- which.max(abs(qtn_effects))
qtn_name <- paste0("QTN", max_effect_qtn)
print(sprintf("分析效应最大的%s (效应=%.3f)：", qtn_name, qtn_effects[max_effect_qtn]))

# 按基因型分组显示表型
for(genotype in 0:2) {
  individuals_with_genotype <- which(qtn_genotypes[, max_effect_qtn] == genotype)
  if(length(individuals_with_genotype) > 0) {
    phenotypes_for_genotype <- final_phenotypes[individuals_with_genotype]
    mean_phenotype <- mean(phenotypes_for_genotype)
    cat(sprintf("基因型%d: 个体数=%d, 平均表型=%.3f\n",
                genotype, length(individuals_with_genotype), mean_phenotype))
  }
}

# =============================================================================
# 教学要点总结
# =============================================================================

print("\n=== 教学要点总结 ===")

print("1. G2P函数模拟了从基因到性状的生物学过程：")
print("   基因型 → 遗传效应 → 表型")

print("\n2. 关键概念的数值体现：")
print(sprintf("   - 遗传力(%.1f): %d%%的表型变异来自遗传因素",
              h2, round(h2*100)))
print(sprintf("   - 多基因效应: %d个QTN共同影响性状", n_qtn))
print("   - 加性模型: 基因效应简单相加")

print("\n3. GWAS分析的'已知答案'：")
print(sprintf("   - 真正的QTN位置: %s", paste(qtn_positions, collapse=", ")))
print("   - 可以验证GWAS方法是否找到正确答案")

print("\n4. 参数调节的生物学意义：")
print("   - 增加h²: 遗传信号更强，GWAS更容易成功")
print("   - 增加NQTN: 每个基因效应更小，检测更困难")
print("   - 不同效应分布: 模拟不同的遗传结构")

print("\n5. 为什么这个函数对GWAS教学如此重要：")
print("   - 提供了可控的实验环境")
print("   - 让抽象概念变得具体可测")
print("   - 支持方法验证和比较")
print("   - 帮助理解真实GWAS的挑战和局限")

# =============================================================================
# 与真实GWAS的对比
# =============================================================================

print("\n=== 与真实GWAS的对比 ===")

print("模拟数据的优势：")
print("✓ 已知QTN位置和效应")
print("✓ 可控制遗传参数")
print("✓ 可重现结果")
print("✓ 方便教学演示")

print("\n真实数据的复杂性：")
print("× 不知道真正的致病基因")
print("× 多种混杂因素")
print("× 群体结构和家系效应")
print("× 基因×环境互作")

print("\n教学策略：")
print("1. 先用G2P建立基本概念")
print("2. 理解理想情况下的GWAS")
print("3. 逐步加入复杂因素")
print("4. 最终分析真实数据")

print("\n这种从简单到复杂的教学路径，帮助学生：")
print("- 建立扎实的理论基础")
print("- 理解方法的适用条件")
print("- 培养解决实际问题的能力")
print("- 发展批判性思维")