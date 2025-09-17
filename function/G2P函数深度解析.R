# G2P函数深度解析：从基因型到表型的模拟
# Genotype to Phenotype Simulation - 核心教学函数详解

# =============================================================================
# 第一部分：G2P函数的重要性和教学价值
# =============================================================================

# 【为什么需要G2P函数？】
#
# 1. **GWAS教学的核心挑战**：
#    - 真实GWAS数据：我们不知道哪些基因真正影响性状
#    - 无法验证方法：无法客观评估GWAS方法的准确性
#    - 黑盒问题：学生难以理解"已知答案"的分析过程
#
# 2. **G2P函数的解决方案**：
#    - 创建"已知答案"的模拟数据
#    - 控制遗传参数，理解其对表型的影响
#    - 客观评估GWAS方法的性能
#    - 教学演示：从基因型如何产生表型
#
# 3. **教学优势**：
#    - 可重现：设置随机种子确保结果一致
#    - 可控制：调整参数观察不同遗传情景
#    - 可验证：已知QTN位置，可检验方法准确性

# =============================================================================
# 第二部分：核心生物学概念详解
# =============================================================================

print("=== 核心生物学概念 ===")

# 【1. QTN (Quantitative Trait Nucleotide)】
print("1. QTN - 数量性状核苷酸")
print("   定义：真正影响数量性状的DNA变异位点")
print("   特点：")
print("   - 与疾病基因不同，QTN影响连续性状（如身高、体重）")
print("   - 通常效应较小，需要大样本才能检测到")
print("   - 多个QTN共同作用决定性状表现")
print("   - 在真实GWAS中，我们试图找到这些QTN")

# 【2. 遗传力 (Heritability, h²)】
print("\n2. 遗传力 - 遗传因素的贡献度")
print("   定义：遗传变异占总表型变异的比例")
print("   公式：h² = V_G / V_P = V_G / (V_G + V_E)")
print("   其中：")
print("   - V_G: 遗传方差（基因效应导致的变异）")
print("   - V_E: 环境方差（环境因素导致的变异）")
print("   - V_P: 表型方差（总变异）")
print("   解释：")
print("   - h² = 0.8: 80%的表型差异由遗传因素决定")
print("   - h² = 0.2: 20%的表型差异由遗传因素决定")

# 【3. 加性效应模型】
print("\n3. 加性效应模型")
print("   基本假设：多个基因的效应可以简单相加")
print("   模型：表型 = Σ(基因效应) + 环境效应")
print("   基因型编码：")
print("   - 0 (AA): 参考纯合子，效应为 0 × β")
print("   - 1 (Aa): 杂合子，效应为 1 × β")
print("   - 2 (aa): 变异纯合子，效应为 2 × β")
print("   其中β是该位点的加性效应值")

# =============================================================================
# 第三部分：G2P函数参数详解
# =============================================================================

print("\n=== G2P函数参数详解 ===")

# 【输入参数】
print("输入参数：")
print("1. X: 基因型矩阵")
print("   - 行：个体（281个玉米品系）")
print("   - 列：SNP标记（筛选后的SNP）")
print("   - 值：0/1/2编码的基因型")

print("\n2. h2: 遗传力 (0-1之间)")
print("   - 0.75表示75%的表型变异由遗传决定")
print("   - 影响信噪比：h2越高，信号越强")

print("\n3. alpha: 分布参数")
print("   - 在几何分布中控制效应衰减速度")
print("   - alpha接近1：效应衰减慢")
print("   - alpha接近0：效应衰减快")

print("\n4. NQTN: QTN数量")
print("   - 模拟多少个SNP真正影响性状")
print("   - 数量越多，每个QTN平均效应越小")

print("\n5. distribution: 效应分布")
print("   - 'norm': 正态分布（对称，大部分效应接近平均值）")
print("   - 'geom': 几何分布（少数大效应，多数小效应）")

# 【输出结果】
print("\n输出结果：")
print("1. y: 模拟的表型值")
print("2. add: 加性遗传效应")
print("3. residual: 环境/随机效应")
print("4. QTN.position: QTN在基因组中的位置")
print("5. addeffect: 每个QTN的效应大小")
print("6. SNPQ: QTN的基因型数据")

# =============================================================================
# 第四部分：算法实现逐步解析
# =============================================================================

print("\n=== 算法实现逐步解析 ===")

# 【步骤1：随机选择QTN】
print("步骤1：随机选择QTN位点")
print("代码：QTN.position <- sample(m, NQTN, replace = FALSE)")
print("说明：")
print("- 从所有SNP中随机选择NQTN个作为真正的QTN")
print("- replace=FALSE确保不重复选择")
print("- 模拟真实情况：我们不知道哪些SNP是真正的致病变异")

# 模拟示例
set.seed(123)
total_snps <- 1000
n_qtn <- 10
example_qtn_pos <- sample(total_snps, n_qtn, replace = FALSE)
print(sprintf("示例：从%d个SNP中选择%d个QTN", total_snps, n_qtn))
print(sprintf("选中的QTN位置：%s", paste(example_qtn_pos, collapse=", ")))

# 【步骤2：生成QTN效应】
print("\n步骤2：生成QTN效应值")

# 正态分布示例
print("正态分布效应：")
set.seed(123)
norm_effects <- rnorm(n_qtn, 0, 1)
print(sprintf("效应值：%s", paste(round(norm_effects, 3), collapse=", ")))
print("特点：效应值对称分布，正负效应约各占一半")

# 几何分布示例
print("\n几何分布效应：")
alpha <- 0.7
geom_effects <- alpha^(1:n_qtn)
print(sprintf("alpha=%.1f时的效应值：%s", alpha, paste(round(geom_effects, 3), collapse=", ")))
print("特点：少数大效应，多数小效应，符合'寡基因模型'")

# 【步骤3：计算遗传值】
print("\n步骤3：计算每个个体的遗传值")
print("代码：effect <- SNPQ %*% addeffect")
print("说明：矩阵乘法计算个体的总遗传效应")

# 模拟基因型矩阵
set.seed(123)
n_individuals <- 5
example_genotypes <- matrix(sample(0:2, n_individuals * n_qtn, replace=TRUE),
                           nrow=n_individuals, ncol=n_qtn)
colnames(example_genotypes) <- paste0("QTN", 1:n_qtn)
rownames(example_genotypes) <- paste0("Ind", 1:n_individuals)

print("示例基因型矩阵（5个体×10个QTN）：")
print(example_genotypes[, 1:5])  # 只显示前5个QTN

# 计算遗传值
genetic_values <- example_genotypes %*% norm_effects
print("\n每个个体的遗传值：")
for(i in 1:n_individuals) {
  cat(sprintf("个体%d: 遗传值 = %.3f\n", i, genetic_values[i]))
}

# 【步骤4：添加环境效应】
print("\n步骤4：根据遗传力添加环境效应")
print("目标：确保遗传方差占总方差的h²比例")

h2 <- 0.75
genetic_var <- var(genetic_values)
environmental_var <- (genetic_var - h2 * genetic_var) / h2
environmental_effects <- rnorm(n_individuals, 0, sqrt(environmental_var))

print(sprintf("遗传方差：%.6f", genetic_var))
print(sprintf("环境方差：%.6f", environmental_var))
print(sprintf("遗传力：%.3f", genetic_var / (genetic_var + environmental_var)))

# 【步骤5：生成最终表型】
final_phenotypes <- genetic_values + environmental_effects
print("\n最终表型 = 遗传值 + 环境效应：")
for(i in 1:n_individuals) {
  cat(sprintf("个体%d: %.3f = %.3f + %.3f\n",
              i, final_phenotypes[i], genetic_values[i], environmental_effects[i]))
}

# =============================================================================
# 第五部分：不同参数设置的生物学意义
# =============================================================================

print("\n=== 不同参数设置的生物学意义 ===")

# 【遗传力的影响】
print("1. 遗传力的影响")
h2_values <- c(0.2, 0.5, 0.8)
for(h2 in h2_values) {
  signal_to_noise <- h2 / (1 - h2)
  print(sprintf("h² = %.1f: 信噪比 = %.2f, GWAS检测%s",
                h2, signal_to_noise,
                ifelse(h2 > 0.6, "容易", ifelse(h2 > 0.3, "中等", "困难"))))
}

# 【QTN数量的影响】
print("\n2. QTN数量的影响")
nqtn_values <- c(5, 20, 100)
total_var <- 1.0  # 假设总遗传方差为1
for(nqtn in nqtn_values) {
  avg_effect <- total_var / nqtn
  print(sprintf("NQTN = %d: 平均效应 = %.3f, 检测难度%s",
                nqtn, avg_effect,
                ifelse(nqtn < 10, "低", ifelse(nqtn < 50, "中", "高"))))
}

# 【效应分布的影响】
print("\n3. 效应分布的影响")
print("正态分布：")
print("- 生物学意义：多基因效应，每个基因贡献相近")
print("- 实例：身高、体重等复杂性状")
print("- GWAS挑战：需要大样本检测小效应")

print("\n几何分布：")
print("- 生物学意义：少数主效基因+多个微效基因")
print("- 实例：某些疾病易感性、品质性状")
print("- GWAS优势：主效基因容易检测")

# =============================================================================
# 第六部分：G2P在GWAS教学中的应用
# =============================================================================

print("\n=== G2P在GWAS教学中的应用 ===")

print("1. 方法验证：")
print("   - 已知QTN位置，可检验GWAS是否找到正确答案")
print("   - 计算真阳性率、假阳性率等性能指标")

print("\n2. 参数效应研究：")
print("   - 不同h²下GWAS的检测能力")
print("   - 样本量对检测功效的影响")
print("   - 不同遗传结构的分析挑战")

print("\n3. 方法比较：")
print("   - 简单方法vs复杂方法的性能对比")
print("   - 在不同场景下的适用性")

print("\n4. 概念理解：")
print("   - 遗传力的直观解释")
print("   - 多基因效应的可视化")
print("   - 环境因素的作用机制")

# =============================================================================
# 第七部分：实际使用示例和常见问题
# =============================================================================

print("\n=== 实际使用示例 ===")

# 【基础调用】
print("基础调用示例：")
print("mySim <- G2P(X = X1to5, h2 = 0.75, alpha = 1, NQTN = 10, distribution = 'norm')")
print("解释：")
print("- 使用1-5号染色体的SNP")
print("- 遗传力75%")
print("- 10个QTN，正态分布效应")

# 【结果解读】
print("\n结果解读：")
print("- mySim$y: 用于GWAS分析的表型")
print("- mySim$QTN.position: 验证GWAS结果的'标准答案'")
print("- mySim$add: 理解遗传贡献")
print("- mySim$residual: 理解环境贡献")

# 【常见问题】
print("\n常见问题及解决：")

print("Q1: 为什么要设置random seed？")
print("A1: 确保每次运行得到相同结果，便于教学演示和结果重现")

print("\nQ2: 为什么只用1-5号染色体模拟QTN？")
print("A2: 6-10号染色体作为'干净'的测试集，评估假阳性率")

print("\nQ3: 如何选择合适的h²值？")
print("A3: 参考真实性状：身高~0.8，BMI~0.7，血压~0.4")

print("\nQ4: NQTN设置多少合适？")
print("A4: 教学用途：5-20个；接近真实：50-100个")

# =============================================================================
# 第八部分：扩展思考和高级应用
# =============================================================================

print("\n=== 扩展思考 ===")

print("1. 模型局限性：")
print("   - 只考虑加性效应，忽略上位性")
print("   - 假设QTN效应独立，实际可能存在连锁")
print("   - 环境效应假设为随机噪声")

print("\n2. 真实世界的复杂性：")
print("   - 基因×环境互作")
print("   - 基因×基因互作")
print("   - 群体分层和家系结构")

print("\n3. 高级扩展：")
print("   - 加入显性效应模拟")
print("   - 模拟群体结构影响")
print("   - 多性状联合分析")

print("\n4. 教学启发：")
print("   - 理解模型假设的重要性")
print("   - 认识统计模型与生物现实的差距")
print("   - 培养批判性思维")

# =============================================================================
# 总结
# =============================================================================

print("\n=== 总结 ===")
print("G2P函数是GWAS教学的核心工具，它：")
print("1. 解决了'已知答案'的问题，使GWAS方法验证成为可能")
print("2. 将抽象的遗传学概念转化为可操作的数值模拟")
print("3. 提供了探索不同遗传参数对GWAS结果影响的平台")
print("4. 帮助学生理解从基因型到表型的分子基础")
print("5. 为更高级的GWAS方法学习奠定坚实基础")

print("\n这个函数体现了生物信息学教学的精髓：")
print("- 理论与实践相结合")
print("- 复杂概念的简化建模")
print("- 可重现的科学研究方法")
print("- 批判性思维的培养")