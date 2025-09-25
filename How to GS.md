# How to GS：从零开始的全基因组选择教程

## 目录
1. [GS是什么？](#1-gs是什么)
2. [数据准备和理解](#2-数据准备和理解)
3. [第一步：数据载入和预处理](#3-第一步数据载入和预处理)
4. [第二步：训练测试集划分](#4-第二步训练测试集划分)
5. [第三步：线性混合模型方法](#5-第三步线性混合模型方法)
6. [第四步：贝叶斯回归方法](#6-第四步贝叶斯回归方法)
7. [第五步：神经网络方法](#7-第五步神经网络方法)
8. [第六步：方法比较和评估](#8-第六步方法比较和评估)
9. [第七步：交叉验证策略](#9-第七步交叉验证策略)
10. [完整代码示例](#10-完整代码示例)

---

## 1. GS是什么？

### 1.1 生物学背景与历史发展

**GS（Genomic Selection，全基因组选择）** 是利用全基因组标记信息预测个体基因组估计育种值（GEBV）的方法。

#### 传统育种的挑战
在传统育种中，我们面临以下问题：
- **表型测定时间长**：很多重要性状（如木材品质、抗病性）需要植物成熟后才能测定
- **测定成本高**：大规模表型调查需要大量人力物力
- **环境干扰**：表型受环境影响大，遗传与环境效应难以分离
- **世代间隔长**：传统选择需要等待个体成熟，育种周期漫长

#### GS的革命性突破
**核心理念**：
- **从"看表型"到"看基因型"**：传统选择依赖表型观测，GS直接利用DNA信息
- **从"等成熟"到"早预测"**：通过幼苗基因型预测成体表现
- **从"少标记"到"全基因组"**：利用覆盖全基因组的密集标记
- **从"找基因"到"算总分"**：不寻找单个大效应基因，而是计算所有基因的累积贡献

#### 发展历程
- **2001年**: Meuwissen等人首次提出GS概念
- **2008年**: 第一个商业化GS应用（奶牛育种）
- **2010年代**: 测序成本下降，GS在植物育种中广泛应用
- **现在**: 成为现代分子育种的标准工具

### 1.2 GS与相关概念的区别

**【GS vs GWAS】**：

| 方面 | GWAS | GS |
|------|------|-----|
| **目标** | 发现与性状关联的基因/QTN | 预测个体的基因组育种值（GEBV） |
| **关注点** | 少数大效应位点的显著性 | 所有标记的累积效应 |
| **统计策略** | 假设检验，控制假阳性 | 预测建模，最大化预测精度 |
| **结果解释** | 哪些基因影响性状 | 个体的遗传潜力有多大 |
| **应用** | 基因发掘，分子标记开发 | 育种选择，杂交组合预测 |

**【GS vs MAS（分子标记辅助选择）】**：
- **MAS**: 使用少数几个与目标性状紧密连锁的标记
- **GS**: 使用覆盖全基因组的大量标记，不依赖连锁关系

**为什么GS比MAS更有效？**
1. **捕捉更多遗传变异**：大多数复杂性状受多基因控制，单个标记难以解释全部遗传变异
2. **减少连锁拖累**：不依赖特定标记与QTN的连锁，避免重组破坏连锁关系
3. **适应性更强**：在不同群体和环境中都能保持较好的预测效果

### 1.3 GS的生物学原理

#### 数量遗传学基础
**基本模型**：P = G + E
- **P（表型值）**：我们观察到的性状表现
- **G（基因型值）**：由基因决定的遗传潜力
- **E（环境效应）**：环境因素的影响

**GS的假设**：
1. **多基因假说**：复杂性状受众多基因影响，每个基因效应相对较小
2. **加性效应**：基因效应可以简单相加（忽略上位性互作）
3. **连锁不平衡**：标记与QTN之间存在关联（即使不知道具体QTN位置）

#### 为什么全基因组标记有效？
想象基因组是一条项链，QTN（真正的致病基因）是项链上的珍珠：
- **传统方法**：只看几颗明显的珍珠（已知基因）
- **GS方法**：看整条项链的图案（全基因组标记）

即使我们不知道珍珠的确切位置，但项链的整体图案能反映珍珠的分布！

### 1.4 GS的优势与挑战

#### 优势
1. **时间效率**：
   - 🌱 **幼苗选择**：DNA提取后即可预测，无需等待表型成熟
   - ⏰ **缩短世代间隔**：从10年缩短到1-2年（如森林树木）

2. **成本效率**：
   - 💰 **降低表型测定成本**：测序成本 < 大规模田间试验成本
   - 📊 **提高选择强度**：可以从更大的候选群体中选择

3. **选择精度**：
   - 🎯 **环境影响小**：基因型不受环境影响
   - 📈 **累积效应**：考虑所有基因的贡献，不遗漏小效应基因

#### 挑战
1. **技术挑战**：
   - 需要高质量的基因型数据
   - 需要强大的计算资源和统计方法

2. **生物学挑战**：
   - 不同环境下基因表达可能不同（G×E互作）
   - 基因间互作（上位性）难以建模

3. **经济挑战**：
   - 初期投入大（建立训练群体）
   - 需要持续更新模型

#### 适用性评估
**GS效果好的情况**：
- ✅ 性状遗传力中等到高（h² > 0.3）
- ✅ 训练群体与选择群体亲缘关系近
- ✅ 标记密度足够高
- ✅ 训练群体足够大

**GS效果有限的情况**：
- ❌ 性状主要受环境影响（h² < 0.2）
- ❌ 性状受少数大效应基因控制
- ❌ 训练群体与选择群体差异很大
- ❌ 标记密度过低或训练群体过小

---

## 2. 数据准备和理解

### 2.1 GS数据的本质特征

#### 为什么需要这些特定的数据？

**【基因型数据的重要性】**
基因型数据是GS的"原料"，就像食谱需要食材一样：

1. **全基因组覆盖**：
   - **为什么需要？** 我们不知道哪些基因真正影响性状，所以需要"撒网捕鱼"
   - **类比理解**：就像在黑暗中找东西，需要打开所有灯光才能看清楚

2. **标记密度**：
   - **理想情况**：每个基因都有标记覆盖
   - **实际考虑**：平衡成本与效果，通常每几千个碱基有一个标记
   - **经验法则**：标记数量 > 个体数量时效果较好

3. **编码方式（0/1/2）**：
   - **0 = AA**：两个相同的参考等位基因
   - **1 = Aa**：一个参考，一个变异等位基因（杂合子）
   - **2 = aa**：两个变异等位基因
   - **为什么这样编码？** 反映等位基因剂量效应，便于加性模型计算

**【表型数据的关键要求】**

1. **准确性**：
   - **为什么重要？** "垃圾进，垃圾出" - 表型数据质量直接影响预测精度
   - **质控措施**：重复测量、标准化方案、异常值检测

2. **代表性**：
   - **环境代表性**：训练数据应覆盖目标环境
   - **群体代表性**：训练群体应包含选择群体的遗传多样性

3. **足够的样本量**：
   - **经验公式**：有效训练群体大小 ≥ (性状数量 × 标记数量) / 遗传力
   - **实际考虑**：通常需要几百到几千个个体

#### 数据质量的"黄金标准"

**【基因型数据质量指标】**
- **缺失率 < 10%**：过多缺失值影响预测精度
- **最小等位基因频率(MAF) > 0.05**：太稀有的变异提供的信息有限
- **Hardy-Weinberg平衡**：偏离可能提示基因分型错误

**【表型数据质量指标】**
- **重复性 > 0.7**：同一个体多次测量的一致性
- **遗传力 > 0.3**：遗传力太低的性状难以预测
- **正态分布**：偏态分布可能需要数据转换

### 2.2 我们的数据集详解

#### 数据来源与背景
我们使用的是**玉米多样性面板（Maize Diversity Panel, MDP）**数据：

1. **mdp_numeric.txt**：玉米多样性面板的基因型数据
   - **来源**：281个玉米自交系
   - **标记类型**：SNP（单核苷酸多态性）
   - **覆盖范围**：全基因组
   - **编码方式**：0/1/2数值编码

2. **mdp_SNP_information.txt**：SNP位置信息
   - **染色体号**：标记分布在哪个染色体
   - **物理位置**：在染色体上的具体位置（bp）
   - **用途**：用于绘制曼哈顿图、计算连锁不平衡

3. **mdp_traits.txt**：真实农艺性状表型数据
   - **EarHT**：穗位高（ear height，cm）
   - **dpoll**：开花期（days to pollen，天）
   - **EarDia**：穗径（ear diameter，mm）

4. **mdp_YRef.txt**：参考表型数据（用于本教程）
   - **SimTrait**：模拟性状，基于真实遗传结构
   - **用途**：教学演示，避免真实数据的复杂性

#### 为什么选择这个数据集？

1. **代表性强**：
   - 涵盖玉米的主要遗传群体
   - 包含重要农艺性状
   - 标记覆盖全基因组

2. **质量可靠**：
   - 经过严格质控
   - 被广泛用于研究
   - 结果可重现

3. **教学适宜**：
   - 数据量适中（不会太大导致计算困难）
   - 性状具有生物学意义
   - 有真实表型可供对比

### 2.3 数据格式深度解析

#### 基因型数据格式（mdp_numeric.txt）

```
Taxa        SNP1  SNP2  SNP3  ...
Individual1   0     1     2   ...
Individual2   1     2     0   ...
...
```

**为什么这样组织？**
- **行 = 个体**：方便按个体提取数据
- **列 = 标记**：方便按标记进行统计
- **数值编码**：便于数学运算和统计分析

**数值含义的生物学解释**：
```
基因型  |  编码  |  生物学含义
--------|--------|-------------
  AA    |   0    |  无变异等位基因
  Aa    |   1    |  一个变异等位基因
  aa    |   2    |  两个变异等位基因
```

#### SNP信息格式（mdp_SNP_information.txt）

```
SNP          Chromosome  Position
SNP1         1           157104
SNP2         1           158291
...
```

**用途解析**：
- **绘图需要**：曼哈顿图需要知道每个SNP的染色体和位置
- **质控需要**：检查标记分布是否均匀
- **生物学解释**：候选基因定位和功能注释

#### 表型数据格式对比

**真实农艺性状（mdp_traits.txt）**：
```
Taxa        EarHT   dpoll   EarDia
Individual1 59.5    NaN     NaN
Individual2 65.5    59.5    32.21933
...
```

**参考性状（mdp_YRef.txt）**：
```
Taxa	SimTrait
38-11	3.14387607839918
4226	2.2367231742118
...
```

**为什么有两种表型数据？**

| 特点 | 真实性状 | 参考性状 |
|------|----------|----------|
| **复杂性** | 高（多种环境因子影响） | 中（已知遗传结构） |
| **缺失值** | 多（测定困难） | 少（模拟生成） |
| **教学价值** | 接近实际应用 | 便于理解原理 |
| **分析难度** | 高 | 适中 |

### 2.4 数据需求和质量标准

#### GS分析的数据需求层次

**【最低要求】**（能跑通分析）
- 训练群体 ≥ 100个个体
- 标记数量 ≥ 1000个SNP
- 表型数据完整性 ≥ 80%
- 性状遗传力 ≥ 0.2

**【推荐配置】**（获得较好结果）
- 训练群体 ≥ 500个个体
- 标记数量 ≥ 10000个SNP
- 表型数据完整性 ≥ 95%
- 性状遗传力 ≥ 0.4

**【理想状态】**（最优预测效果）
- 训练群体 ≥ 2000个个体
- 标记数量 ≥ 50000个SNP
- 无缺失值，多环境重复
- 性状遗传力 ≥ 0.6

#### 数据质量评估清单

**【基因型数据】**
- [ ] **完整性检查**：缺失率是否可接受？
- [ ] **多态性检查**：是否存在单态SNP？
- [ ] **频率检查**：MAF是否在合理范围？
- [ ] **分布检查**：标记是否均匀分布在基因组？

**【表型数据】**
- [ ] **数值检查**：是否存在异常值？
- [ ] **分布检查**：是否近似正态分布？
- [ ] **重复性检查**：多次测量是否一致？
- [ ] **匹配检查**：个体名是否与基因型数据匹配？

**【关联性检查】**
- [ ] **样本匹配**：基因型与表型个体是否一致？
- [ ] **群体结构**：是否需要考虑亚群分化？
- [ ] **亲缘关系**：个体间是否存在密切亲缘关系？

这些检查看似繁琐，但每一步都关系到GS分析的成败！

---

## 3. 第一步：数据载入和预处理

### 3.1 设置工作环境 - 为什么这些步骤很重要？

#### 环境清理的必要性
```r
# 清理环境
rm(list = ls())
```

**为什么要清理环境？**
- **避免变量污染**：之前的分析可能留下同名变量，导致意外错误
- **内存释放**：R会保留所有对象在内存中，清理可释放内存空间
- **重现性保证**：确保每次运行都是从干净状态开始

**类比理解**：就像做实验前要清洁实验台一样，保证实验的纯净性！

#### 工作目录设置的重要性
```r
# 设置工作目录
setwd("/Users/your_path/statistical genomics")
```

**为什么要设置工作目录？**
- **文件定位**：R需要知道在哪里找数据文件
- **结果保存**：分析结果会保存在工作目录中
- **路径简化**：使用相对路径，代码更简洁

**检查方法**：
```r
# 查看当前工作目录
getwd()

# 查看目录中的文件
list.files()

# 检查数据文件是否存在
file.exists("data/mdp_numeric.txt")
```

#### 依赖包管理的策略
```r
# 载入必要的函数库
source("function/gapit_functions.R")

# 智能包管理
required_packages <- c("rrBLUP", "BGLR", "caret", "ggplot2", "corrplot")
for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
```

**为什么这样写代码？**

1. **自动化检查**：`if(!require(...))` 先尝试载入，失败才安装
2. **避免重复安装**：已安装的包不会重复安装，节省时间
3. **错误处理**：缺少包时自动安装，避免程序中断
4. **可移植性**：代码在不同机器上都能正常运行

**各包的作用**：
- **rrBLUP**: 岭回归BLUP和gBLUP
- **BGLR**: 贝叶斯回归方法
- **caret**: 交叉验证和模型评估
- **ggplot2**: 高质量图形绘制
- **corrplot**: 相关矩阵可视化

### 3.2 数据载入 - 读取策略与错误处理

#### 安全的数据读取方法
```r
# 带错误检查的数据载入
safe_read_table <- function(file_path, description) {
  if(!file.exists(file_path)) {
    stop(sprintf("错误：找不到%s文件：%s", description, file_path))
  }

  cat(sprintf("载入%s...", description))

  # 尝试读取数据
  tryCatch({
    data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    cat(" 成功!\n")
    return(data)
  }, error = function(e) {
    stop(sprintf("读取%s时出错：%s", description, e$message))
  })
}

# 使用安全读取函数
genotype_data <- safe_read_table("data/mdp_numeric.txt", "基因型数据")
snp_info <- safe_read_table("data/mdp_SNP_information.txt", "SNP信息")
trait_data <- safe_read_table("data/mdp_YRef.txt", "表型数据")
```

**代码解析**：

1. **`file.exists()`检查**：确保文件存在后再读取
2. **`tryCatch()`错误捕获**：优雅地处理读取错误
3. **`stringsAsFactors = FALSE`**：避免字符串被转换为因子
4. **分隔符指定**：明确指定制表符分隔

**为什么使用这种方法？**
- **减少调试时间**：明确的错误信息帮助快速定位问题
- **提高代码稳定性**：处理各种异常情况
- **便于维护**：统一的读取接口

#### 数据读取参数详解

```r
read.table(file_path,
          header = TRUE,        # 第一行是列名
          sep = "\t",          # 制表符分隔
          stringsAsFactors = FALSE,  # 字符串不转换为因子
          check.names = TRUE,  # 检查列名合法性
          na.strings = c("NA", "NaN", ""))  # 缺失值标识
```

**参数说明**：
- **header = TRUE**: 数据文件第一行包含列名
- **sep = "\\t"**: 使用制表符作为分隔符（适合从Excel导出的数据）
- **stringsAsFactors = FALSE**: 保持字符串为字符类型，避免意外的类型转换
- **na.strings**: 指定哪些字符串应被识别为缺失值

### 3.3 数据基本检查 - 质量评估的第一步

#### 维度检查：数据规模评估
```r
cat("=== 数据基本信息 ===\n")
cat(sprintf("个体数量: %d\n", nrow(genotype_data)))
cat(sprintf("SNP数量: %d\n", ncol(genotype_data) - 1))  # 减1是因为第一列是个体名
cat(sprintf("表型数据个体数: %d\n", nrow(trait_data)))
```

**为什么要检查维度？**
- **确认数据完整性**：是否有数据丢失或截断
- **评估分析复杂度**：大数据集需要更多内存和时间
- **验证数据格式**：行列数是否符合期望

#### 数据结构检查：理解数据组织
```r
# 检查数据结构
str(genotype_data)
str(trait_data)

# 检查数据类型
class(genotype_data)
sapply(genotype_data[, 1:5], class)
```

**`str()`函数的价值**：
- **数据类型**：每列的数据类型（数值型、字符型等）
- **数据维度**：行数和列数
- **缺失值概况**：NA的数量
- **数据预览**：前几个值的预览

#### 数据质量初步评估
```r
# 基因型数据质量检查
cat("\n=== 基因型数据质量 ===\n")

# 检查基因型编码是否正确（应该只有0, 1, 2）
genotype_matrix <- as.matrix(genotype_data[, -1])
unique_values <- unique(as.vector(genotype_matrix))
cat(sprintf("基因型编码: %s\n", paste(sort(unique_values), collapse = ", ")))

# 检查缺失值
missing_count <- sum(is.na(genotype_matrix))
total_count <- length(genotype_matrix)
cat(sprintf("缺失值: %d (%.2f%%)\n", missing_count, missing_count/total_count*100))

# 表型数据质量检查
cat("\n=== 表型数据质量 ===\n")
phenotype_col <- 2  # 假设表型在第2列
phenotype_values <- trait_data[, phenotype_col]

cat(sprintf("表型观测数: %d\n", sum(!is.na(phenotype_values))))
cat(sprintf("缺失值: %d\n", sum(is.na(phenotype_values))))
cat(sprintf("均值: %.3f\n", mean(phenotype_values, na.rm = TRUE)))
cat(sprintf("标准差: %.3f\n", sd(phenotype_values, na.rm = TRUE)))
cat(sprintf("范围: %.3f - %.3f\n",
            min(phenotype_values, na.rm = TRUE),
            max(phenotype_values, na.rm = TRUE)))
```

**质量检查的重要性**：
- **数据格式验证**：确保编码正确（0/1/2）
- **缺失值评估**：过多缺失值影响分析质量
- **异常值识别**：极端值可能是错误数据
- **分布特征**：正态性影响某些统计方法的有效性

#### 个体名匹配检查
```r
# 检查个体名匹配
genotype_individuals <- genotype_data[, 1]
trait_individuals <- trait_data[, 1]

cat("\n=== 个体匹配检查 ===\n")
cat(sprintf("基因型数据个体数: %d\n", length(genotype_individuals)))
cat(sprintf("表型数据个体数: %d\n", length(trait_individuals)))

# 找到共同个体
common_individuals <- intersect(genotype_individuals, trait_individuals)
cat(sprintf("共同个体数: %d\n", length(common_individuals)))

# 检查匹配率
match_rate <- length(common_individuals) / max(length(genotype_individuals),
                                               length(trait_individuals))
cat(sprintf("匹配率: %.2f%%\n", match_rate * 100))

# 列出不匹配的个体（前10个）
only_genotype <- setdiff(genotype_individuals, trait_individuals)
only_trait <- setdiff(trait_individuals, genotype_individuals)

if(length(only_genotype) > 0) {
  cat("只有基因型的个体（前10个）:\n")
  print(head(only_genotype, 10))
}

if(length(only_trait) > 0) {
  cat("只有表型的个体（前10个）:\n")
  print(head(only_trait, 10))
}
```

**为什么个体匹配如此重要？**
- **分析前提**：GS需要同一个体的基因型和表型数据
- **错误预防**：不匹配的个体名会导致分析失败
- **数据理解**：了解数据收集的完整性

#### 数据分布可视化检查
```r
# 表型分布检查
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# 直方图
hist(phenotype_values, breaks = 30, col = "lightblue",
     main = "表型分布", xlab = "表型值", ylab = "频数")

# 箱线图
boxplot(phenotype_values, col = "lightgreen",
        main = "表型分布", ylab = "表型值")

# QQ图（正态性检验）
qqnorm(phenotype_values, main = "正态性检验")
qqline(phenotype_values, col = "red")

# 密度图
plot(density(phenotype_values, na.rm = TRUE),
     main = "表型密度分布", col = "blue", lwd = 2)

# 正态性统计检验
shapiro_result <- shapiro.test(sample(phenotype_values[!is.na(phenotype_values)],
                                     min(5000, sum(!is.na(phenotype_values)))))
cat(sprintf("\n正态性检验 (Shapiro-Wilk): p-value = %.2e\n", shapiro_result$p.value))
```

**可视化检查的价值**：
- **分布形状**：正态、偏态、多峰分布
- **异常值识别**：箱线图中的离群点
- **数据转换需求**：严重偏态可能需要对数转换

### 3.4 高级数据检查

#### 基因型数据的遗传学检查
```r
# 等位基因频率检查
calculate_allele_freq <- function(genotypes) {
  # 计算变异等位基因频率
  freq <- mean(genotypes, na.rm = TRUE) / 2
  return(freq)
}

# 计算所有SNP的等位基因频率
allele_freqs <- apply(genotype_matrix, 2, calculate_allele_freq)

cat("=== 等位基因频率分布 ===\n")
cat(sprintf("频率范围: %.3f - %.3f\n",
            min(allele_freqs, na.rm = TRUE),
            max(allele_freqs, na.rm = TRUE)))
cat(sprintf("平均频率: %.3f\n", mean(allele_freqs, na.rm = TRUE)))

# MAF分布
maf <- pmin(allele_freqs, 1 - allele_freqs)
cat(sprintf("MAF < 0.05的SNP数量: %d (%.1f%%)\n",
            sum(maf < 0.05, na.rm = TRUE),
            sum(maf < 0.05, na.rm = TRUE) / length(maf) * 100))

# MAF分布直方图
hist(maf, breaks = 50, col = "orange",
     main = "最小等位基因频率(MAF)分布",
     xlab = "MAF", ylab = "SNP数量")
abline(v = 0.05, col = "red", lty = 2, lwd = 2)
text(0.1, max(hist(maf, plot = FALSE)$counts) * 0.8,
     "MAF = 0.05", col = "red")
```

**等位基因频率检查的意义**：
- **信息量评估**：极低频率的SNP提供的信息有限
- **统计功效**：MAF太低影响关联检测的统计功效
- **群体代表性**：频率分布反映群体的遗传多样性

#### 染色体分布检查
```r
# 检查SNP在染色体上的分布
if(ncol(snp_info) >= 2) {
  chr_distribution <- table(snp_info[, 2])  # 假设第2列是染色体号

  cat("\n=== 染色体分布 ===\n")
  print(chr_distribution)

  # 绘制染色体分布图
  barplot(chr_distribution,
          col = rainbow(length(chr_distribution)),
          main = "SNP在染色体上的分布",
          xlab = "染色体", ylab = "SNP数量")

  # 检查分布均匀性
  chi_test <- chisq.test(chr_distribution)
  cat(sprintf("分布均匀性检验: p-value = %.2e\n", chi_test$p.value))
}
```

**染色体分布检查的作用**：
- **覆盖度评估**：确保全基因组都有标记覆盖
- **密度均匀性**：过度稀疏的区域可能遗漏重要基因
- **质控指标**：严重不均匀可能提示技术问题

通过这些详细的检查步骤，我们可以：
1. **提前发现问题**：避免在后续分析中遇到错误
2. **了解数据特性**：为选择合适的分析方法提供依据
3. **建立信心**：确认数据质量后可放心进行分析
4. **记录基线信息**：为结果解释提供背景

---

## 4. 第二步：训练测试集划分

### 4.1 数据匹配和整理

```r
# 匹配基因型和表型数据
genotype_individuals <- genotype_data[, 1]
trait_individuals <- trait_data[, 1]
common_individuals <- intersect(genotype_individuals, trait_individuals)

cat(sprintf("共同个体数: %d\n", length(common_individuals)))

# 根据共同个体筛选数据
genotype_matched_idx <- match(common_individuals, genotype_individuals)
trait_matched_idx <- match(common_individuals, trait_individuals)

matched_genotype_data <- genotype_data[genotype_matched_idx, ]
matched_trait_data <- trait_data[trait_matched_idx, ]

# 验证匹配
if(!all(matched_genotype_data[, 1] == matched_trait_data[, 1])) {
  stop("个体名匹配失败！")
}
cat("个体名匹配成功！\n")
```

### 4.2 划分训练和测试集

```r
# 设置随机种子确保结果可重现
set.seed(99164)

# 划分策略：80%训练，20%测试
n_total <- nrow(matched_genotype_data)
test_size <- round(n_total / 5)  # 20%用于测试
testing_idx <- sample(n_total, test_size, replace = FALSE)
training_idx <- setdiff(1:n_total, testing_idx)

cat(sprintf("总个体数: %d\n", n_total))
cat(sprintf("训练集个体数: %d\n", length(training_idx)))
cat(sprintf("测试集个体数: %d\n", length(testing_idx)))

# 提取训练和测试数据
# 基因型矩阵
X_all <- as.matrix(matched_genotype_data[, -1])
rownames(X_all) <- matched_genotype_data[, 1]

X_train <- X_all[training_idx, ]
X_test <- X_all[testing_idx, ]

# 表型向量
y_all <- matched_trait_data[, 2]
names(y_all) <- matched_trait_data[, 1]

y_train <- y_all[training_idx]
y_test <- y_all[testing_idx]

# 个体名称
taxa_train <- matched_genotype_data[training_idx, 1]
taxa_test <- matched_genotype_data[testing_idx, 1]
```

---

## 5. 第三步：线性混合模型方法

### 5.1 Ridge Regression BLUP (RR-BLUP)

```r
cat("=== Ridge Regression BLUP (RR-BLUP) ===\n")

library(rrBLUP)

# 方法1：直接使用标记效应模型
# 模型：y = Xβ + Zg + e，其中 g ~ N(0, σ²gI)
rr_blup_result <- mixed.solve(y = y_train, Z = X_train)

# 预测测试集
y_pred_rrblup <- X_test %*% rr_blup_result$u

# 计算预测精度
accuracy_rrblup <- cor(y_test, y_pred_rrblup)
cat(sprintf("RR-BLUP预测精度: %.3f\n", accuracy_rrblup))
cat(sprintf("RR-BLUP预测准确性 (R²): %.3f\n", accuracy_rrblup^2))
```

### 5.2 基因组BLUP (gBLUP)

```r
cat("\n=== 基因组BLUP (gBLUP) ===\n")

# 方法2：使用基因组关系矩阵
# 计算基因组关系矩阵 K = XX'/p，其中p为标记数
K_train <- tcrossprod(X_train) / ncol(X_train)

# gBLUP模型：y = 1μ + Kg + e，其中 g ~ N(0, σ²gK)
gblup_result <- mixed.solve(y = y_train, K = K_train)

# 预测测试集需要计算训练集与测试集间的关系矩阵
K_test_train <- tcrossprod(X_test, X_train) / ncol(X_train)
y_pred_gblup <- K_test_train %*% gblup_result$u

# 计算预测精度
accuracy_gblup <- cor(y_test, y_pred_gblup)
cat(sprintf("gBLUP预测精度: %.3f\n", accuracy_gblup))
cat(sprintf("gBLUP预测准确性 (R²): %.3f\n", accuracy_gblup^2))

# 验证两种方法的等价性
gebv_rrblup <- X_all %*% rr_blup_result$u
gebv_gblup_all <- rbind(gblup_result$u, y_pred_gblup)
equivalence_cor <- cor(gebv_rrblup[training_idx], gblup_result$u)
cat(sprintf("RR-BLUP与gBLUP的等价性验证: %.6f\n", equivalence_cor))
```

### 5.3 使用GAPIT进行gBLUP

```r
cat("\n=== 使用GAPIT进行gBLUP ===\n")

# 准备GAPIT所需的数据格式
train_GD <- data.frame(Taxa = taxa_train, X_train)
train_Y <- data.frame(Taxa = taxa_train, Trait = y_train)

# 使用GAPIT进行gBLUP分析
gapit_gblup <- GAPIT(
  Y = train_Y,
  GD = train_GD,
  GM = snp_info,
  model = "gBLUP",
  SNP.test = FALSE,  # 不进行SNP检验，只做预测
  PCA.total = 3      # 使用3个主成分控制群体结构
)

# 提取预测结果
pred_results <- gapit_gblup$Pred
colnames(pred_results)[ncol(pred_results)] <- "GEBV"

# 匹配测试集个体的预测值
test_pred_idx <- match(taxa_test, pred_results$Taxa)
y_pred_gapit <- pred_results$GEBV[test_pred_idx]

# 计算预测精度
accuracy_gapit <- cor(y_test, y_pred_gapit)
cat(sprintf("GAPIT gBLUP预测精度: %.3f\n", accuracy_gapit))
cat(sprintf("GAPIT gBLUP预测准确性 (R²): %.3f\n", accuracy_gapit^2))
```

---

## 6. 第四步：贝叶斯回归方法

### 6.1 贝叶斯Ridge回归 (BRR)

```r
cat("=== 贝叶斯Ridge回归 (BRR) ===\n")

library(BGLR)

# 设置MCMC参数
nIter <- 2000   # 迭代次数
burnIn <- 500   # 老化期

# 准备协变量矩阵（主成分）
pca_result <- prcomp(X_train, center = TRUE, scale. = FALSE)
pc_scores <- pca_result$x[, 1:3]

# 贝叶斯Ridge回归
set.seed(99164)
brr_model <- BGLR(
  y = y_train,
  ETA = list(
    list(X = pc_scores, model = 'FIXED'),
    list(X = X_train, model = 'BRR')
  ),
  nIter = nIter,
  burnIn = burnIn,
  verbose = FALSE
)

# 预测测试集
pc_scores_test <- predict(pca_result, X_test)[, 1:3]
y_pred_brr <- pc_scores_test %*% brr_model$ETA[[1]]$b +
              X_test %*% brr_model$ETA[[2]]$b

# 计算预测精度
accuracy_brr <- cor(y_test, y_pred_brr)
cat(sprintf("BRR预测精度: %.3f\n", accuracy_brr))
cat(sprintf("BRR预测准确性 (R²): %.3f\n", accuracy_brr^2))
```

### 6.2 贝叶斯LASSO (BL)

```r
cat("\n=== 贝叶斯LASSO (BL) ===\n")

set.seed(99164)
bl_model <- BGLR(
  y = y_train,
  ETA = list(
    list(X = pc_scores, model = 'FIXED'),
    list(X = X_train, model = 'BL')
  ),
  nIter = nIter,
  burnIn = burnIn,
  verbose = FALSE
)

# 预测测试集
y_pred_bl <- pc_scores_test %*% bl_model$ETA[[1]]$b +
             X_test %*% bl_model$ETA[[2]]$b

# 计算预测精度
accuracy_bl <- cor(y_test, y_pred_bl)
cat(sprintf("BL预测精度: %.3f\n", accuracy_bl))
cat(sprintf("BL预测准确性 (R²): %.3f\n", accuracy_bl^2))
```

### 6.3 BayesA方法

```r
cat("\n=== BayesA方法 ===\n")

set.seed(99164)
bayesa_model <- BGLR(
  y = y_train,
  ETA = list(
    list(X = pc_scores, model = 'FIXED'),
    list(X = X_train, model = 'BayesA')
  ),
  nIter = nIter,
  burnIn = burnIn,
  verbose = FALSE
)

# 预测测试集
y_pred_bayesa <- pc_scores_test %*% bayesa_model$ETA[[1]]$b +
                 X_test %*% bayesa_model$ETA[[2]]$b

# 计算预测精度
accuracy_bayesa <- cor(y_test, y_pred_bayesa)
cat(sprintf("BayesA预测精度: %.3f\n", accuracy_bayesa))
cat(sprintf("BayesA预测准确性 (R²): %.3f\n", accuracy_bayesa^2))
```

---

## 7. 第五步：神经网络方法

### 7.1 数据预处理

```r
cat("=== 神经网络方法 ===\n")

# 注意：神经网络方法需要安装keras和tensorflow
# install.packages(c("keras", "tensorflow"))
# library(keras)
# install_tensorflow()

# 如果无法安装keras，可以跳过此部分，使用其他方法

# 数据标准化
X_train_scaled <- scale(X_train)
X_test_scaled <- scale(X_test,
                       center = attr(X_train_scaled, "scaled:center"),
                       scale = attr(X_train_scaled, "scaled:scale"))

y_train_scaled <- scale(y_train)
y_test_scaled <- scale(y_test,
                       center = attr(y_train_scaled, "scaled:center"),
                       scale = attr(y_train_scaled, "scaled:scale"))

# 处理缺失值
X_train_scaled[is.na(X_train_scaled)] <- 0
X_test_scaled[is.na(X_test_scaled)] <- 0
```

### 7.2 简单的线性神经网络

```r
# 由于神经网络需要复杂的环境配置，这里提供一个简化的实现思路
# 实际应用中建议参考NN-GS.R文件中的完整实现

# 替代方案：使用岭回归模拟神经网络的线性部分
library(glmnet)

# 使用交叉验证选择最优lambda
cv_ridge <- cv.glmnet(X_train_scaled, y_train_scaled, alpha = 0)
optimal_lambda <- cv_ridge$lambda.min

# 训练岭回归模型
ridge_model <- glmnet(X_train_scaled, y_train_scaled,
                      alpha = 0, lambda = optimal_lambda)

# 预测
y_pred_ridge <- predict(ridge_model, newx = X_test_scaled)

# 反标准化
y_pred_ridge_original <- y_pred_ridge * attr(y_train_scaled, "scaled:scale") +
                        attr(y_train_scaled, "scaled:center")

# 计算预测精度
accuracy_ridge <- cor(y_test, y_pred_ridge_original)
cat(sprintf("Ridge回归预测精度: %.3f\n", accuracy_ridge))
cat(sprintf("Ridge回归预测准确性 (R²): %.3f\n", accuracy_ridge^2))
```

---

## 8. 第六步：方法比较和评估

### 8.1 汇总预测结果

```r
cat("\n=== 方法性能比较 ===\n")

# 汇总所有方法的预测结果
methods <- c("RR-BLUP", "gBLUP", "GAPIT_gBLUP", "BRR", "BL", "BayesA", "Ridge")
predictions <- list(y_pred_rrblup, y_pred_gblup, y_pred_gapit,
                   y_pred_brr, y_pred_bl, y_pred_bayesa, y_pred_ridge_original)
accuracies <- c(accuracy_rrblup, accuracy_gblup, accuracy_gapit,
               accuracy_brr, accuracy_bl, accuracy_bayesa, accuracy_ridge)

# 创建结果汇总表
results_summary <- data.frame(
  Method = methods,
  Accuracy = round(accuracies, 3),
  R_squared = round(accuracies^2, 3)
)

print(results_summary)

# 找出最佳方法
best_method_idx <- which.max(accuracies)
cat(sprintf("\n最佳方法: %s (精度: %.3f)\n",
            methods[best_method_idx], accuracies[best_method_idx]))
```

### 8.2 可视化比较

```r
# 创建预测精度比较图
library(ggplot2)

# 准备绘图数据
plot_data <- data.frame(
  Method = factor(methods, levels = methods[order(accuracies, decreasing = TRUE)]),
  Accuracy = accuracies,
  R_squared = accuracies^2
)

# 预测精度柱状图
p1 <- ggplot(plot_data, aes(x = Method, y = Accuracy, fill = Method)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  geom_text(aes(label = round(Accuracy, 3)), vjust = -0.3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "不同GS方法的预测精度比较",
       y = "预测精度 (相关系数)",
       x = "方法") +
  guides(fill = FALSE)

print(p1)

# 最佳方法的散点图
best_pred <- predictions[[best_method_idx]]
p2 <- ggplot(data.frame(Observed = y_test, Predicted = as.vector(best_pred)),
            aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.6, color = "blue") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(title = paste("最佳方法 (", methods[best_method_idx], ") 预测效果"),
       x = "观测值",
       y = "预测值") +
  annotate("text", x = min(y_test), y = max(as.vector(best_pred)),
           label = sprintf("r = %.3f", accuracies[best_method_idx]),
           hjust = 0, vjust = 1) +
  theme_minimal()

print(p2)
```

---

## 9. 第七步：交叉验证策略

### 9.1 K折交叉验证

```r
cat("=== 5折交叉验证评估 ===\n")

library(caret)

# 设置5折交叉验证
set.seed(123)
k_folds <- 5
folds <- createFolds(y_all, k = k_folds, returnTrain = TRUE)

# 存储交叉验证结果
cv_results <- data.frame(
  Fold = integer(),
  Method = character(),
  Accuracy = numeric(),
  stringsAsFactors = FALSE
)

# 只对主要方法进行交叉验证（RR-BLUP, gBLUP, BRR）
cv_methods <- c("RR-BLUP", "gBLUP", "BRR")

for(fold in 1:k_folds) {
  cat(sprintf("处理第%d折...\n", fold))

  # 获取训练和验证索引
  train_idx <- folds[[fold]]
  val_idx <- setdiff(1:length(y_all), train_idx)

  # 提取训练和验证数据
  X_fold_train <- X_all[train_idx, ]
  y_fold_train <- y_all[train_idx]
  X_fold_val <- X_all[val_idx, ]
  y_fold_val <- y_all[val_idx]

  # RR-BLUP
  rrblup_fold <- mixed.solve(y = y_fold_train, Z = X_fold_train)
  pred_rrblup_fold <- X_fold_val %*% rrblup_fold$u
  acc_rrblup_fold <- cor(y_fold_val, pred_rrblup_fold)

  cv_results <- rbind(cv_results, data.frame(
    Fold = fold, Method = "RR-BLUP", Accuracy = acc_rrblup_fold
  ))

  # gBLUP
  K_fold_train <- tcrossprod(X_fold_train) / ncol(X_fold_train)
  gblup_fold <- mixed.solve(y = y_fold_train, K = K_fold_train)
  K_fold_val <- tcrossprod(X_fold_val, X_fold_train) / ncol(X_fold_train)
  pred_gblup_fold <- K_fold_val %*% gblup_fold$u
  acc_gblup_fold <- cor(y_fold_val, pred_gblup_fold)

  cv_results <- rbind(cv_results, data.frame(
    Fold = fold, Method = "gBLUP", Accuracy = acc_gblup_fold
  ))

  # BRR（简化版）
  pca_fold <- prcomp(X_fold_train, center = TRUE, scale. = FALSE)
  pc_fold_train <- pca_fold$x[, 1:3]
  pc_fold_val <- predict(pca_fold, X_fold_val)[, 1:3]

  brr_fold <- BGLR(
    y = y_fold_train,
    ETA = list(
      list(X = pc_fold_train, model = 'FIXED'),
      list(X = X_fold_train, model = 'BRR')
    ),
    nIter = 1000, burnIn = 200, verbose = FALSE
  )

  pred_brr_fold <- pc_fold_val %*% brr_fold$ETA[[1]]$b +
                   X_fold_val %*% brr_fold$ETA[[2]]$b
  acc_brr_fold <- cor(y_fold_val, pred_brr_fold)

  cv_results <- rbind(cv_results, data.frame(
    Fold = fold, Method = "BRR", Accuracy = acc_brr_fold
  ))
}

# 汇总交叉验证结果
cv_summary <- aggregate(Accuracy ~ Method, data = cv_results,
                       FUN = function(x) c(Mean = mean(x), SD = sd(x)))
cv_summary <- do.call(data.frame, cv_summary)
colnames(cv_summary) <- c("Method", "Mean_Accuracy", "SD_Accuracy")

print(cv_summary)

# 可视化交叉验证结果
library(ggplot2)
p3 <- ggplot(cv_results, aes(x = Method, y = Accuracy, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "red") +
  theme_minimal() +
  labs(title = "5折交叉验证结果",
       y = "预测精度",
       x = "方法") +
  guides(fill = FALSE)

print(p3)
```

---

## 10. 完整代码示例

```r
# ======================================================================
# 完整的GS分析流程脚本
# ======================================================================

# 环境设置
rm(list = ls())
setwd("/Users/your_path/statistical genomics")

# 载入必要库
source("function/gapit_functions.R")
library(rrBLUP)
library(BGLR)
library(ggplot2)
library(corrplot)

# 数据载入
genotype_data <- read.table("data/mdp_numeric.txt", header = TRUE, sep = "\t")
snp_info <- read.table("data/mdp_SNP_information.txt", header = TRUE, sep = "\t")
trait_data <- read.table("data/mdp_YRef.txt", header = TRUE, sep = "\t")

# 数据匹配和预处理
genotype_individuals <- genotype_data[, 1]
trait_individuals <- trait_data[, 1]
common_individuals <- intersect(genotype_individuals, trait_individuals)

genotype_matched_idx <- match(common_individuals, genotype_individuals)
trait_matched_idx <- match(common_individuals, trait_individuals)

matched_genotype_data <- genotype_data[genotype_matched_idx, ]
matched_trait_data <- trait_data[trait_matched_idx, ]

# 划分训练测试集
set.seed(99164)
n_total <- nrow(matched_genotype_data)
test_size <- round(n_total / 5)
testing_idx <- sample(n_total, test_size, replace = FALSE)
training_idx <- setdiff(1:n_total, testing_idx)

# 准备数据矩阵
X_all <- as.matrix(matched_genotype_data[, -1])
y_all <- matched_trait_data[, 2]
rownames(X_all) <- matched_genotype_data[, 1]
names(y_all) <- matched_trait_data[, 1]

X_train <- X_all[training_idx, ]
X_test <- X_all[testing_idx, ]
y_train <- y_all[training_idx]
y_test <- y_all[testing_idx]

# GS方法比较
methods <- character()
accuracies <- numeric()

# 1. RR-BLUP
rrblup_model <- mixed.solve(y = y_train, Z = X_train)
pred_rrblup <- X_test %*% rrblup_model$u
acc_rrblup <- cor(y_test, pred_rrblup)
methods <- c(methods, "RR-BLUP")
accuracies <- c(accuracies, acc_rrblup)

# 2. gBLUP
K_train <- tcrossprod(X_train) / ncol(X_train)
gblup_model <- mixed.solve(y = y_train, K = K_train)
K_test_train <- tcrossprod(X_test, X_train) / ncol(X_train)
pred_gblup <- K_test_train %*% gblup_model$u
acc_gblup <- cor(y_test, pred_gblup)
methods <- c(methods, "gBLUP")
accuracies <- c(accuracies, acc_gblup)

# 3. 贝叶斯Ridge回归
pca_result <- prcomp(X_train, center = TRUE, scale. = FALSE)
pc_scores <- pca_result$x[, 1:3]
pc_scores_test <- predict(pca_result, X_test)[, 1:3]

set.seed(99164)
brr_model <- BGLR(
  y = y_train,
  ETA = list(
    list(X = pc_scores, model = 'FIXED'),
    list(X = X_train, model = 'BRR')
  ),
  nIter = 1500, burnIn = 300, verbose = FALSE
)

pred_brr <- pc_scores_test %*% brr_model$ETA[[1]]$b +
            X_test %*% brr_model$ETA[[2]]$b
acc_brr <- cor(y_test, pred_brr)
methods <- c(methods, "BRR")
accuracies <- c(accuracies, acc_brr)

# 结果汇总
results <- data.frame(
  Method = methods,
  Accuracy = round(accuracies, 3),
  R_squared = round(accuracies^2, 3)
)

cat("=== GS方法性能比较 ===\n")
print(results)

# 最佳方法
best_idx <- which.max(accuracies)
cat(sprintf("\n最佳方法: %s (精度: %.3f)\n",
            methods[best_idx], accuracies[best_idx]))

cat("\n=== GS分析完成 ===\n")
cat("主要发现:\n")
cat("- 比较了多种基因组选择方法\n")
cat("- 评估了不同方法的预测精度\n")
cat("- 提供了完整的GS分析流程\n")
cat("- 为实际育种应用提供了参考\n")

cat("\n在实际应用中，还需要考虑:\n")
cat("- 训练群体大小和遗传多样性\n")
cat("- 标记密度和基因组覆盖度\n")
cat("- 性状遗传力和遗传结构\n")
cat("- 计算效率和实施成本\n")
```

---

## 总结

这个教程展示了全基因组选择的完整分析流程：

### 核心概念
1. **GS原理**：利用全基因组标记预测育种值
2. **方法多样性**：线性混合模型、贝叶斯方法、机器学习
3. **预测精度**：评估不同方法的预测性能
4. **实际应用**：为育种决策提供科学依据

### 方法特点
- **RR-BLUP/gBLUP**：计算效率高，适合大群体
- **贝叶斯方法**：灵活建模，考虑先验信息
- **神经网络**：捕捉复杂非线性关系
- **交叉验证**：客观评估预测性能

### 实践要点
1. 训练群体与目标群体的相关性是关键
2. 标记密度需要平衡成本与精度
3. 性状遗传结构影响方法选择
4. 需要持续更新训练数据

这个框架为基因组选择实践提供了坚实基础，可以根据具体育种目标进行调整和优化。