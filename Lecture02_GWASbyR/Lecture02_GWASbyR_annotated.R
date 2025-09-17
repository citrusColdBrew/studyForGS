# =============================================================================
# Lecture02: 全基因组关联分析 (GWAS) 入门教程
# Genome-Wide Association Study with R
# =============================================================================
#
# 【课程目标】
# 1. 理解GWAS的生物学原理和统计学基础
# 2. 学习基本的关联分析方法：相关性分析和线性回归
# 3. 了解群体结构对GWAS结果的影响
# 4. 掌握亲缘关系矩阵（Kinship Matrix）的计算和应用
#
# 【生物学背景知识】
# GWAS是一种在全基因组范围内寻找与复杂性状（如身高、疾病易感性）
# 相关的遗传变异的方法。基本原理是：
# 1. 基因型-表型关联：不同基因型的个体在表型上存在差异
# 2. 连锁不平衡：致病变异附近的标记与致病变异存在关联
# 3. 统计关联：通过统计方法检测基因型与表型的关联强度
# =============================================================================

# =============================================================================
# 第一部分：基于相关性的简单GWAS方法
# =============================================================================

# 【方法原理】
# 相关性分析是最简单的关联分析方法，直接计算每个SNP与表型的相关系数
# 优点：计算简单，易于理解
# 缺点：无法控制混杂因素（如群体结构、亲缘关系）

# 定义基于相关性的GWAS函数
# 【统计学原理】
# 1. 计算Pearson相关系数：r = cov(X,Y) / (sd(X) * sd(Y))
# 2. 将相关系数转换为t统计量：t = r * sqrt((n-2)/(1-r²))
# 3. 计算双尾p值：p = 2 * P(T > |t|)，其中T服从自由度为n-2的t分布
GWASbyCor <- function(X, y) {
     # 【输入参数说明】
     # X: 基因型矩阵 (行=个体, 列=SNP, 数值=0/1/2表示等位基因数量)
     # y: 表型向量 (连续型性状值)

     # 【生物学意义】
     # 基因型编码：0=AA(参考纯合子), 1=Aa(杂合子), 2=aa(变异纯合子)
     # 表型：可测量的生物学性状，如身高、产量等

     n <- nrow(X) # 样本个体数量
     r <- cor(y, X) # 计算表型与每个SNP的相关系数
     n <- nrow(X) # 重复定义（代码冗余）
     t <- r / sqrt((1 - r^2) / (n - 2)) # 相关系数的t统计量
     p <- 2 * (1 - pt(abs(t), n - 2)) # 双尾t检验的p值

     # 【统计学处理】
     # 避免p值为0的情况（在实际计算中可能出现）
     zeros <- p == 0 # 找出p值为0的位置
     p[zeros] <- 1e-10 # 将0替换为极小值
     return(p) # 返回每个SNP的p值向量
}

# =============================================================================
# 第二部分：数据载入和预处理
# =============================================================================

# 【数据来源说明】
# 使用的是玉米多样性面板（Maize Diversity Panel）数据集
# - 包含281个玉米自交系的基因型数据
# - 每个个体有3093个SNP标记
# - 提供了表型数据和SNP位置信息

# 设置工作目录（需要根据实际路径调整）
setwd("/Users/Jiabo/Documents/China/SWUN/Conference/2025西安培训/lecture/Code_Materials")

# 载入GAPIT功能函数
source("/Users/sameen/workspace/statistical genomics/function/gapit_functions.R")

# 【数据读取】
# 读取基因型数据：第一列为品系名称，其余列为SNP基因型(0/1/2编码)
myGD <- read.table("mdp_numeric.txt", head = T)
# 读取SNP信息：包含SNP名称、染色体号、物理位置
myGM <- read.table("mdp_SNP_information.txt", head = T)

# 载入自定义函数
source("G2P.R") # 基因型到表型模拟函数
source("GWASbyCor.R") # 相关性GWAS函数

# =============================================================================
# 第三部分：表型模拟和简单关联分析
# =============================================================================

# 【数据预处理】
X <- myGD[, -1] # 提取基因型矩阵（去除第一列品系名）
index1to5 <- myGM[, 2] < 6 # 选择1-5号染色体的SNP索引
X1to5 <- X[, index1to5] # 提取1-5号染色体的基因型数据

# 【表型模拟】
# 【生物学意义】在真实GWAS中，我们通常有观测的表型数据
# 这里为了教学目的，我们模拟一个已知QTN（数量性状位点）的表型
set.seed(99164) # 设置随机种子，确保结果可重现
mySim <- G2P(X = X1to5, h2 = .75, alpha = 1, NQTN = 10, distribution = "norm")
# 参数说明：
# h2 = 0.75: 遗传力，表示遗传因素解释75%的表型变异
# NQTN = 10: 模拟10个QTN
# distribution = "norm": QTN效应服从正态分布

# 【第一次GWAS分析：相关性方法】
p <- GWASbyCor(X = X, y = mySim$y) # 使用全部SNP进行关联分析

# 【结果可视化：曼哈顿图】
# 【图形解释】曼哈顿图是GWAS结果的标准展示方式
# x轴：SNP在基因组上的位置
# y轴：-log10(p值)，值越大表示关联越显著
color.vector <- rep(c("deepskyblue", "orange", "forestgreen", "indianred3"), 10) # 不同染色体用不同颜色
m <- nrow(myGM) # SNP总数
plot(t(-log10(p)) ~ seq(1:m),
     col = color.vector[myGM[, 2]],
     main = "Manhattan Plot - Correlation Method",
     xlab = "SNP Index", ylab = "-log10(P-value)"
)
# 添加真实QTN位置的垂直线（已知答案）
abline(v = mySim$QTN.position, lty = 2, lwd = 2, col = "black")

# =============================================================================
# 第四部分：QQ图诊断和假阳性控制
# =============================================================================

# 【QQ图原理】
# QQ图比较观测到的p值分布与期望的均匀分布
# 如果没有关联信号且无系统偏差，点应该落在y=x直线上
# 偏离直线表示：1) 真实关联信号 2) 系统偏差（如群体结构）

# 提取测试集数据（6-10号染色体，不包含模拟的QTN）
p.obs <- p[!index1to5] # 观测的p值
m2 <- length(p.obs) # 测试SNP数量
p.uni <- runif(m2, 0, 1) # 生成均匀分布的期望p值
order.obs <- order(p.obs) # 观测p值的排序索引
order.uni <- order(p.uni) # 期望p值的排序索引

# 绘制QQ图
plot(-log10(p.uni[order.uni]), -log10(p.obs[order.obs]),
     main = "QQ Plot - Testing for Inflation",
     xlab = "Expected -log10(P)", ylab = "Observed -log10(P)"
)
abline(a = 0, b = 1, col = "red") # 添加y=x参考线

# 【解释】如果点严重偏离红线上方，说明存在假阳性问题

# =============================================================================
# 第五部分：单SNP效应分析
# =============================================================================

# 【目的】分析最显著SNP的基因型-表型关系

# 找到最显著的SNP
order.obs <- order(p.obs) # 按p值排序
X6to10 <- X[, !index1to5] # 6-10号染色体的基因型
Xtop <- X6to10[, order.obs[1]] # 最显著SNP的基因型

# 【箱线图分析】
# 展示不同基因型组的表型分布
boxplot(mySim$y ~ Xtop,
     main = "Phenotype by Genotype",
     xlab = "Genotype (0=AA, 1=Aa, 2=aa)",
     ylab = "Phenotype Value"
)

# =============================================================================
# 第六部分：主成分分析（PCA）和群体结构控制
# =============================================================================

# 【生物学背景】
# 群体结构是GWAS中的主要混杂因素
# 不同亚群的个体在等位基因频率和表型均值上都可能存在差异
# 这会导致虚假关联（假阳性）

# 【PCA原理】
# 主成分分析可以捕捉基因型数据中的主要变异模式
# 前几个主成分通常反映群体结构和亚群分化

# 对基因型矩阵进行PCA
PCA <- prcomp(X) # 计算主成分

# 【群体结构可视化】
# 查看表型与第二主成分的关系
plot(mySim$y, PCA$x[, 2],
     main = "Phenotype vs PC2",
     xlab = "Phenotype", ylab = "PC2 Score"
)
cor(mySim$y, PCA$x[, 2]) # 计算相关性

# 【小样本演示】
# 随机选择10个个体，观察在小样本下的相关性
set.seed(99164)
s <- sample(length(mySim$y), 10) # 随机抽样10个个体
plot(mySim$y[s], PCA$x[s, 2],
     main = "Small Sample: Phenotype vs PC2",
     xlab = "Phenotype", ylab = "PC2 Score"
)
cor(mySim$y[s], PCA$x[s, 2]) # 小样本相关性

# =============================================================================
# 第七部分：广义线性模型（GLM）方法
# =============================================================================

# 【统计学原理】
# GLM方法通过多元线性回归分析基因型-表型关联
# 模型：Y = β₀ + β₁×PC2 + β₂×SNP + ε
# 其中PC2控制群体结构，SNP是待检验的遗传标记

# 【手动实现GLM分析】
y <- mySim$y # 表型向量
X <- cbind(1, PCA$x[, 2], Xtop) # 设计矩阵：截距项+PC2+最显著SNP

# 【矩阵运算】
# 正规方程：β = (X'X)⁻¹X'y
LHS <- t(X) %*% X # 左侧矩阵 X'X
C <- solve(LHS) # 逆矩阵 (X'X)⁻¹
RHS <- t(X) %*% y # 右侧向量 X'y
b <- C %*% RHS # 回归系数估计

# 【统计推断】
yb <- X %*% b # 拟合值
e <- y - yb # 残差
n <- length(y) # 样本量
ve <- sum(e^2) / (n - 1) # 残差方差估计
vt <- C * ve # 回归系数的方差-协方差矩阵
t <- b / sqrt(diag(vt)) # t统计量
p <- 2 * (1 - pt(abs(t), n - 2)) # p值

# 【结果整理】
LM <- cbind(b, t, sqrt(diag(vt)), p) # 合并结果
rownames(LM) <- cbind("Mean", "PC2", "Xtop") # 行名
colnames(LM) <- cbind("b", "t", "SD", "p") # 列名
print("GLM Analysis Results:")
print(LM)

# =============================================================================
# 第八部分：全基因组GLM分析（控制一个主成分）
# =============================================================================

# 【批量分析】
# 对每个SNP进行GLM分析，控制第二主成分
G <- myGD[, -1] # 基因型矩阵
n <- nrow(G) # 个体数
m <- ncol(G) # SNP数
P <- matrix(NA, 1, m) # 初始化p值矩阵

# 【循环分析每个SNP】
for (i in 1:m) {
     x <- G[, i] # 当前SNP的基因型

     # 【质量控制】检查SNP是否为单态（所有个体基因型相同）
     if (max(x) == min(x)) {
          p <- 1 # 单态SNP设p值为1
     } else {
          # 【GLM分析】
          X <- cbind(1, PCA$x[, 2], x) # 设计矩阵：截距+PC2+SNP
          LHS <- t(X) %*% X # X'X
          C <- solve(LHS) # (X'X)⁻¹
          RHS <- t(X) %*% y # X'y
          b <- C %*% RHS # 回归系数
          yb <- X %*% b # 拟合值
          e <- y - yb # 残差
          n <- length(y) # 样本量
          ve <- sum(e^2) / (n - 1) # 残差方差
          vt <- C * ve # 方差-协方差矩阵
          t <- b / sqrt(diag(vt)) # t统计量
          p <- 2 * (1 - pt(abs(t), n - 2)) # p值计算
     }
     P[i] <- p[length(p)] # 存储SNP的p值（最后一个是SNP效应的p值）
}

# 【结果评估：QQ图】
p.obs <- P[!index1to5] # 测试集p值
m2 <- length(p.obs) # 测试SNP数量
p.uni <- runif(m2, 0, 1) # 期望p值
order.obs <- order(p.obs) # 排序
order.uni <- order(p.uni)

plot(-log10(p.uni[order.uni]), -log10(p.obs[order.obs]),
     ylim = c(0, 7),
     main = "QQ Plot - GLM with 1 PC",
     xlab = "Expected -log10(P)", ylab = "Observed -log10(P)"
)
abline(a = 0, b = 1, col = "red")

# =============================================================================
# 第九部分：改进的GLM分析（控制三个主成分）
# =============================================================================

# 【分析改进】
# 增加更多主成分可以更好地控制群体结构
# 但也会降低检验功效（减少自由度）

# 重置分析矩阵
G <- myGD[, -1]
n <- nrow(G)
m <- ncol(G)
P <- matrix(NA, 1, m)

# 【使用前3个主成分的GLM分析】
for (i in 1:m) {
     x <- G[, i]
     if (max(x) == min(x)) {
          p <- 1
     } else {
          X <- cbind(1, PCA$x[, 1:3], x) # 控制前3个主成分
          LHS <- t(X) %*% X
          C <- solve(LHS)
          RHS <- t(X) %*% y
          b <- C %*% RHS
          yb <- X %*% b
          e <- y - yb
          n <- length(y)
          ve <- sum(e^2) / (n - 1)
          vt <- C * ve
          t <- b / sqrt(diag(vt))
          p <- 2 * (1 - pt(abs(t), n - 2))
     }
     P[i] <- p[length(p)]
}

# 【QQ图比较】
p.obs <- P[!index1to5]
m2 <- length(p.obs)
p.uni <- runif(m2, 0, 1)
order.obs <- order(p.obs)
order.uni <- order(p.uni)

plot(-log10(p.uni[order.uni]), -log10(p.obs[order.obs]),
     ylim = c(0, 7),
     main = "QQ Plot - GLM with 3 PCs",
     xlab = "Expected -log10(P)", ylab = "Observed -log10(P)"
)
abline(a = 0, b = 1, col = "red")

# 【全基因组QQ图】
p.obs <- P # 使用全部SNP
m2 <- length(p.obs)
p.uni <- runif(m2, 0, 1)
order.obs <- order(p.obs)
order.uni <- order(p.uni)

plot(-log10(p.uni[order.uni]), -log10(p.obs[order.obs]),
     main = "QQ Plot - All SNPs with 3 PCs",
     xlab = "Expected -log10(P)", ylab = "Observed -log10(P)"
)
abline(a = 0, b = 1, col = "red")

# 【最终曼哈顿图】
color.vector <- rep(c("deepskyblue", "orange", "forestgreen", "indianred3"), 10)
m <- nrow(myGM)
plot(t(-log10(P)) ~ seq(1:m),
     col = color.vector[myGM[, 2]],
     main = "Manhattan Plot - GLM with Population Structure Control",
     xlab = "SNP Index", ylab = "-log10(P-value)"
)
abline(v = mySim$QTN.position, lty = 2, lwd = 2, col = "black")

# =============================================================================
# 第十部分：亲缘关系矩阵（Kinship Matrix）
# =============================================================================

# 【生物学背景】
# 亲缘关系矩阵量化个体间的遗传相似性
# 在GWAS中用于控制由亲缘关系导致的虚假关联
# 特别重要在近交或具有复杂家系结构的群体中

# 【手动计算亲缘关系矩阵】
myGD <- read.table(file = "mdp_numeric.txt", head = T)
X <- myGD[, -1] # 基因型矩阵
p <- colMeans(X) / 2 # 等位基因频率
M <- X - 1 # 中心化基因型矩阵
Z <- t(M) - 2 * (p - .5) # 标准化基因型矩阵
K <- crossprod((Z), (Z)) # 计算亲缘关系矩阵
adj <- 2 * sum(p * (1 - p)) # 标准化因子
K <- K / adj # 标准化

# 【比较知名品系】
taxa <- myGD[, 1] # 品系名称
favorite <- c("33-16", "38-11", "B73", "B73HTRHM", "CM37", "CML333", "MO17", "YU796NS")
index <- taxa %in% favorite # 找到知名品系索引
snps <- myGD[, -1]

# 显示手动计算的亲缘关系
print("Manual Kinship Calculation:")
print(K[index, index])

# 【使用GAPIT函数计算亲缘关系矩阵】

# VanRaden方法（最常用）
K1 <- GAPIT.kinship.VanRaden(snps)
print("VanRaden Kinship Method:")
print(K1[index, index])

# Zhang方法
K2 <- GAPIT.kinship.Zhang(snps)
print("Zhang Kinship Method:")
print(K2[index, index])

# 【亲缘关系矩阵可视化】
# 使用热图展示亲缘关系模式
heatmap.2(K1,
     cexRow = .2, cexCol = 0.2, col = rev(heat.colors(256)),
     scale = "none", symkey = FALSE, trace = "none",
     main = "Kinship Matrix - VanRaden Method"
)

quartz() # 新建图形窗口
heatmap.2(K2,
     cexRow = .2, cexCol = 0.2, col = rev(heat.colors(256)),
     scale = "none", symkey = FALSE, trace = "none",
     main = "Kinship Matrix - Zhang Method"
)

# 【比较不同计算方法】
n <- nrow(myGD)
ind.a <- seq(1:(n * n)) # 所有矩阵元素索引
i <- 1:n
j <- (i - 1) * n
ind.d <- i + j # 对角线元素索引

# 三个子图比较
par(mfrow = c(1, 3))
# 所有元素的比较
plot(K2[ind.a], K1[ind.a],
     main = "All elements",
     xlab = "Zhang", ylab = "VanRaden"
)
lines(K2[ind.d], K1[ind.d], col = "red", type = "p") # 突出显示对角线元素

# 对角线元素（自身亲缘关系）
plot(K2[ind.d], K1[ind.d],
     main = "Diagonals",
     xlab = "Zhang", ylab = "VanRaden"
)

# 非对角线元素（个体间亲缘关系）
plot(K2[-ind.d], K1[-ind.d],
     main = "Off diag",
     xlab = "Zhang", ylab = "VanRaden"
)

# =============================================================================
# 总结和下一步
# =============================================================================
#
# 【本节课学到的内容】
# 1. GWAS的基本原理和统计方法
# 2. 相关性分析：最简单但容易产生假阳性的方法
# 3. 线性模型（GLM）：通过协变量控制混杂因素
# 4. 主成分分析：识别和控制群体结构
# 5. QQ图：诊断分析质量和假阳性问题
# 6. 亲缘关系矩阵：量化个体间遗传相似性
#
# 【下一步学习方向】
# 1. 混合线性模型（MLM）：同时控制群体结构和亲缘关系
# 2. 多重检验校正：控制假发现率
# 3. 效应量估计：基因解释的表型变异比例
# 4. 连锁不平衡分析：精细定位候选基因
#
# 【思考题】
# 1. 为什么需要控制群体结构？
# 2. 主成分分析和亲缘关系矩阵有什么区别？
# 3. 如何判断GWAS分析的质量？
# 4. 在不同的群体中（如人类、植物、动物），GWAS分析有什么差异？
# =============================================================================
