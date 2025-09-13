#----Part 1: 基础概念模拟 - 二项分布------#
quartz()
# 在Mac环境下打开一个新的图形窗口。在RStudio等集成环境中通常不是必需的。
par(mfrow = c(4, 4), mar = c(3, 4, 1, 1))
# 设置图形参数:
# - mfrow = c(4, 4): 创建一个4x4的网格，后续的图形将按行填充到这个网格中。
# - mar = c(3, 4, 1, 1): 设置图形的边距（下、左、上、右），单位是行。

n <- 10000
m <- 100
# n: 定义样本量，即模拟10000个个体。
# m: 定义每个个体事件的数量，可以理解为有100个基因座。

x <- rbinom(n, m, .5)
# - rbinom(): 从二项分布中生成随机数。
# - n: 生成10000个随机数。
# - m: 每次试验的总次数，即100个基因座。
# - .5: 每次试验成功的概率，可以理解为等位基因频率为0.5。
# 这行代码模拟了10000个个体，每个个体有100个基因，每个基因有两种等位基因（如A/a），频率均为0.5。x向量存储了每个个体拥有"A"等位基因的总数。
hist(x)
# Histogram-直方图
# 绘制变量x的直方图，展示这10000个个体等位基因总数的分布情况。
# 根据中心极限定理，这个分布会非常接近正态分布。
plot(density(x))
# 绘制变量x的核密度估计图，这是对真实分布的平滑估计，可以更清晰地看出分布的形状。


#------Part 2: 手动模拟基因型和表型-------#
x <- rep(1, m)
# 生成一个包含100个1的向量。
# x=runif(m)

gene <- matrix(x, n, m, byrow = T)
# 创建一个 n x m (10000 x 100) 的矩阵，并将向量x按行填充。
# 结果是一个所有元素都为1的矩阵，可以理解为初始状态下所有个体在所有位点上的基因型都是纯合显性（AA）。
head(gene)
galton <- matrix(runif(n * m), n, m)
# 生成一个 10000 x 100 的矩阵，其元素为0到1之间的均匀分布随机数。
head(galton)
galton.binary <- galton < .5
# 将galton矩阵中的每个元素与0.5进行比较。如果小于0.5，则为TRUE，否则为FALSE。
# 这相当于为每个个体的每个基因座随机分配一个等位基因（抛硬币决定）。
head(galton.binary)
gene[galton.binary] <- 0
# 利用逻辑索引。在`galton.binary`矩阵中为TRUE的位置，将`gene`矩阵中对应位置的元素赋值为0。
# 这一步完成了基因型模拟：初始全为1的矩阵，经过随机的0/1替换，最终形成了一个代表两种纯合基因型（比如AA=1, aa=0）的矩阵。
head(gene)
y <- rowSums(gene)
# 对`gene`矩阵的每一行求和。这模拟了一个简单的加性模型，即个体的表型值（y）等于其所有位点“优良”等位基因的总和。

y[1:6]
hist(y)
plot(density(y))
# 同样地，绘制模拟出的表型值y的直方图和密度图，其分布也应接近正态分布。


#----Part 3: 利用真实基因型数据模拟复杂性状------#
myGD <- read.table(file = "/Users/sameen/workspace/statistical genomics/data/mdp_numeric.txt", head = T)
# 读取基因型数据（Genotype Data, GD）。
myGM <- read.table(file = "/Users/sameen/workspace/statistical genomics/data/mdp_SNP_information.txt", head = T)
# 读取SNP的遗传图谱信息（Genetic Map, GM）。

# Sampling QTN
NQTN <- 10
# 设置QTN的数量为10。QTN是真正影响目标性状的遗传变异位点（causal variants）。
X <- myGD[, -1]
# 从`myGD`中提取纯基因型矩阵，通过 `[,-1]` 排除了第一列（通常是个体ID）。
n <- nrow(X)
m <- ncol(X)
# 获取基因型矩阵的维度：n为个体数，m为SNP标记总数。

QTN.position <- sample(m, NQTN, replace = F)
# 从总共m=3091个SNP中，无放回地（replace = F）随机抽取NQTN（10）个SNP，作为我们预设的QTN。
# `QTN.position` 存储了这些QTN在矩阵X中的列索引。

SNPQ <- as.matrix(X[, QTN.position])
# 从完整基因型矩阵X中，提取出这10个QTN的基因型数据，存为`SNPQ`矩阵。
QTN.position
head(SNPQ)


plot(myGM[, c(2, 3)])
# 绘制一个散点图，展示所有SNP在基因组上的物理位置。通常第2列是染色体，第3列是物理位置（bp）。这构成了曼哈顿图的背景。

lines(myGM[QTN.position, c(2, 3)], type = "p", col = "red")
points(myGM[QTN.position, c(2, 3)], type = "p", col = "blue", cex = 5)
# 在图上高亮标记出我们抽取的QTN的位置。
# `lines`和`points`在这里效果类似，都是在指定位置画点。用蓝色大点（cex=5）醒目地标出了QTN。

# --- 模拟QTN效应和表型 ---
addeffect <- rnorm(NQTN, 0, 1)
# 为10个QTN分别模拟一个加性效应值。这些效应值从均值为0，标准差为1的正态分布中抽取。
addeffect

effect <- SNPQ %*% addeffect
# 这是模拟的核心步骤，计算每个个体的真实遗传值（Genetic value）。
# 通过矩阵乘法 (`%*%`)，将每个个体的QTN基因型矩阵 (`SNPQ`) 与QTN效应值向量 (`addeffect`) 相乘。
# 对于每个个体，其遗传值 = sum(genotype_qtn_i * effect_qtn_i) for i in 1 to NQTN。

head(effect)

h2 <- .7
# 设定目标遗传力（narrow-sense heritability, h²）为0.7。遗传力表示表型变异中可以由加性遗传效应解释的比例。

effectvar <- var(effect)
# 计算由QTN效应产生的遗传方差 (Variance of Genetic effect, Vg or Va)。

effectvar
residualvar <- (effectvar - h2 * effectvar) / h2
# 这是根据遗传力公式 h² = Vg / (Vg + Ve) 推导出来的。
# 我们已知 h² 和 Vg (`effectvar`)，需要计算出相应的环境方差 (Variance of error, Ve, 即`residualvar`)。
# 推导: Vp = Vg + Ve => Vp = Vg / h² => Ve = Vp - Vg = (Vg / h²) - Vg = Vg * (1/h² - 1) = Vg * (1-h²)/h²
# 代码中的 `(effectvar - h2 * effectvar) / h2` 等价于 `effectvar * (1-h2) / h2`，公式正确。

residualvar
residual <- rnorm(n, 0, sqrt(residualvar))
# 模拟环境效应（或称残差效应）。从均值为0，标准差为`sqrt(residualvar)`的正态分布中为每个个体抽取一个残差值。
head(residual)

y <- effect + residual
# 最终的表型值y由遗传效应`effect`和环境效应`residual`相加得到，即 P = G + E 模型。

# --- 可视化模拟的表型 ---
par(mfrow = c(2, 2))
plot(y)
# 绘制表型值的散点图（按个体顺序）。
hist(y)
# 绘制表型值的直方图。
boxplot(y)
# 绘制表型值的箱线图。
plot(ecdf(y))
# 绘制表型值的经验累积分布函数图(ECDF)。


plot(density(y), ylim = c(0, .5))
lines(density(effect), col = "blue")
lines(density(residual), col = "red")

# Check on heritability探究不同遗传力的影响
source("/Users/sameen/workspace/statistical genomics/function/G2P.R")
par(mfrow = c(3, 1), mar = c(3, 4, 1, 1))

# 案例一: 高遗传力 h² = 0.75
myG2P <- G2P(X, .75, 1, 10, "norm")
str(myG2P)
va <- var(myG2P$add)
# 计算加性遗传方差
ve <- var(myG2P$residual)
# 计算环境方差
vp <- var(myG2P$y)
# 计算总表型方差
v <- matrix(c(va, ve, vp), 1, 3)
colnames(v) <- c("A", "E", "P")
barplot(v, col = "gray")
# 绘制方差组分条形图

# 案例二: 中等遗传力 h² = 0.5
myG2P <- G2P(X, .5, 1, 10, "norm")
va <- var(myG2P$add)
ve <- var(myG2P$residual)
vp <- var(myG2P$y)
v <- matrix(c(va, ve, vp), 1, 3)
colnames(v) <- c("A", "E", "P")
barplot(v, col = "gray")

# 案例三: 低遗传力 h² = 0.25
myG2P <- G2P(X, .25, 1, 10, "norm")
va <- var(myG2P$add)
ve <- var(myG2P$residual)
vp <- var(myG2P$y)
v <- matrix(c(va, ve, vp), 1, 3)
colnames(v) <- c("A", "E", "P")
barplot(v, col = "gray")

# Check on number of genes
par(mfrow = c(3, 3), mar = c(3, 4, 1, 1))

# 案例一: NQTN = 2 (寡基因)
myG2P <- G2P(X, .75, 1, 2, "norm")
hist(myG2P$add)
# 遗传效应的直方图
hist(myG2P$residual)
# 环境效应的直方图
hist(myG2P$y)
# 最终表型的直方图

# 案例二: NQTN = 10
myG2P <- G2P(X, .75, 1, 10, "norm")
hist(myG2P$add)
hist(myG2P$residual)
hist(myG2P$y)

# 案例三: NQTN = 100 (多基因)
myG2P <- G2P(X, .75, 1, 100, "norm")
hist(myG2P$add)
hist(myG2P$residual)
hist(myG2P$y)

# Desect phenotype可视化不同QTN数量对表型组分分布的影响
par(mfrow = c(3, 1))

# NQTN = 2
myG2P <- G2P(X, .5, 1, 2, "norm")
plot(density(myG2P$y), ylim = c(0, 1.5))
lines(density(myG2P$add), col = "blue")
lines(density(myG2P$residual), col = "red")

# NQTN = 10
myG2P <- G2P(X, .5, 1, 10, "norm")
plot(density(myG2P$y), ylim = c(0, .3))
lines(density(myG2P$add), col = "blue")
lines(density(myG2P$residual), col = "red")

# NQTN = 100探究不同QTN效应分布的影响
myG2P <- G2P(X, .5, 1, 100, "norm")
plot(density(myG2P$y), ylim = c(0, .1))
lines(density(myG2P$add), col = "blue")
lines(density(myG2P$residual), col = "red")

# Check gene effect distribution
par(mfrow = c(3, 3), mar = c(3, 4, 1, 1))

# 几何分布，参数可能为1 (极端情况)
myG2P <- G2P(X, .75, 1, 10, "geom")
hist(myG2P$add)
hist(myG2P$residual)
hist(myG2P$y)

# 几何分布，参数可能为0.95
myG2P <- G2P(X, .75, .95, 10, "geom")
hist(myG2P$add)
hist(myG2P$residual)
hist(myG2P$y)

# 几何分布，参数可能为0.5
myG2P <- G2P(X, .75, .5, 10, "geom")
hist(myG2P$add)
hist(myG2P$residual)
hist(myG2P$y)
