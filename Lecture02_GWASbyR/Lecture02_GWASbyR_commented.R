# GWASbyCor：一个R函数，通过相关性分析进行全基因组关联研究（GWAS）。
# 输入：
#   X: 基因型数值矩阵，行代表个体，列代表单核苷酸多态性（SNP）。
#   y: 每个个体的表型数值向量。
# 输出：
#   p: 每个SNP的p值数值向量，表示其与表型关联的统计显著性。
GWASbyCor <- function(X, y) {
  # 'n' 是个体数量，由基因型矩阵'X'的行数确定。
  n <- nrow(X)
  # 'r' 是表型'y'与'X'中每个SNP之间的相关系数。
  r <- cor(y, X)
  # 'n' 是个体数量，由基因型矩阵'X'的行数确定。
  n <- nrow(X)
  # 't' 为每个SNP计算t统计量，以检验相关的显著性。
  t <- r / sqrt((1 - r^2) / (n - 2))
  # 'p' 从t分布计算双尾p值。
  p <- 2 * (1 - pt(abs(t), n - 2))
  # 'zeros' 是一个逻辑向量，用于识别p值为零的情况。
  zeros <- p == 0
  # 这行代码将0的p值替换为一个非常小的数（1e-10），以避免后续计算（如-log10转换）中出现问题。
  p[zeros] <- 1e-10
  # 函数返回p值向量。
  return(p)
}


# 将工作目录设置为指定路径。这是R将查找文件和保存输出的地方。
# 请将此路径替换为您系统上的实际路径。
setwd("/Users/sameen/workspace/statistical genomics/Lecture02_GWASbyR/output")

# 'source' 加载并执行 "gapit_functions.txt" 文件中的R代码，使其函数在当前会话中可用。
source("/Users/sameen/workspace/statistical genomics/function/gapit_functions.R")

# 'myGD' 从 "mdp_numeric.txt" 读取一个表格到数据框中。'head = T' 表示文件的第一行是标题。
myGD <- read.table("/Users/sameen/workspace/statistical genomics/data/mdp_numeric.txt", head = T)
# 'myGM' 从 "mdp_SNP_information.txt" 读取一个表格到数据框中。'head = T' 表示文件的第一行是标题。
myGM <- read.table("/Users/sameen/workspace/statistical genomics/data/mdp_SNP_information.txt", head = T)

# 'source' 加载并执行 "G2P.R" 文件中的R代码。
source("/Users/sameen/workspace/statistical genomics/function/G2P.R")
# 'source' 加载并执行 "GWASbyCor.R" 文件中的R代码。
source("/Users/sameen/workspace/statistical genomics/function/GWASbyCor.R")
# 'X' 通过从'myGD'中移除第一列（个体ID）来创建基因型矩阵。
X <- myGD[, -1]
# 'index1to5' 为位于1到5号染色体上的SNP创建一个逻辑索引。
index1to5 <- myGM[, 2] < 6
# 'X1to5' 对基因型矩阵'X'进行子集化，只包括来自1到5号染色体的SNP。
X1to5 <- X[, index1to5]
# 'set.seed' 将随机数生成器的种子设置为一个特定值（99164），以确保模拟的可复现性。
set.seed(99164)
# 'mySim' 使用'G2P'函数模拟一个表型，参数包括遗传力（h2）、效应大小（alpha）和数量性状核苷酸（NQTN）的数量。
mySim <- G2P(X = X1to5, h2 = .75, alpha = 1, NQTN = 10, distribution = "norm")
# 'p' 使用'GWASbyCor'函数计算'X'中所有SNP与模拟表型'mySim$y'之间关联的p值。
p <- GWASbyCor(X = X, y = mySim$y)

# 'color.vector' 创建一个颜色向量，用于根据SNP所在的染色体进行绘图。
color.vector <- rep(c("deepskyblue", "orange", "forestgreen", "indianred3"), 10)
# 'm' 是SNP的总数，由SNP信息矩阵'myGM'的行数确定。
m <- nrow(myGM)
# 这行代码创建一个曼哈顿图，通过绘制每个SNP的-log10转换后的p值来可视化GWAS结果。
plot(t(-log10(p)) ~ seq(1:m), col = color.vector[myGM[, 2]])
# 'abline' 在曼哈顿图上添加垂直虚线，以指示模拟的数量性状核苷酸（QTN）的位置。
abline(v = mySim$QTN.position, lty = 2, lwd = 2, col = "black")

# 'p.obs' 提取不在1到5号染色体上的SNP的p值（即无效SNP）。
p.obs <- p[!index1to5]
# 'm2' 是这些无效SNP的数量。
m2 <- length(p.obs)
# 'p.uni' 生成一组与无效SNP同样大小的均匀分布的随机p值，用于创建Q-Q图。
p.uni <- runif(m2, 0, 1)
# 'order.obs' 获取将无效SNP的观测p值按升序排序的索引。
order.obs <- order(p.obs)
# 'order.uni' 获取将均匀分布的p值按升序排序的索引。
order.uni <- order(p.uni)

# 这会创建一个分位数-分位数（Q-Q）图，该图比较观测p值的分布与在零假设下预期的均匀分布。
plot(
  -log10(p.uni[order.uni]),
  -log10(p.obs[order.obs])
)
# 'abline' 在Q-Q图上添加一条红色参考线（y=x）。与这条线的偏离表明一些观测到的p值比偶然预期的要小。
abline(a = 0, b = 1, col = "red")

# 'order.obs' 获取将无效SNP的观测p值按升序排序的索引。
order.obs <- order(p.obs)
# 'X6to10' 对基因型矩阵'X'进行子集化，只包括不在1到5号染色体上的SNP。
X6to10 <- X[, !index1to5]
# 'Xtop' 选择无效SNP中p值最小的SNP的基因型数据。
Xtop <- X6to10[, order.obs[1]]

# 这会创建一个箱线图，以可视化与表型关联最强的SNP（'Xtop'）的基因型与模拟表型（'mySim$y'）之间的关系。
boxplot(mySim$y ~ Xtop)

# 'PCA' 对基因型矩阵'X'进行主成分分析（PCA），以分析群体结构。
PCA <- prcomp(X)
# 这会创建模拟表型（'mySim$y'）与第二主成分（'PCA$x[, 2]'）的散点图。
plot(mySim$y, PCA$x[, 2])
# 'cor' 计算模拟表型与第二主成分之间的相关性。
cor(mySim$y, PCA$x[, 2])

# 'set.seed' 重新设置随机种子以确保可复现性。
set.seed(99164)
# 's' 从数据集中随机抽取10个个体。
s <- sample(length(mySim$y), 10)
# 这为10个个体的子集创建一个散点图，显示其表型与第二主成分之间的关系。
plot(mySim$y[s], PCA$x[s, 2])
# 'cor' 计算这10个个体子集的相关性。
cor(mySim$y[s], PCA$x[s, 2])

# 'y' 将模拟表型赋给一个新变量。
y <- mySim$y
# 'X' 为线性模型创建一个设计矩阵，包括截距、第二主成分和关联最强的SNP。
X <- cbind(1, PCA$x[, 2], Xtop)
# 'LHS' 计算正规方程的左侧（X'X）。
LHS <- t(X) %*% X
# 'C' 计算LHS矩阵的逆，这是求解模型系数所必需的。
C <- solve(LHS)
# 'RHS' 计算正规方程的右侧（X'y）。
RHS <- t(X) %*% y
# 'b' 计算线性模型的估计系数（beta-hat）。
b <- C %*% RHS
# 'yb' 根据拟合模型计算预测的表型值。
yb <- X %*% b
# 'e' 计算残差（观测值与预测表型值之差）。
e <- y - yb
# 'n' 是个体数量。
n <- length(y)
# 've' 计算残差的估计方差。
ve <- sum(e^2) / (n - 1)
# 'vt' 计算估计系数的方差-协方差矩阵。
vt <- C * ve
# 't' 为每个系数计算t统计量，以检验其显著性。
t <- b / sqrt(diag(vt))
# 'p' 从t分布为每个系数计算双尾p值。
p <- 2 * (1 - pt(abs(t), n - 2))

# 'LM' 将估计的系数、t统计量、标准差和p值合并到一个矩阵中。
LM <- cbind(b, t, sqrt(diag(vt)), p)
# 'rownames' 为'LM'矩阵的行分配名称。
rownames(LM) <- cbind("Mean", "PC2", "Xtop")
# 'colnames' 为'LM'矩阵的列分配名称。
colnames(LM) <- cbind("b", "t", "SD", "p")
# 这会打印'LM'矩阵，显示线性模型的结果。
LM

# 'G' 将基因型数据（不包括第一列）赋给一个新变量。
G <- myGD[, -1]
# 'n' 是个体数量。
n <- nrow(G)
# 'm' 是SNP的数量。
m <- ncol(G)
# 'P' 初始化一个矩阵，用于存储GWAS循环中的p值。
P <- matrix(NA, 1, m)

# 这个循环遍历每个SNP，使用包含第二主成分作为协变量的线性模型进行GWAS，以校正群体结构。
for (i in 1:m) {
  # 'x' 是当前SNP的基因型数据。
  x <- G[, i]
  # 这会检查SNP是否有任何变异；如果没有，p值设置为1。
  if (max(x) == min(x)) {
    p <- 1
  } else {
    # 'X' 为线性模型创建设计矩阵，包括截距、第二主成分和当前SNP。
    X <- cbind(1, PCA$x[, 2], x)
    # 'LHS' 计算正规方程的左侧。
    LHS <- t(X) %*% X
    # 'C' 计算LHS矩阵的逆。
    C <- solve(LHS)
    # 'RHS' 计算正规方程的右侧。
    RHS <- t(X) %*% y
    # 'b' 计算估计的系数。
    b <- C %*% RHS
    # 'yb' 计算预测的表型值。
    yb <- X %*% b
    # 'e' 计算残差。
    e <- y - yb
    # 'n' 是个体数量。
    n <- length(y)
    # 've' 计算残差的估计方差。
    ve <- sum(e^2) / (n - 1)
    # 'vt' 计算系数的方差-协方差矩阵。
    vt <- C * ve
    # 't' 为每个系数计算t统计量。
    t <- b / sqrt(diag(vt))
    # 'p' 为每个系数计算双尾p值。
    p <- 2 * (1 - pt(abs(t), n - 2))
  } # 结束变异测试
  # 这会将当前SNP的p值（'p'向量的最后一个元素）存储在'P'矩阵中。
  P[i] <- p[length(p)]
} # 结束标记循环

# 'p.obs' 提取无效SNP（不在1-5号染色体上）的p值。
p.obs <- P[!index1to5]
# 'm2' 是这些无效SNP的数量。
m2 <- length(p.obs)
# 'p.uni' 为Q-Q图生成一组均匀分布的随机p值。
p.uni <- runif(m2, 0, 1)
# 'order.obs' 获取对观测p值进行排序的索引。
order.obs <- order(p.obs)
# 'order.uni' 获取对均匀p值进行排序的索引。
order.uni <- order(p.uni)

# 这为GWAS结果创建一个Q-Q图，y轴限制为7。
plot(-log10(p.uni[order.uni]),
  -log10(p.obs[order.obs]),
  ylim = c(0, 7)
)
# 'abline' 在Q-Q图上添加一条红色参考线。
abline(a = 0, b = 1, col = "red")

# 'G' 重新分配基因型数据。
G <- myGD[, -1]
# 'n' 是个体数量。
n <- nrow(G)
# 'm' 是SNP的数量。
m <- ncol(G)
# 'P' 重新初始化p值矩阵。
P <- matrix(NA, 1, m)

# 这个循环与前一个类似，但它在线性模型中包括前三个主成分作为协变量，以更稳健地校正群体结构。
for (i in 1:m) {
  # 'x' 是当前SNP的基因型数据。
  x <- G[, i]
  # 这会检查SNP的变异。
  if (max(x) == min(x)) {
    p <- 1
  } else {
    # 'X' 创建包含截距和前三个主成分的设计矩阵。
    X <- cbind(1, PCA$x[, 1:3], x)
    # 以下行执行线性模型回归，与之前相同。
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
  } # 结束变异测试
  # 这会存储当前SNP的p值。
  P[i] <- p[length(p)]
} # 结束标记循环

# 以下行为使用三个主成分的GWAS结果生成另一个Q-Q图。
p.obs <- P[!index1to5]
m2 <- length(p.obs)
p.uni <- runif(m2, 0, 1)
order.obs <- order(p.obs)
order.uni <- order(p.uni)

plot(-log10(p.uni[order.uni]),
  -log10(p.obs[order.obs]),
  ylim = c(0, 7)
)
abline(a = 0, b = 1, col = "red")

# 以下行为所有SNP（不仅仅是无效SNP）生成一个Q-Q图。
p.obs <- P
m2 <- length(p.obs)
p.uni <- runif(m2, 0, 1)
order.obs <- order(p.obs)
order.uni <- order(p.uni)

plot(
  -log10(p.uni[order.uni]),
  -log10(p.obs[order.obs]),
)
abline(a = 0, b = 1, col = "red")

# 以下行为使用线性模型（包括三个主成分）获得的GWAS结果生成一个新的曼哈顿图。
color.vector <- rep(c("deepskyblue", "orange", "forestgreen", "indianred3"), 10)
m <- nrow(myGM)
plot(t(-log10(P)) ~ seq(1:m), col = color.vector[myGM[, 2]])
# 'abline' 添加垂直线以指示模拟QTN的位置。
abline(v = mySim$QTN.position, lty = 2, lwd = 2, col = "black")

# 'myGD' 再次读取数值基因型数据。
myGD <- read.table(file = "data/mdp_numeric.txt", head = T)
# 'X' 提取基因型矩阵。
X <- myGD[, -1]
# 'p' 计算每个SNP的等位基因频率。
p <- colMeans(X) / 2
# 'M' 中心化基因型矩阵。
M <- X - 1
# 'Z' 标准化基因型矩阵。
Z <- t(M) - 2 * (p - .5)
# 'K' 使用标准化的基因型计算基因组关系矩阵（亲缘关系矩阵）。
K <- crossprod((Z), (Z))
# 'adj' 计算基于等位基因频率的调整因子。
adj <- 2 * sum(p * (1 - p))
# 'K' 标准化亲缘关系矩阵。
K <- K / adj

# 这些行被注释掉了，但如果取消注释，它们将加载额外的GAPIT函数。
# source("/Users/Jiabo/Dropbox/GAPIT/Functions/GAPIT.library.R")
# source("/Users/Jiabo/Dropbox/GAPIT/Functions/gapit_functions.txt")

# 'myGD' 再次读取数值基因型数据。
myGD <- read.table(file = "data/mdp_numeric.txt", head = T)
# 'taxa' 提取个体ID。
taxa <- myGD[, 1]
# 'favorite' 创建一个特定个体ID的向量。
favorite <- c("33-16", "38-11", "B73", "B73HTRHM", "CM37", "CML333", "MO17", "YU796NS")
# 'index' 为'favorite'中列出的个体创建一个逻辑索引。
index <- taxa %in% favorite
# 'snps' 提取基因型矩阵。
snps <- myGD[, -1]

# 这行被注释掉了，但它会使用GAPIT的Loiselle方法计算一个亲缘关系矩阵。
# K=GAPIT.kinship.loiselle(t(myGD[,-1]), method="additive", use="all")
# 这将显示'favorite'个体的亲缘关系值。
K[index, index]

# 'K1' 使用GAPIT中实现的VanRaden方法计算亲缘关系矩阵。
K1 <- GAPIT.kinship.VanRaden(snps)
# 这使用VanRaden方法显示'favorite'个体的亲缘关系值。
K1[index, index]

# 'K2' 使用GAPIT中实现的Zhang方法计算亲缘关系矩阵。
K2 <- GAPIT.kinship.Zhang(snps)
# 这使用Zhang方法显示'favorite'个体的亲缘关系值。
K2[index, index]

# 'heatmap.2' 生成VanRaden亲缘关系矩阵的热图。
heatmap.2(K1, cexRow = .2, cexCol = 0.2, col = rev(heat.colors(256)), scale = "none", symkey = FALSE, trace = "none")
# 'quartz()' 打开一个新的图形设备（在macOS上）。
quartz()
# 'heatmap.2' 生成Zhang亲缘关系矩阵的热图。
heatmap.2(K2, cexRow = .2, cexCol = 0.2, col = rev(heat.colors(256)), scale = "none", symkey = FALSE, trace = "none")

# 'n' 是个体数量。
n <- nrow(myGD)
# 'ind.a' 为亲缘关系矩阵中的所有元素创建一个索引序列。
ind.a <- seq(1:(n * n))
# 'i' 是从1到n的序列。
i <- 1:n
# 'j' 计算扁平化矩阵中每行的起始索引。
j <- (i - 1) * n
# 'ind.d' 为亲缘关系矩阵的对角线元素创建索引。
ind.d <- i + j
# 'par(mfrow = c(1, 3))' 设置绘图区域为1行3列。
par(mfrow = c(1, 3))
# 这会创建一个散点图，比较Zhang和VanRaden亲缘关系矩阵的所有元素。
plot(K2[ind.a], K1[ind.a], main = "All elements", xlab = "Zhang", ylab = "VanRaden")
# 'lines' 以红色点的形式在图上添加对角线元素。
lines(K2[ind.d], K1[ind.d], main = "All elements", xlab = "Zhang", ylab = "VanRaden", col = "red", type = "p")
# 这会创建一个散点图，只比较两个亲缘关系矩阵的对角线元素。
plot(K2[ind.d], K1[ind.d], main = "Diagonals", xlab = "Zhang", ylab = "VanRaden")
# 这会创建一个散点图，只比较两个亲缘关系矩阵的非对角线元素。
plot(K2[-ind.d], K1[-ind.d], main = "Off diag", xlab = "Zhang", ylab = "VanRaden")
