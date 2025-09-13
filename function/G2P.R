# Genotype to Phenotype Simulation
# This function is to simulate phenotype based on genotype
# Input:
# GD: Genotype in numeric format, pure 0,1,2 matrix. Need to be imputed
# 纯数字（0, 1, 2）编码的基因型矩阵。
# h2: Heritability 遗传力（Heritability）
# D: Dominance显性度, not implemented yet(尚未实现)
# NQTN: Number of QTN 数量性状核苷酸的数量。
# dist: Distribution of QTN effects, options are "norm" and "geom"
# QTN效应的分布。可选值为 "norm" (正态分布) 和 "geom" (几何分布)。
# a: parameter for distribution 分布的参数。
# 这个在函数定义中并没有出现，但在调用时被放在了第3个参数的位置，我们可以理解为上次分析的参数 D 实际上就是这里的 a。
# Output:
# y: phenotype 模拟的表型值
# add: additive effect 加性效应值
# dom: dominance effect, not implemented yet  显性效应值 (尚未实现)
# residual: residual effect 残差效应值
# Author: Zhiwu Zhang
# Last update: Febuary 26, 2016

# Genotype to Phenotype Simulation (功能：基因型到表型的模拟)。
G2P <- function(X, h2, alpha, NQTN, distribution) {
  # 函数定义:
  # 定义了名为 G2P 的函数。
  # GD: 基因型数据，这是唯一一个没有默认值的参数，必须由用户提供。
  # h2 = .75: 遗传力，默认值为0.75。
  # D = 1: 第三个参数，我们现在知道它在"geom"分布时起作用，默认值为1。
  # NQTN = 10: QTN数量，默认值为10。
  # dist = "norm": QTN效应分布，默认为正态分布。
  n <- nrow(X) # 个体数量 (行数)。
  m <- ncol(X) # SNP标记总数 (列数)。
  # Sampling QTN
  QTN.position <- sample(m, NQTN, replace = FALSE)
  # 从总共 m 个SNP中，无放回地随机抽取 NQTN 个作为QTN。
  SNPQ <- as.matrix(X[, QTN.position])
  # 从基因型矩阵GD中提取出这些QTN所在的列，确保提取出的QTN基因型数据 SNPQ 是矩阵格式
  QTN.position

  # QTN effects
  if (distribution == "norm") {
    # 如果分布是正态分布
    addeffect <- rnorm(NQTN, 0, 1)
    # 从均值为0，标准差为1的标准正态分布中生成 NQTN 个随机数，作为每个QTN的加性效应。
  } else {
    # 如果分布是几何分布
    addeffect <- alpha^(1:NQTN)
    # 1:NQTN 会生成一个从1到NQTN的序列 (1, 2, 3, ...)。
    # D^(1:NQTN) 会计算 D 的1次方, 2次方, 3次方...
    # 这模拟了一种效应值按几何级数递减的模式。
    # 当 D (0 < D < 1) 接近1时（如0.95），效应值下降缓慢；当 D 较小时（如0.5），效应值下降非常快。
    # 这实现了“少数大效应、多数微效”的遗传结构。例如，D=0.5, NQTN=10时，效应值为 (0.5, 0.25, 0.125, ...)。
  }
  # Simulate phenotype
  effect <- SNPQ %*% addeffect
  # 通过矩阵乘法 (%*%) 计算每个个体的总加性遗传值 effect
  effectvar <- var(effect)
  # 计算遗传值的方差 (Vg)。
  residualvar <- (effectvar - h2 * effectvar) / h2
  # 计算与目标 h2 和 effectvar 相匹配的环境方差 (Ve)
  residual <- rnorm(n, 0, sqrt(residualvar))
  # 从均值为0，标准差为 sqrt(residualvar) 的正态分布中，为每个个体生成一个残差值。
  y <- effect + residual
  # 将遗传效应 effect 和残差效应 residual 相加，得到最终的表型 y。
  return(
    list(
      addeffect = addeffect,
      y = y,
      add = effect,
      residual = residual,
      QTN.position = QTN.position,
      SNPQ = SNPQ
    )
  )
}
