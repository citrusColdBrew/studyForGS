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

### 生物学背景
GS（Genomic Selection，全基因组选择）是利用全基因组标记信息预测个体基因组估计育种值（GEBV）的方法。

**核心理念**：
- 传统选择依赖表型观测，GS依赖基因组信息
- 通过训练群体建立"基因型-表型"预测模型
- 利用模型预测候选群体的育种价值

**与GWAS的区别**：
- **GWAS目标**：发现与性状关联的基因/QTN
- **GS目标**：预测个体的基因组育种值（GEBV）
- **GWAS关注**：少数大效应位点的显著性检验
- **GS关注**：所有标记的累积效应预测

### GS的优势
1. **时间效率**：无需等待表型成熟即可选择
2. **成本效率**：减少大规模表型测定需求
3. **选择强度**：可在幼体阶段进行选择
4. **遗传增益**：缩短世代间隔，提高选择效率

---

## 2. 数据准备和理解

### 我们的数据集
1. **mdp_numeric.txt**：玉米多样性面板的基因型数据
2. **mdp_SNP_information.txt**：SNP位置信息
3. **mdp_traits.txt**：真实农艺性状表型数据
4. **mdp_YRef.txt**：参考表型数据（可用于模拟研究）

### GS数据需求
**必需数据**：
- **训练群体**：有基因型和表型的个体
- **候选群体**：只有基因型的待选个体

**数据特点**：
- 标记密度：通常需要覆盖全基因组的高密度标记
- 群体大小：训练群体越大，预测精度通常越高
- 亲缘关系：训练群体与候选群体的亲缘关系影响预测精度

---

## 3. 第一步：数据载入和预处理

### 3.1 设置工作环境

```r
# 清理环境
rm(list = ls())

# 设置工作目录
setwd("/Users/your_path/statistical genomics")

# 载入必要的函数库
source("function/gapit_functions.R")

# 检查必要的包
required_packages <- c("rrBLUP", "BGLR", "caret", "ggplot2", "corrplot")
for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
```

### 3.2 载入数据

```r
# 载入基因型数据
cat("载入基因型数据...\n")
genotype_data <- read.table("data/mdp_numeric.txt", header = TRUE, sep = "\t")

# 载入SNP信息
cat("载入SNP位置信息...\n")
snp_info <- read.table("data/mdp_SNP_information.txt", header = TRUE, sep = "\t")

# 载入表型数据（选择真实性状或参考性状）
# 选项1：使用真实农艺性状
# trait_data <- read.table("data/mdp_traits.txt", header = TRUE, sep = "\t")

# 选项2：使用参考表型（用于本教程演示）
trait_data <- read.table("data/mdp_YRef.txt", header = TRUE, sep = "\t")

cat("数据载入完成！\n")
```

### 3.3 数据基本检查

```r
# 检查数据维度
cat("=== 数据基本信息 ===\n")
cat(sprintf("个体数量: %d\n", nrow(genotype_data)))
cat(sprintf("SNP数量: %d\n", ncol(genotype_data) - 1))
cat(sprintf("表型数据个体数: %d\n", nrow(trait_data)))

# 检查基因型数据前几行
cat("\n=== 基因型数据示例 ===\n")
print(genotype_data[1:5, 1:6])

# 检查表型数据
cat("\n=== 表型数据示例 ===\n")
print(head(trait_data))

# 表型数据基本统计
phenotype_values <- trait_data[, 2]
cat("\n=== 表型基本统计 ===\n")
cat(sprintf("个体数: %d\n", length(phenotype_values)))
cat(sprintf("均值: %.3f\n", mean(phenotype_values, na.rm = TRUE)))
cat(sprintf("标准差: %.3f\n", sd(phenotype_values, na.rm = TRUE)))
cat(sprintf("范围: %.3f - %.3f\n",
            min(phenotype_values, na.rm = TRUE),
            max(phenotype_values, na.rm = TRUE)))
```

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