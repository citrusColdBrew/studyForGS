# 基因组选择(GS)实战代码
# 基于玉米多样性面板数据的完整GS分析流程

# =============================================================================
# 第一步：环境准备和数据载入
# =============================================================================

# 清理环境
rm(list = ls())

# 载入数据（请根据实际路径调整）
cat("载入数据...\n")
genotype_data <- read.table("data/mdp_numeric.txt", header = TRUE)
snp_info <- read.table("data/mdp_SNP_information.txt", header = TRUE)
trait_data <- read.table("data/mdp_YRef.txt", header = TRUE)

# 载入必要的包
required_packages <- c("rrBLUP", "BGLR", "ggplot2", "caret", "corrplot", "glmnet")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    cat(sprintf("安装包: %s\n", pkg))
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 载入GAPIT函数
source("function/gapit_functions.R")

# 查看数据基本信息
cat("=== 数据概况 ===\n")
cat(sprintf("个体数: %d\n", nrow(genotype_data)))
cat(sprintf("SNP数: %d\n", ncol(genotype_data) - 1))
cat(sprintf("表型数据个体数: %d\n", nrow(trait_data)))

# =============================================================================
# 第二步：数据预处理和匹配
# =============================================================================

cat("\n=== 数据预处理 ===\n")

# 匹配基因型和表型数据的个体
genotype_individuals <- genotype_data[, 1]
trait_individuals <- trait_data[, 1]
common_individuals <- intersect(genotype_individuals, trait_individuals)

cat(sprintf("基因型数据个体数: %d\n", length(genotype_individuals)))
cat(sprintf("表型数据个体数: %d\n", length(trait_individuals)))
cat(sprintf("共同个体数: %d\n", length(common_individuals)))

# 根据共同个体筛选数据
genotype_matched_idx <- match(common_individuals, genotype_individuals)
trait_matched_idx <- match(common_individuals, trait_individuals)

matched_genotype_data <- genotype_data[genotype_matched_idx, ]
matched_trait_data <- trait_data[trait_matched_idx, ]

# 验证个体名是否匹配
if(!all(matched_genotype_data[, 1] == matched_trait_data[, 1])) {
  stop("个体名匹配失败！")
}
cat("个体名匹配成功！\n")

# 提取基因型矩阵和表型向量
X_all <- as.matrix(matched_genotype_data[, -1])
y_all <- matched_trait_data[, 2]
rownames(X_all) <- matched_genotype_data[, 1]
names(y_all) <- matched_trait_data[, 1]

# 基本质量控制：移除单态SNP
polymorphic_mask <- apply(X_all, 2, function(snp) {
  unique_values <- unique(snp[!is.na(snp)])
  length(unique_values) > 1
})

X_filtered <- X_all[, polymorphic_mask]
snp_info_filtered <- snp_info[polymorphic_mask, ]

cat(sprintf("质控后SNP数: %d\n", ncol(X_filtered)))

# 表型数据基本统计
cat("\n=== 表型数据统计 ===\n")
cat(sprintf("个体数: %d\n", length(y_all)))
cat(sprintf("均值: %.3f\n", mean(y_all, na.rm = TRUE)))
cat(sprintf("标准差: %.3f\n", sd(y_all, na.rm = TRUE)))
cat(sprintf("范围: %.3f - %.3f\n", min(y_all, na.rm = TRUE), max(y_all, na.rm = TRUE)))

# =============================================================================
# 第三步：训练测试集划分
# =============================================================================

cat("\n=== 训练测试集划分 ===\n")

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
X_train <- X_filtered[training_idx, ]
X_test <- X_filtered[testing_idx, ]
y_train <- y_all[training_idx]
y_test <- y_all[testing_idx]

# 个体名称
taxa_train <- matched_genotype_data[training_idx, 1]
taxa_test <- matched_genotype_data[testing_idx, 1]

# =============================================================================
# 第四步：基因组选择方法
# =============================================================================

cat("\n=== 基因组选择方法比较 ===\n")

# 存储所有方法的结果
gs_methods <- character()
gs_predictions <- list()
gs_accuracies <- numeric()

# 方法1：Ridge Regression BLUP (RR-BLUP)
cat("1. Ridge Regression BLUP (RR-BLUP)...\n")
rrblup_model <- mixed.solve(y = y_train, Z = X_train)
pred_rrblup <- X_test %*% rrblup_model$u
acc_rrblup <- cor(y_test, pred_rrblup)

gs_methods <- c(gs_methods, "RR-BLUP")
gs_predictions[[length(gs_predictions) + 1]] <- as.vector(pred_rrblup)
gs_accuracies <- c(gs_accuracies, acc_rrblup)

cat(sprintf("   预测精度: %.3f\n", acc_rrblup))

# 方法2：基因组BLUP (gBLUP)
cat("2. 基因组BLUP (gBLUP)...\n")
K_train <- tcrossprod(X_train) / ncol(X_train)
gblup_model <- mixed.solve(y = y_train, K = K_train)
K_test_train <- tcrossprod(X_test, X_train) / ncol(X_train)
pred_gblup <- K_test_train %*% gblup_model$u
acc_gblup <- cor(y_test, pred_gblup)

gs_methods <- c(gs_methods, "gBLUP")
gs_predictions[[length(gs_predictions) + 1]] <- as.vector(pred_gblup)
gs_accuracies <- c(gs_accuracies, acc_gblup)

cat(sprintf("   预测精度: %.3f\n", acc_gblup))

# 验证RR-BLUP和gBLUP的等价性
equiv_cor <- cor(X_filtered %*% rrblup_model$u,
                c(gblup_model$u, pred_gblup)[match(rownames(X_filtered),
                                                  c(names(gblup_model$u), taxa_test))])
cat(sprintf("   RR-BLUP与gBLUP等价性验证: %.6f\n", equiv_cor))

# 方法3：使用GAPIT的gBLUP
cat("3. GAPIT gBLUP...\n")
train_GD <- data.frame(Taxa = taxa_train, X_train)
train_Y <- data.frame(Taxa = taxa_train, Trait = y_train)

# 创建完整的基因型数据（包含测试集）
all_GD <- data.frame(Taxa = matched_genotype_data[, 1], X_filtered)

gapit_gblup <- GAPIT(
  Y = train_Y,
  GD = all_GD,  # 使用完整数据集以便预测所有个体
  GM = snp_info_filtered,
  model = "gBLUP",
  SNP.test = FALSE,
  PCA.total = 3
)

# 提取测试集的预测值
pred_results <- gapit_gblup$Pred
test_pred_idx <- match(taxa_test, pred_results$Taxa)
pred_gapit <- pred_results[test_pred_idx, ncol(pred_results)]
acc_gapit <- cor(y_test, pred_gapit)

gs_methods <- c(gs_methods, "GAPIT_gBLUP")
gs_predictions[[length(gs_predictions) + 1]] <- as.vector(pred_gapit)
gs_accuracies <- c(gs_accuracies, acc_gapit)

cat(sprintf("   预测精度: %.3f\n", acc_gapit))

# 方法4：贝叶斯Ridge回归 (BRR)
cat("4. 贝叶斯Ridge回归 (BRR)...\n")
# 准备主成分协变量
pca_result <- prcomp(X_train, center = TRUE, scale. = FALSE)
pc_scores_train <- pca_result$x[, 1:3]
pc_scores_test <- predict(pca_result, X_test)[, 1:3]

# 设置MCMC参数
nIter <- 1500
burnIn <- 300

set.seed(99164)
brr_model <- BGLR(
  y = y_train,
  ETA = list(
    list(X = pc_scores_train, model = 'FIXED'),
    list(X = X_train, model = 'BRR')
  ),
  nIter = nIter,
  burnIn = burnIn,
  verbose = FALSE
)

pred_brr <- pc_scores_test %*% brr_model$ETA[[1]]$b +
            X_test %*% brr_model$ETA[[2]]$b
acc_brr <- cor(y_test, pred_brr)

gs_methods <- c(gs_methods, "BRR")
gs_predictions[[length(gs_predictions) + 1]] <- as.vector(pred_brr)
gs_accuracies <- c(gs_accuracies, acc_brr)

cat(sprintf("   预测精度: %.3f\n", acc_brr))

# 方法5：贝叶斯LASSO (BL)
cat("5. 贝叶斯LASSO (BL)...\n")
set.seed(99164)
bl_model <- BGLR(
  y = y_train,
  ETA = list(
    list(X = pc_scores_train, model = 'FIXED'),
    list(X = X_train, model = 'BL')
  ),
  nIter = nIter,
  burnIn = burnIn,
  verbose = FALSE
)

pred_bl <- pc_scores_test %*% bl_model$ETA[[1]]$b +
           X_test %*% bl_model$ETA[[2]]$b
acc_bl <- cor(y_test, pred_bl)

gs_methods <- c(gs_methods, "BL")
gs_predictions[[length(gs_predictions) + 1]] <- as.vector(pred_bl)
gs_accuracies <- c(gs_accuracies, acc_bl)

cat(sprintf("   预测精度: %.3f\n", acc_bl))

# 方法6：BayesA
cat("6. BayesA...\n")
set.seed(99164)
bayesa_model <- BGLR(
  y = y_train,
  ETA = list(
    list(X = pc_scores_train, model = 'FIXED'),
    list(X = X_train, model = 'BayesA')
  ),
  nIter = nIter,
  burnIn = burnIn,
  verbose = FALSE
)

pred_bayesa <- pc_scores_test %*% bayesa_model$ETA[[1]]$b +
               X_test %*% bayesa_model$ETA[[2]]$b
acc_bayesa <- cor(y_test, pred_bayesa)

gs_methods <- c(gs_methods, "BayesA")
gs_predictions[[length(gs_predictions) + 1]] <- as.vector(pred_bayesa)
gs_accuracies <- c(gs_accuracies, acc_bayesa)

cat(sprintf("   预测精度: %.3f\n", acc_bayesa))

# 方法7：岭回归 (作为神经网络的替代)
cat("7. 岭回归 (Ridge Regression)...\n")
# 数据标准化
X_train_scaled <- scale(X_train)
X_test_scaled <- scale(X_test,
                       center = attr(X_train_scaled, "scaled:center"),
                       scale = attr(X_train_scaled, "scaled:scale"))
y_train_scaled <- scale(y_train)

# 处理缺失值
X_train_scaled[is.na(X_train_scaled)] <- 0
X_test_scaled[is.na(X_test_scaled)] <- 0

# 使用交叉验证选择最优lambda
cv_ridge <- cv.glmnet(X_train_scaled, y_train_scaled, alpha = 0)
ridge_model <- glmnet(X_train_scaled, y_train_scaled,
                     alpha = 0, lambda = cv_ridge$lambda.min)

# 预测并反标准化
pred_ridge_scaled <- predict(ridge_model, newx = X_test_scaled)
pred_ridge <- pred_ridge_scaled * attr(y_train_scaled, "scaled:scale") +
              attr(y_train_scaled, "scaled:center")
acc_ridge <- cor(y_test, pred_ridge)

gs_methods <- c(gs_methods, "Ridge")
gs_predictions[[length(gs_predictions) + 1]] <- as.vector(pred_ridge)
gs_accuracies <- c(gs_accuracies, acc_ridge)

cat(sprintf("   预测精度: %.3f\n", acc_ridge))

# =============================================================================
# 第五步：结果汇总和比较
# =============================================================================

cat("\n=== 方法性能比较 ===\n")

# 创建结果汇总表
results_summary <- data.frame(
  Method = gs_methods,
  Accuracy = round(gs_accuracies, 3),
  R_squared = round(gs_accuracies^2, 3),
  stringsAsFactors = FALSE
)

# 按精度排序
results_summary <- results_summary[order(results_summary$Accuracy, decreasing = TRUE), ]
print(results_summary)

# 找出最佳方法
best_method_idx <- which.max(gs_accuracies)
cat(sprintf("\n最佳方法: %s (精度: %.3f, R² = %.3f)\n",
            gs_methods[best_method_idx],
            gs_accuracies[best_method_idx],
            gs_accuracies[best_method_idx]^2))

# =============================================================================
# 第六步：结果可视化
# =============================================================================

cat("\n=== 结果可视化 ===\n")

# 1. 预测精度比较柱状图
plot_data <- data.frame(
  Method = factor(results_summary$Method,
                  levels = results_summary$Method),
  Accuracy = results_summary$Accuracy
)

p1 <- ggplot(plot_data, aes(x = Method, y = Accuracy, fill = Method)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  geom_text(aes(label = Accuracy), vjust = -0.3, size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "不同基因组选择方法的预测精度比较",
       y = "预测精度 (相关系数)",
       x = "方法") +
  ylim(0, max(gs_accuracies) * 1.1)

print(p1)

# 2. 最佳方法的预测效果散点图
best_prediction <- gs_predictions[[best_method_idx]]

p2 <- ggplot(data.frame(Observed = y_test, Predicted = best_prediction),
            aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.7, color = "blue", size = 2) +
  geom_smooth(method = "lm", color = "red", se = TRUE, alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(title = paste("最佳方法预测效果:", gs_methods[best_method_idx]),
       x = "观测值",
       y = "预测值") +
  annotate("text",
           x = min(y_test) + 0.1 * (max(y_test) - min(y_test)),
           y = max(best_prediction) - 0.1 * (max(best_prediction) - min(best_prediction)),
           label = sprintf("r = %.3f\nR² = %.3f",
                          gs_accuracies[best_method_idx],
                          gs_accuracies[best_method_idx]^2),
           hjust = 0, vjust = 1, size = 4,
           color = "darkred") +
  theme_minimal()

print(p2)

# 3. 方法间相关性热图
if(length(gs_predictions) >= 3) {
  # 创建预测结果矩阵
  pred_matrix <- do.call(cbind, gs_predictions)
  colnames(pred_matrix) <- gs_methods

  # 计算相关矩阵
  cor_matrix <- cor(pred_matrix)

  # 绘制热图
  corrplot(cor_matrix, method = "color", type = "upper",
           order = "hclust", tl.cex = 0.8, tl.col = "black",
           title = "不同GS方法预测结果的相关性",
           mar = c(0,0,1,0))
}

# =============================================================================
# 第七步：简化的交叉验证
# =============================================================================

cat("\n=== 3折交叉验证评估 ===\n")

# 选择主要方法进行交叉验证
cv_methods <- c("RR-BLUP", "gBLUP", "BRR")
k_folds <- 3

set.seed(123)
folds <- createFolds(y_all, k = k_folds, returnTrain = TRUE)

cv_results <- data.frame(
  Fold = integer(),
  Method = character(),
  Accuracy = numeric(),
  stringsAsFactors = FALSE
)

for(fold in 1:k_folds) {
  cat(sprintf("处理第%d折...\n", fold))

  # 获取训练和验证索引
  train_fold_idx <- folds[[fold]]
  val_fold_idx <- setdiff(1:length(y_all), train_fold_idx)

  # 提取数据
  X_fold_train <- X_filtered[train_fold_idx, ]
  y_fold_train <- y_all[train_fold_idx]
  X_fold_val <- X_filtered[val_fold_idx, ]
  y_fold_val <- y_all[val_fold_idx]

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

  # BRR (快速版)
  pca_fold <- prcomp(X_fold_train, center = TRUE, scale. = FALSE)
  pc_fold_train <- pca_fold$x[, 1:3]
  pc_fold_val <- predict(pca_fold, X_fold_val)[, 1:3]

  brr_fold <- BGLR(
    y = y_fold_train,
    ETA = list(
      list(X = pc_fold_train, model = 'FIXED'),
      list(X = X_fold_train, model = 'BRR')
    ),
    nIter = 800, burnIn = 200, verbose = FALSE
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
cv_summary$Mean_Accuracy <- round(cv_summary$Mean_Accuracy, 3)
cv_summary$SD_Accuracy <- round(cv_summary$SD_Accuracy, 3)

cat("\n交叉验证结果汇总:\n")
print(cv_summary)

# 交叉验证可视化
p3 <- ggplot(cv_results, aes(x = Method, y = Accuracy, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4,
               fill = "red", color = "black") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "3折交叉验证结果",
       y = "预测精度",
       x = "方法")

print(p3)

# =============================================================================
# 第八步：分析总结
# =============================================================================

cat("\n=== 基因组选择分析总结 ===\n")
cat("本次GS分析完成了以下内容:\n")
cat("1. 数据预处理和训练测试集划分\n")
cat("2. 多种基因组选择方法比较\n")
cat("3. 预测精度评估和可视化\n")
cat("4. 交叉验证性能验证\n")
cat("5. 方法间相关性分析\n")

cat("\n主要发现:\n")
cat(sprintf("- 测试了%d种基因组选择方法\n", length(gs_methods)))
cat(sprintf("- 最佳方法: %s (精度: %.3f)\n",
            gs_methods[best_method_idx], gs_accuracies[best_method_idx]))
cat("- 不同方法在预测精度上存在差异\n")
cat("- 贝叶斯方法能够有效整合先验信息\n")
cat("- 线性混合模型方法计算效率较高\n")

cat("\n实际应用建议:\n")
cat("- 根据群体特征选择合适的GS方法\n")
cat("- 考虑计算资源和预测精度的平衡\n")
cat("- 定期更新训练群体以维持预测精度\n")
cat("- 结合表型测定进行模型验证\n")
cat("- 在育种决策中综合考虑多种信息\n")

cat("\n技术要点:\n")
cat("- 训练群体与候选群体的遗传关系是关键\n")
cat("- 标记密度影响预测精度但存在边际递减效应\n")
cat("- 性状遗传力直接影响GS效果\n")
cat("- 交叉验证是评估GS方法的重要手段\n")

cat("\n分析完成！\n")
cat("========================================\n")