# Genomic Prediction using Neural Network with 5-fold Cross-Validation

# Load required libraries
library(keras)        # For building and training the neural network
library(tensorflow)   # Backend for Keras
library(data.table)   # For efficient data handling
library(ggplot2)      # For plotting training history
library(dplyr)        # For data manipulation
library(caret)        # For creating cross-validation folds

library(reticulate)
use_python("/opt/miniconda3/envs/forGS/bin/python")
use_condaenv("forGS", required = TRUE)
py_config()
# Step 1: Load and preprocess data
# Load genotype (SNP) data
setwd("/Users/sameen/workspace/statistical genomics/code_day3")
# reading data and source functions
source("gapit_functions.txt")
myGD=read.table(file="mdp_numeric.txt",head=T) 
myGM=read.table(file="mdp_SNP_information.txt",head=T)
set.seed(99164)
n=nrow(myGD)
testing=sample(n,round(n/5),replace=F)
training=-testing

#Simultate 20 QTN on the first half chromosomes
X=myGD[,-1]
taxa=myGD[,1]

set.seed(99164)

mySim<-GAPIT(h2=0.75,NQTN=20,GD= myGD,GM=myGM,PCA.total=3)
myQTN=mySim$QTN.position
myY=mySim$Y
myCV=mySim$PCA

# Prepare genotype data
genotypes <- myGD
genotypes <- genotypes[, -1]  # Remove ID column

# Perform PCA on genotype data to capture population structure
pcs <- prcomp(myGD[,-1])  # Exclude ID column
pca <- as.data.frame(pcs$x[, 1:3])  # First three PCA components

covariates <- pca

# Combine genotypes and covariates into input matrix X
X <- as.matrix(cbind(genotypes, covariates))

# Normalize the input data
X <- scale(X)
X[is.na(X)] <- 0  # Replace NAs with 0

# Load and normalize phenotype data
y <- myY$Sim
y <- scale(y)

# Step 2: Set up 5-fold cross-validation
set.seed(123)
folds <- createFolds(y, k = 5, returnTrain = TRUE)

# Initialize lists to store predictions and metrics
predictions <- list()
correlations <- numeric(5)
rmse_values <- numeric(5)

# Step 3: Train and predict with cross-validation
for (fold_idx in seq_along(folds)) {
  cat("Training fold", fold_idx, "of 5\n")
  
  # Split data into training and validation sets
  train_idx <- folds[[fold_idx]]
  val_idx <- setdiff(1:nrow(X), train_idx)
  
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_val <- X[val_idx, ]
  y_val <- y[val_idx]
  
  # Step 4: Build the neural network model
  model <- keras_model_sequential() %>%
    layer_dense(units = 48, activation = "relu", input_shape = ncol(X), 
                kernel_regularizer = regularizer_l2(0.08)) %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = 0.25) %>%
    layer_dense(units = 24, activation = "relu", 
                kernel_regularizer = regularizer_l2(0.08)) %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = 0.25) %>%
    layer_dense(units = 1)
  
  # Step 5: Compile the model
  model %>% keras::compile(
    optimizer = optimizer_adam(learning_rate = 0.0003),
    loss = "mean_squared_error"
  )
  
  # Step 6: Define callbacks
  lr_scheduler <- callback_reduce_lr_on_plateau(
    monitor = "val_loss",
    factor = 0.3,
    patience = 5,
    min_lr = 1e-6
  )
  
  early_stopping <- callback_early_stopping(
    monitor = "val_loss",
    patience = 5,
    restore_best_weights = TRUE
  )
  
  # Step 7: Train the model
  history <- model %>% fit(
    X_train, y_train,
    epochs = 500,
    batch_size = 154,
    validation_data = list(X_val, y_val),
    callbacks = list(early_stopping, lr_scheduler),
    verbose = 1
  )
  
  # Step 8: Predict on validation set
  y_pred <- model %>% predict(X_val, batch_size = 154)
  
  # Store predictions
  predictions[[fold_idx]] <- data.frame(
    Fold = fold_idx,
    True = y_val,
    Predicted = as.vector(y_pred)
  )
  
  # Compute evaluation metrics
  correlations[fold_idx] <- cor(y_val, y_pred)^2
  rmse_values[fold_idx] <- sqrt(mean((y_val - y_pred)^2))
  
  # Plot training history
##Delete
  #png(filename = paste0("training_history_fold_", fold_idx, ".png"))
  #plot(history)
  #dev.off()
## Delete
 
##add
library(ggplot2)

# Transform training history into data.frame
df_history <- as.data.frame(history)

# Plot
p <- ggplot(df_history, aes(x = epoch, y = value, color = data)) +
  geom_line(size = 1.2) +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(
    title = paste("Training History - Fold", fold_idx),
    x = "Epoch",
    y = "Metric"
       )

# Save image
ggsave(
  filename = paste0("training_history_fold_", fold_idx, ".png"),
  plot = p,
  width = 6,
  height = 5,
  dpi = 150
      )
##add  
  
  
}

# Step 9: Combine predictions and compute average metrics
predictions_df <- do.call(rbind, predictions)
# write.csv(predictions_df, "genomic_predictions.csv", row.names = FALSE)
write.csv(predictions_df, "genomic_predictions_epi.csv", row.names = FALSE)

# Summary of prediction performance
cat("Prediction performance across folds:\n")
cat("Mean correlation:", mean(correlations), "\n")
cat("Correlation SD:", sd(correlations), "\n")
cat("Mean RMSE:", mean(rmse_values), "\n")
cat("RMSE SD:", sd(rmse_values), "\n")

#-----------------------------------------------
# Correlation Plots for Validation Population in Each Fold
# Requires predictions from genomic_prediction.R

# Load required library
library(ggplot2)

# Assume predictions_df is available from genomic_prediction.R
# Alternatively, load from file if running separately
# predictions_df <- read.csv("genomic_predictions.csv")
predictions_df <- read.csv("genomic_predictions_epi.csv")

# Create correlation plots for each fold
##add
corr <- numeric(5)
##add

for (fold_idx in 1:5) {
  # Subset data for the current fold
  fold_data <- predictions_df[predictions_df$Fold == fold_idx, ]
  
  # Extract true and predicted values
  y_true <- fold_data$True
  y_pred <- fold_data$Predicted
  
  # Calculate Pearson correlation （add 164）
  corr[fold_idx] <- cor(y_true, y_pred)^2
  
  # Create scatter plot
  p <- ggplot(fold_data, aes(x = True, y = Predicted)) +
    geom_point(color = "blue", alpha = 0.6) +
    geom_smooth(method = "lm", color = "red", se = FALSE) +  # Add regression line
    labs(
      title = paste("Fold", fold_idx, ": Predicted vs. True Phenotypes"),
      x = "True Phenotype",
      y = "Predicted Phenotype"
    ) +
    annotate(
      "text",
      x = min(y_true) + 0.1 * (max(y_true) - min(y_true)),
      y = max(y_pred) - 0.1 * (max(y_pred) - min(y_pred)),
      label = sprintf("Correlation: %.3f", corr[fold_idx]),
      hjust = 0,
      vjust = 1
    ) +
    theme_minimal()
  
  # Save the plot
  ggsave(
    # filename = paste0("correlation_plot_fold_", fold_idx, ".png"),
    filename = paste0("correlation_plot_fold_epi_", fold_idx, ".png"),
    plot = p,
    width = 6,
    height = 5,
    dpi = 150
  )
}

cat("Correlation across folds:\n")
cat("Mean Cor:", mean(corr), "\n")
#-------------------------------------------------


# # Step 10: Train final model on all data (optional)
# model <- keras_model_sequential() %>%
#   layer_dense(units = 48, activation = "relu", input_shape = ncol(X), 
#               kernel_regularizer = regularizer_l2(0.08)) %>%
#   layer_batch_normalization() %>%
#   layer_dropout(rate = 0.25) %>%
#   layer_dense(units = 24, activation = "relu", 
#               kernel_regularizer = regularizer_l2(0.08)) %>%
#   layer_batch_normalization() %>%
#   layer_dropout(rate = 0.25) %>%
#   layer_dense(units = 1)
# 
# model %>% compile(
#   optimizer = optimizer_adam(learning_rate = 0.0003),
#   loss = "mean_squared_error"
# )
# 
# history <- model %>% fit(
#   X, y,
#   epochs = 500,
#   batch_size = 281,
#   callbacks = list(
#     callback_early_stopping(monitor = "loss", patience = 5, restore_best_weights = TRUE),
#     callback_reduce_lr_on_plateau(monitor = "loss", factor = 0.3, patience = 5, min_lr = 1e-6)
#   ),
#   verbose = 1
# )

# Save the final model
save_model_hdf5(model, "final_genomic_prediction_model.h5")

# # Load the final model and predict phenotypes for new genomic data
# model <- load_model_hdf5("final_genomic_prediction_model.h5")
# X_new <- scale(read.table("new_genomic_data.txt")[,-1])  # Preprocess new data
# y_pred_new <- model %>% predict(X_new)
# write.csv(y_pred_new, "new_predictions.csv")