# =====================================================================
# Script 02: Multi-Dimensional Machine Learning Pipeline
# Description: LASSO, Random Forest, SVM-RFE consensus, and XGBoost
# =====================================================================

library(glmnet)
library(randomForest)
library(e1071)
library(caret)
library(xgboost)
library(pROC)
library(ggplot2)
library(VennDiagram)

# 1. Data Preparation
# Assuming all_sig_genes is extracted from previous virtual knockout results
hub_candidates <- unique(all_sig_genes)
expr_mat_small <- count_matrix_filtered[rownames(count_matrix_filtered) %in% hub_candidates, ]
X_small <- t(as.matrix(expr_mat_small))
X_small <- as.data.frame(X_small)

# Generate binary labels based on RPL34 median expression
set.seed(123)
Y_small <- ifelse(X_small[, "RPL34"] > median(X_small[, "RPL34"]), 1, 0)
Y_small <- as.factor(Y_small)

# 2. LASSO Regression
set.seed(123)
cv_fit <- cv.glmnet(as.matrix(X_small), Y_small, family = "binomial", alpha = 1, nfolds = 5)
lasso_coef <- coef(cv_fit, s = "lambda.min")
lasso_genes <- rownames(lasso_coef)[which(lasso_coef != 0)]
lasso_genes <- setdiff(lasso_genes, "(Intercept)")

# 3. Random Forest
set.seed(123)
rf_model <- randomForest(x = X_small, y = Y_small, ntree = 300, importance = TRUE)
importance_df <- as.data.frame(importance(rf_model))
importance_df$Gene <- rownames(importance_df)
importance_df <- importance_df[order(importance_df$MeanDecreaseGini, decreasing = TRUE), ]
rf_genes <- head(importance_df$Gene, 15)

# 4. SVM-RFE (Custom implementation)
svm_data <- as.data.frame(cbind(Y = as.factor(Y_small), X_small))
gene_names <- colnames(X_small)
while (length(gene_names) > 5) {
  model <- svm(Y ~ ., data = svm_data[, c("Y", gene_names)], kernel = "linear", cost = 1)
  weights <- t(model$coefs) %*% model$SV
  worst_gene <- names(which.min(abs(weights)))
  gene_names <- setdiff(gene_names, worst_gene)
}
svm_genes <- gene_names

# 5. Consensus Intersection (Hub Gene Identification)
final_hub <- Reduce(intersect, list(lasso_genes, rf_genes, svm_genes))
print(paste("Identified Hub Gene(s):", final_hub))

# 6. XGBoost Modeling & ROC Evaluation
dtrain <- xgb.DMatrix(data = as.matrix(X_small), label = as.numeric(as.character(Y_small)))
params <- list(objective = "binary:logistic", eval_metric = "auc", max_depth = 3, eta = 0.1)
xgb_model <- xgb.train(params = params, data = dtrain, nrounds = 50, verbose = 0)
preds <- predict(xgb_model, dtrain)
roc_obj <- roc(as.numeric(as.character(Y_small)), preds)