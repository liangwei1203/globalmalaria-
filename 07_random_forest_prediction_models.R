# ============================================================
# File: 07_random_forest_prediction_models.R
# Original file: 2026.3.16新数据库随机森林代码.R
# Purpose: Random-forest prediction models, performance metrics, variable importance, and diagnostic plots
# Study: Imported malaria analysis
#
# Inputs:
# - introduced260318_imputed.csv
#
# Outputs:
# - 新再传播模型变量重要性3.26.csv
# - 新死亡模型4.7默认model_full_metrics.csv
# - 新再传播模型3.26评估指标_SeSp_cutoff_AUC95CI.csv
# - 新再传播模型3.26评估指标_PPV_NPV_cutoff_AUC95CI.csv
# - 新死亡模型3.26评估指标_Youden指数_cutoff_AUC95CI.csv
# - 新再传播3.26样本预测结果.csv
# - Train_Sensitivity_Specificity_Accuracy.pdf
# - 死亡Train_PPV_NPV_Accuracy.pdf
# - 死亡Train_Youden_Accuracy.pdf
# - 0414再传播Train_Sensitivity_Specificity_Accuracy.pdf
# - 0414再传播Train_PPV_NPV_Accuracy.pdf
# - 0414再传播Train_Youden_Accuracy.pdf
#
# Notes:
# - Personal working-directory paths have been removed.
# - Run from the repository root using source("R/07_random_forest_prediction_models.R").
# - Original unmodified script is retained in original_scripts/.
# ============================================================

source(file.path("R", "00_setup.R"))

# setwd() removed; project paths are defined in R/00_setup.R
mydata <- read.csv(input_file("introduced260318_imputed.csv"))
# View(mydata)  # Interactive inspection removed for reproducibility
library(dplyr)
library(caret)
library(randomForest)
library(pROC)

# 查看每列的缺失值数量
missing_values <- colSums(is.na(mydata))
cat("\n每列的缺失值数量:\n")
print(missing_values)

# 删除包含缺失值的行
mydata <- na.omit(mydata)

# 首先确保 Temp 和 Precip 是数值类型
mydata$Temp <- as.numeric(mydata$Temp)
mydata$Precip <- as.numeric(mydata$Precip)
mydata$Mosquito_borne_species <- as.numeric(mydata$Mosquito_borne_species)

# 使用 preProcess 进行归一化
preprocess_params <- preProcess(mydata[, c("Temp", "Precip","Mosquito_borne_species")], method = "range")
mydata[, c("Temp", "Precip","Mosquito_borne_species")] <- predict(preprocess_params, mydata[, c("Temp", "Precip","Mosquito_borne_species")])

# 检查归一化后的数值范围
cat("Temp 归一化后的范围:", range(mydata$Temp, na.rm = TRUE), "\n")
cat("Precip 归一化后的范围:", range(mydata$Precip, na.rm = TRUE), "\n")
cat("Mosquito_borne_species 归一化后的范围:", range(mydata$Mosquito_borne_species, na.rm = TRUE), "\n")

# 因子转换
mydata$Outcome = as.factor(mydata$Outcome)
mydata$Region_of_acquisition = as.factor(mydata$Region_of_acquisition)
mydata$Country_of_destination = as.factor(mydata$Country_of_destination)
mydata$Plasmodium_species = as.factor(mydata$Plasmodium_species)
mydata$Occupation = as.factor(mydata$Occupation)
mydata$Groups_of_economies = as.factor(mydata$Groups_of_economies)

# 指定需要筛选的变量
selected_columns <- c("Region_of_acquisition", "Country_of_destination", 
                      "Plasmodium_species", "Temp", "Precip", 
                      "Mosquito_borne_species","Occupation","Groups_of_economies")

# 选择指定的变量
mydata <- mydata %>% select(all_of(selected_columns), Outcome)

# 分离特征和目标变量
X <- mydata %>% select(all_of(selected_columns))
y <- mydata$Outcome

# 使用分层抽样划分训练集和测试集
set.seed(42)
train_indices <- createDataPartition(y, p = 0.8, list = FALSE, times = 1)

# 划分训练集和测试集
X_train <- X[train_indices, ]
y_train <- y[train_indices]
X_test <- X[-train_indices, ]
y_test <- y[-train_indices]

# 检查类别分布
cat("原始训练集类别分布:\n")
print(table(y_train))
cat("原始测试集类别分布:\n")
print(table(y_test))

# 递归特征消除(RFE)函数
rfe_feature_selection <- function(X_train, y_train, initial_vars) {
  current_vars <- initial_vars
  results <- list()
  n_iter <- 0
  
  # 迭代直到只剩2个变量
  while(length(current_vars) > 2) {
    n_iter <- n_iter + 1
    cat("迭代", n_iter, ": 当前变量数 =", length(current_vars), "\n")
    
    # 设置树的数量（第一次5000，后续2000）
    ntree <- ifelse(n_iter == 1, 5000, 2000)
    
    # 训练随机森林模型
    set.seed(42)
    rf_model <- randomForest(
      x = X_train[, current_vars, drop = FALSE],
      y = y_train,
      ntree = ntree,
      importance = TRUE,
      keep.forest = TRUE
    )
    
    # 获取变量重要性（基于准确率下降）
    importance_df <- as.data.frame(rf_model$importance)
    importance_df$Variable <- rownames(importance_df)
    
    # 使用MeanDecreaseAccuracy作为重要性指标
    if ("MeanDecreaseAccuracy" %in% colnames(importance_df)) {
      importance_df <- importance_df[order(importance_df$MeanDecreaseAccuracy), ]
    } else if ("MeanDecreaseGini" %in% colnames(importance_df)) {
      importance_df <- importance_df[order(importance_df$MeanDecreaseGini), ]
    } else {
      # 如果上述列名都不存在，使用第一列
      importance_df <- importance_df[order(importance_df[, 1]), ]
    }
    
    # 计算要移除的变量数量（20%）
    n_remove <- max(1, round(length(current_vars) * 0.2))
    vars_to_remove <- head(importance_df$Variable, n_remove)
    
    # 存储结果
    results[[n_iter]] <- list(
      iteration = n_iter,
      n_vars = length(current_vars),
      variables = current_vars,
      oob_error = tail(rf_model$err.rate[, "OOB"], 1),
      importance = importance_df,
      removed_vars = vars_to_remove
    )
    
    # 移除最不重要的变量
    current_vars <- setdiff(current_vars, vars_to_remove)
  }
  
  return(results)
}

# 执行RFE
initial_vars <- selected_columns
rfe_results <- rfe_feature_selection(X_train, y_train, initial_vars)

# 找到最优变量子集（一个标准误差法则）
find_optimal_model <- function(rfe_results) {
  # 提取所有模型的OOB误差
  oob_errors <- sapply(rfe_results, function(x) x$oob_error)
  n_vars <- sapply(rfe_results, function(x) x$n_vars)
  
  # 找到最小误差
  min_error <- min(oob_errors)
  min_error_index <- which.min(oob_errors)
  
  # 计算最小误差模型的标准误差
  # 这里我们使用所有模型误差的标准差作为近似的标准误差
  se <- sd(oob_errors) / sqrt(length(oob_errors))
  
  # 找到变量数最少且误差在最小误差一个标准误差范围内的模型
  threshold <- min_error + se
  candidate_models <- which(oob_errors <= threshold)
  
  if (length(candidate_models) > 0) {
    # 选择变量最少的模型
    optimal_index <- candidate_models[which.min(n_vars[candidate_models])]
    optimal_model <- rfe_results[[optimal_index]]
    
    cat("最小OOB误差:", round(min_error, 4), "\n")
    cat("标准误差:", round(se, 4), "\n")
    cat("阈值:", round(threshold, 4), "\n")
    cat("最优模型: 迭代", optimal_model$iteration, "\n")
    cat("变量数量:", optimal_model$n_vars, "\n")
    cat("OOB误差:", round(optimal_model$oob_error, 4), "\n")
    
    return(optimal_model)
  } else {
    # 如果没有满足条件的模型，返回最小误差模型
    cat("没有找到满足一个标准误差法则的模型，返回最小误差模型\n")
    return(rfe_results[[min_error_index]])
  }
}

# 找到最优模型
optimal_model <- find_optimal_model(rfe_results)

# 输出最终选择的变量
final_variables <- optimal_model$variables
cat("\n最终选择的变量 (", length(final_variables), "个):\n")
print(final_variables)

# 使用最终变量训练模型
final_rf_model <- randomForest(
  x = X_train[, final_variables, drop = FALSE],
  y = y_train,
  ntree = 2000,
  importance = TRUE
)

#==================== 重要性导出代码 ====================
##提取最优模型的变量重要性
importance_data <- as.data.frame(importance(final_rf_model))
importance_data$Variable <- rownames(importance_data)

# 检查并重命名重要性指标列
colnames(importance_data) <- gsub(" ", "", colnames(importance_data))  # 移除空格

# 确保列名正确
if ("MeanDecreaseAccuracy" %in% colnames(importance_data)) {
  # 列名已经是正确的
} else if ("MeanDecreaseAccuracy" %in% colnames(importance_data)) {
  # 处理可能的列名变体
  colnames(importance_data)[colnames(importance_data) == "MeanDecreaseAccuracy"] <- "MeanDecreaseAccuracy"
}

if ("MeanDecreaseGini" %in% colnames(importance_data)) {
  # 列名已经是正确的
} else if ("MeanDecreaseGini" %in% colnames(importance_data)) {
  # 处理可能的列名变体
  colnames(importance_data)[colnames(importance_data) == "MeanDecreaseGini"] <- "MeanDecreaseGini"
}

# 打印当前列名以便检查
cat("\n重要性数据列名:\n")
print(colnames(importance_data))
# 按MeanDecreaseAccuracy从大到小排序
if ("MeanDecreaseAccuracy" %in% colnames(importance_data)) {
  importance_sorted <- importance_data[order(-importance_data$MeanDecreaseAccuracy), ]
} else {
  # 如果MeanDecreaseAccuracy列不存在，使用第一列数值列排序
  numeric_cols <- sapply(importance_data, is.numeric)
  if (sum(numeric_cols) > 0) {
    first_numeric_col <- names(importance_data)[numeric_cols][1]
    importance_sorted <- importance_data[order(-importance_data[[first_numeric_col]]), ]
    cat("使用", first_numeric_col, "列进行排序\n")
  } else {
    importance_sorted <- importance_data
    cat("没有找到数值列进行排序\n")
  }
}

# 选择需要的列并重命名
if ("MeanDecreaseAccuracy" %in% colnames(importance_sorted) && "MeanDecreaseGini" %in% colnames(importance_sorted)) {
  final_importance <- data.frame(
    变量名 = importance_sorted$Variable,
    MeanDecreaseAccuracy = importance_sorted$MeanDecreaseAccuracy,
    MeanDecreaseGini = importance_sorted$MeanDecreaseGini
  )
} else if ("MeanDecreaseAccuracy" %in% colnames(importance_sorted)) {
  final_importance <- data.frame(
    变量名 = importance_sorted$Variable,
    MeanDecreaseAccuracy = importance_sorted$MeanDecreaseAccuracy
  )
} else if ("MeanDecreaseGini" %in% colnames(importance_sorted)) {
  final_importance <- data.frame(
    变量名 = importance_sorted$Variable,
    MeanDecreaseGini = importance_sorted$MeanDecreaseGini
  )
} else {
  # 如果两个标准列名都不存在，使用所有数值列
  numeric_cols <- sapply(importance_sorted, is.numeric)
  final_importance <- data.frame(
    变量名 = importance_sorted$Variable
  )
  for (col in names(importance_sorted)[numeric_cols]) {
    final_importance[[col]] <- importance_sorted[[col]]
  }
}

# 打印排序后的重要性数据
cat("\n最优模型变量重要性 (按MeanDecreaseAccuracy降序排列):\n")
print(final_importance)

# 导出到CSV文件
write.csv(final_importance, output_table_file("新再传播模型变量重要性3.26.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("\n变量重要性数据已保存到 '最优模型变量重要性.csv'\n")


library(pROC)
library(caret)
library(PRROC)
library(randomForest)

# -------------------------------
# 1️⃣ 完整性能评估函数
# -------------------------------
evaluate_full_metrics <- function(model, X_data, y_data, dataset_name, final_vars = NULL, cutoff = 0.5) {
  
  if (!is.null(final_vars)) {
    X_data <- X_data[, final_vars, drop = FALSE]
  }
  
  # 预测概率和类别
  predictions_prob <- predict(model, X_data, type = "prob")[,2]
  predictions_class <- ifelse(predictions_prob >= cutoff, levels(y_data)[2], levels(y_data)[1])
  predictions_class <- factor(predictions_class, levels = levels(y_data))
  
  # 混淆矩阵
  cm <- confusionMatrix(predictions_class, y_data, positive = levels(y_data)[2])
  
  TP <- cm$table[2,2]; TN <- cm$table[1,1]; FP <- cm$table[1,2]; FN <- cm$table[2,1]
  
  # 常规指标
  sensitivity <- as.numeric(cm$byClass["Sensitivity"])
  specificity <- as.numeric(cm$byClass["Specificity"])
  accuracy <- as.numeric(cm$overall["Accuracy"])
  ppv <- as.numeric(cm$byClass["Pos Pred Value"])
  npv <- as.numeric(cm$byClass["Neg Pred Value"])
  kappa <- as.numeric(cm$overall["Kappa"])
  
  # 平衡准确率 & Youden指数
  balanced_acc <- (sensitivity + specificity)/2
  youden_index <- sensitivity + specificity - 1
  
  # F1 / F2
  F1 <- ifelse((ppv + sensitivity)==0, 0, 2*ppv*sensitivity/(ppv + sensitivity))
  F2 <- ifelse((ppv + sensitivity)==0, 0, (1 + 2^2)*ppv*sensitivity/(2^2*ppv + sensitivity))
  
  # MCC
  MCC <- ifelse((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)==0, NA, (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
  
  # PR AUC
  pr <- pr.curve(scores.class0 = predictions_prob[y_data==levels(y_data)[2]], 
                 scores.class1 = predictions_prob[y_data==levels(y_data)[1]], curve=FALSE)
  pr_auc <- pr$auc.integral
  
  # Prevalence
  prevalence <- mean(y_data == levels(y_data)[2])
  
  # ROC AUC + 95% CI
  roc_obj <- roc(y_data, predictions_prob)
  auc_val <- as.numeric(auc(roc_obj))
  auc_ci <- ci.auc(roc_obj, method="bootstrap")
  
  # 返回完整指标
  results <- data.frame(
    Dataset = dataset_name,
    阈值 = cutoff,
    灵敏度 = round(sensitivity,3),
    特异度 = round(specificity,3),
    准确率 = round(accuracy,3),
    平衡准确率 = round(balanced_acc,3),
    Youden指数 = round(youden_index,3),
    PPV = round(ppv,3),
    NPV = round(npv,3),
    F1 = round(F1,3),
    F2 = round(F2,3),
    MCC = round(MCC,3),
    PR_AUC = round(pr_auc,3),
    Kappa = round(kappa,3),
    Prevalence = round(prevalence,3),
    AUC = round(auc_val,3),
    AUC_95CI_lower = round(auc_ci[1],3),
    AUC_95CI_upper = round(auc_ci[3],3)
  )
  
  return(results)
}

# -------------------------------
# 2️⃣ 训练集和测试集评估
# -------------------------------
train_metrics <- evaluate_full_metrics(final_rf_model, X_train, y_train, "训练集", final_variables, cutoff = 0.5)
test_metrics <- evaluate_full_metrics(final_rf_model, X_test, y_test, "测试集", final_variables, cutoff = 0.5)

all_metrics <- rbind(train_metrics, test_metrics)

# -------------------------------
# 3️⃣ 导出 CSV
# -------------------------------
write.csv(all_metrics, output_table_file("新死亡模型4.7默认model_full_metrics.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("\n✅ 模型完整性能指标已成功导出: model_full_metrics.csv\n")

library(pROC)
library(caret)
library(PRROC)

# -------------------------------
# 1️⃣ 训练集 ROC & Se=Sp cutoff
# -------------------------------
train_probs <- predict(final_rf_model, X_train[, final_variables, drop=FALSE], type="prob")[,2]

roc_train <- roc(y_train, train_probs)
roc_coords <- coords(roc_train, x = "all", ret=c("threshold","sensitivity","specificity"))
roc_coords$diff <- abs(roc_coords$sensitivity - roc_coords$specificity)

best_cutoff_row <- roc_coords[which.min(roc_coords$diff), ]
best_cutoff <- best_cutoff_row$threshold

cat("训练集 Se=Sp 最优 cutoff:", round(best_cutoff,3), "\n")
cat("此时 Sensitivity:", round(best_cutoff_row$sensitivity,3), 
    "Specificity:", round(best_cutoff_row$specificity,3), "\n")

# -------------------------------
# 2️⃣ 完整指标计算函数（新增 Prevalence & AUC & AUC95CI）
# -------------------------------
calculate_full_metrics <- function(y_true, probs, cutoff) {
  # 分类结果
  pred_class <- ifelse(probs >= cutoff, levels(y_true)[2], levels(y_true)[1])
  pred_class <- factor(pred_class, levels=levels(y_true))
  
  # 混淆矩阵
  cm <- confusionMatrix(pred_class, y_true, positive=levels(y_true)[2])
  TP <- cm$table[2,2]; TN <- cm$table[1,1]; FP <- cm$table[1,2]; FN <- cm$table[2,1]
  
  # 常规指标
  sensitivity <- as.numeric(cm$byClass["Sensitivity"])
  specificity <- as.numeric(cm$byClass["Specificity"])
  accuracy <- as.numeric(cm$overall["Accuracy"])
  ppv <- as.numeric(cm$byClass["Pos Pred Value"])
  npv <- as.numeric(cm$byClass["Neg Pred Value"])
  kappa <- as.numeric(cm$overall["Kappa"])
  
  # 平衡准确率 & Youden指数
  balanced_acc <- (sensitivity + specificity)/2
  youden_index <- sensitivity + specificity -1
  
  # F1 / F2
  F1 <- ifelse((ppv + sensitivity)==0, 0, 2*ppv*sensitivity/(ppv + sensitivity))
  F2 <- ifelse((ppv + sensitivity)==0, 0, (1 + 2^2)*ppv*sensitivity/(2^2*ppv + sensitivity))
  
  # MCC
  MCC <- ifelse((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)==0, NA, (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
  
  # PR AUC
  pr <- pr.curve(scores.class0 = probs[y_true==levels(y_true)[2]], 
                 scores.class1 = probs[y_true==levels(y_true)[1]], curve=FALSE)
  pr_auc <- pr$auc.integral
  
  # Prevalence
  prevalence <- mean(y_true == levels(y_true)[2])
  
  # ROC AUC + 95% CI
  roc_obj <- roc(y_true, probs)
  auc_val <- as.numeric(auc(roc_obj))
  auc_ci <- ci.auc(roc_obj, method="bootstrap")
  
  data.frame(
    阈值 = cutoff,
    灵敏度 = round(sensitivity,3),
    特异度 = round(specificity,3),
    准确率 = round(accuracy,3),
    平衡准确率 = round(balanced_acc,3),
    Youden指数 = round(youden_index,3),
    PPV = round(ppv,3),
    NPV = round(npv,3),
    F1 = round(F1,3),
    F2 = round(F2,3),
    MCC = round(MCC,3),
    PR_AUC = round(pr_auc,3),
    Kappa = round(kappa,3),
    Prevalence = round(prevalence,3),
    AUC = round(auc_val,3),
    AUC_95CI_lower = round(auc_ci[1],3),
    AUC_95CI_upper = round(auc_ci[3],3)
  )
}

# -------------------------------
# 3️⃣ 训练集 & 验证集指标
# -------------------------------
test_probs <- predict(final_rf_model, X_test[, final_variables, drop=FALSE], type="prob")[,2]

train_metrics <- calculate_full_metrics(y_train, train_probs, best_cutoff)
train_metrics$数据集 <- "训练集"

test_metrics <- calculate_full_metrics(y_test, test_probs, best_cutoff)
test_metrics$数据集 <- "验证集"

all_metrics <- rbind(train_metrics, test_metrics)
all_metrics <- all_metrics[, c("数据集", names(train_metrics)[1:17])]  # 包含 AUC CI

# -------------------------------
# 4️⃣ 导出 CSV
# -------------------------------
write.csv(all_metrics, output_table_file("新再传播模型3.26评估指标_SeSp_cutoff_AUC95CI.csv"), 
          row.names=FALSE, fileEncoding="UTF-8")
cat("\n✅ 已导出 CSV: 训练验证集完整评估指标_SeSp_cutoff_AUC95CI.csv\n")


# -------------------------------
# 1️⃣ 训练集 ROC & PPV=NPV cutoff
# -------------------------------
train_probs <- predict(final_rf_model, X_train[, final_variables, drop=FALSE], type="prob")[,2]

roc_train <- roc(y_train, train_probs)

# 获取所有阈值下的指标
roc_coords <- coords(roc_train, x = "all", ret=c("threshold","ppv","npv"))

# 计算 PPV 与 NPV 的差值
roc_coords$diff <- abs(roc_coords$ppv - roc_coords$npv)

# 找到差值最小的阈值
best_cutoff_row <- roc_coords[which.min(roc_coords$diff), ]
best_cutoff <- best_cutoff_row$threshold

cat("训练集 PPV=NPV 最优 cutoff:", round(best_cutoff,3), "\n")
cat("此时 PPV:", round(best_cutoff_row$ppv,3), 
    "NPV:", round(best_cutoff_row$npv,3), "\n")

# -------------------------------
# 1️⃣ 训练集 PPV=NPV 最优 cutoff
# -------------------------------
train_probs <- predict(final_rf_model, X_train[, final_variables, drop=FALSE], type="prob")[,2]

roc_train <- roc(y_train, train_probs)
roc_coords <- coords(roc_train, x = "all", ret=c("threshold","sensitivity","specificity","ppv","npv"))
roc_coords$diff <- abs(roc_coords$ppv - roc_coords$npv)

best_cutoff_row <- roc_coords[which.min(roc_coords$diff), ]
best_cutoff <- best_cutoff_row$threshold

cat("训练集 PPV=NPV 最优 cutoff:", round(best_cutoff,3), "\n")
cat("此时 PPV:", round(best_cutoff_row$ppv,3), 
    "NPV:", round(best_cutoff_row$npv,3), "\n")

# -------------------------------
# 2️⃣ 完整指标计算函数（新增 Prevalence & AUC & AUC95CI）
# -------------------------------
calculate_full_metrics <- function(y_true, probs, cutoff) {
  # 分类结果
  pred_class <- ifelse(probs >= cutoff, levels(y_true)[2], levels(y_true)[1])
  pred_class <- factor(pred_class, levels=levels(y_true))
  
  # 混淆矩阵
  cm <- confusionMatrix(pred_class, y_true, positive=levels(y_true)[2])
  TP <- cm$table[2,2]; TN <- cm$table[1,1]; FP <- cm$table[1,2]; FN <- cm$table[2,1]
  
  # 常规指标
  sensitivity <- as.numeric(cm$byClass["Sensitivity"])
  specificity <- as.numeric(cm$byClass["Specificity"])
  accuracy <- as.numeric(cm$overall["Accuracy"])
  ppv <- as.numeric(cm$byClass["Pos Pred Value"])
  npv <- as.numeric(cm$byClass["Neg Pred Value"])
  kappa <- as.numeric(cm$overall["Kappa"])
  
  # 平衡准确率 & Youden指数
  balanced_acc <- (sensitivity + specificity)/2
  youden_index <- sensitivity + specificity -1
  
  # F1 / F2
  F1 <- ifelse((ppv + sensitivity)==0, 0, 2*ppv*sensitivity/(ppv + sensitivity))
  F2 <- ifelse((ppv + sensitivity)==0, 0, (1 + 2^2)*ppv*sensitivity/(2^2*ppv + sensitivity))
  
  # MCC
  MCC <- ifelse((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)==0, NA, (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
  
  # PR AUC
  if(!requireNamespace("PRROC", quietly = TRUE)) install.packages("PRROC")
  library(PRROC)
  pr <- pr.curve(scores.class0 = probs[y_true==levels(y_true)[2]], 
                 scores.class1 = probs[y_true==levels(y_true)[1]], curve=FALSE)
  pr_auc <- pr$auc.integral
  
  # Prevalence
  prevalence <- mean(y_true == levels(y_true)[2])
  
  # ROC AUC + 95% CI
  roc_obj <- roc(y_true, probs)
  auc_val <- as.numeric(auc(roc_obj))
  auc_ci <- ci.auc(roc_obj, method="bootstrap")
  
  data.frame(
    阈值 = cutoff,
    灵敏度 = round(sensitivity,3),
    特异度 = round(specificity,3),
    准确率 = round(accuracy,3),
    平衡准确率 = round(balanced_acc,3),
    Youden指数 = round(youden_index,3),
    PPV = round(ppv,3),
    NPV = round(npv,3),
    F1 = round(F1,3),
    F2 = round(F2,3),
    MCC = round(MCC,3),
    PR_AUC = round(pr_auc,3),
    Kappa = round(kappa,3),
    Prevalence = round(prevalence,3),
    AUC = round(auc_val,3),
    AUC_95CI_lower = round(auc_ci[1],3),
    AUC_95CI_upper = round(auc_ci[3],3)
  )
}

# -------------------------------
# 3️⃣ 训练集 & 验证集指标（PPV=NPV cutoff）
# -------------------------------
test_probs <- predict(final_rf_model, X_test[, final_variables, drop=FALSE], type="prob")[,2]

train_metrics <- calculate_full_metrics(y_train, train_probs, best_cutoff)
train_metrics$数据集 <- "训练集"

test_metrics <- calculate_full_metrics(y_test, test_probs, best_cutoff)
test_metrics$数据集 <- "验证集"

all_metrics <- rbind(train_metrics, test_metrics)
all_metrics <- all_metrics[, c("数据集", names(train_metrics)[1:17])]  # 包含 AUC CI

# -------------------------------
# 4️⃣ 导出 CSV
# -------------------------------
write.csv(all_metrics, output_table_file("新再传播模型3.26评估指标_PPV_NPV_cutoff_AUC95CI.csv"), 
          row.names=FALSE, fileEncoding="UTF-8")
cat("\n✅ 已导出 CSV: 777训练验证集完整评估指标_PPV_NPV_cutoff_AUC95CI.csv\n")

###最优约登指数
# -------------------------------
# 1️⃣ 训练集 ROC & 最佳约登指数 cutoff
# -------------------------------
train_probs <- predict(final_rf_model, X_train[, final_variables, drop=FALSE], type="prob")[,2]

roc_train <- roc(y_train, train_probs)

# 获取所有阈值下的指标
roc_coords <- coords(roc_train, x = "all", ret=c("threshold","sensitivity","specificity","ppv","npv"))

# 计算约登指数（Youden index = sensitivity + specificity - 1）
roc_coords$youden_index <- roc_coords$sensitivity + roc_coords$specificity - 1

# 找到约登指数最大的阈值
best_cutoff_row <- roc_coords[which.max(roc_coords$youden_index), ]
best_cutoff <- best_cutoff_row$threshold

cat("训练集 最佳约登指数 cutoff:", round(best_cutoff,3), "\n")
cat("此时约登指数:", round(best_cutoff_row$youden_index,3), "\n")
cat("此时灵敏度:", round(best_cutoff_row$sensitivity,3), 
    "特异度:", round(best_cutoff_row$specificity,3), "\n")

# -------------------------------
# 2️⃣ 完整指标计算函数（新增 Prevalence & AUC & AUC95CI）
# -------------------------------
calculate_full_metrics <- function(y_true, probs, cutoff) {
  # 分类结果
  pred_class <- ifelse(probs >= cutoff, levels(y_true)[2], levels(y_true)[1])
  pred_class <- factor(pred_class, levels=levels(y_true))
  
  # 混淆矩阵
  cm <- confusionMatrix(pred_class, y_true, positive=levels(y_true)[2])
  TP <- cm$table[2,2]; TN <- cm$table[1,1]; FP <- cm$table[1,2]; FN <- cm$table[2,1]
  
  # 常规指标
  sensitivity <- as.numeric(cm$byClass["Sensitivity"])
  specificity <- as.numeric(cm$byClass["Specificity"])
  accuracy <- as.numeric(cm$overall["Accuracy"])
  ppv <- as.numeric(cm$byClass["Pos Pred Value"])
  npv <- as.numeric(cm$byClass["Neg Pred Value"])
  kappa <- as.numeric(cm$overall["Kappa"])
  
  # 平衡准确率 & Youden指数
  balanced_acc <- (sensitivity + specificity)/2
  youden_index <- sensitivity + specificity -1
  
  # F1 / F2
  F1 <- ifelse((ppv + sensitivity)==0, 0, 2*ppv*sensitivity/(ppv + sensitivity))
  F2 <- ifelse((ppv + sensitivity)==0, 0, (1 + 2^2)*ppv*sensitivity/(2^2*ppv + sensitivity))
  
  # MCC
  MCC <- ifelse((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)==0, NA, (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
  
  # PR AUC
  if(!requireNamespace("PRROC", quietly = TRUE)) install.packages("PRROC")
  library(PRROC)
  pr <- pr.curve(scores.class0 = probs[y_true==levels(y_true)[2]], 
                 scores.class1 = probs[y_true==levels(y_true)[1]], curve=FALSE)
  pr_auc <- pr$auc.integral
  
  # Prevalence
  prevalence <- mean(y_true == levels(y_true)[2])
  
  # ROC AUC + 95% CI
  roc_obj <- roc(y_true, probs)
  auc_val <- as.numeric(auc(roc_obj))
  auc_ci <- ci.auc(roc_obj, method="bootstrap")
  
  data.frame(
    阈值 = cutoff,
    灵敏度 = round(sensitivity,3),
    特异度 = round(specificity,3),
    准确率 = round(accuracy,3),
    平衡准确率 = round(balanced_acc,3),
    Youden指数 = round(youden_index,3),
    PPV = round(ppv,3),
    NPV = round(npv,3),
    F1 = round(F1,3),
    F2 = round(F2,3),
    MCC = round(MCC,3),
    PR_AUC = round(pr_auc,3),
    Kappa = round(kappa,3),
    Prevalence = round(prevalence,3),
    AUC = round(auc_val,3),
    AUC_95CI_lower = round(auc_ci[1],3),
    AUC_95CI_upper = round(auc_ci[3],3)
  )
}

# -------------------------------
# 3️⃣ 训练集 & 验证集指标（最佳约登指数 cutoff）
# -------------------------------
test_probs <- predict(final_rf_model, X_test[, final_variables, drop=FALSE], type="prob")[,2]

train_metrics <- calculate_full_metrics(y_train, train_probs, best_cutoff)
train_metrics$数据集 <- "训练集"

test_metrics <- calculate_full_metrics(y_test, test_probs, best_cutoff)
test_metrics$数据集 <- "验证集"

all_metrics <- rbind(train_metrics, test_metrics)
all_metrics <- all_metrics[, c("数据集", names(train_metrics)[1:17])]  # 包含 AUC CI

# -------------------------------
# 4️⃣ 导出 CSV
# -------------------------------
write.csv(all_metrics, output_table_file("新死亡模型3.26评估指标_Youden指数_cutoff_AUC95CI.csv"), 
          row.names=FALSE, fileEncoding="UTF-8")
cat("\n✅ 已导出 CSV: 新再传播模型3.26评估指标_Youden指数_cutoff_AUC95CI.csv\n")

# ==================== 新增代码：导出每个样本的预测值 ====================

# 获取训练集和测试集的预测概率
train_predictions_prob <- predict(final_rf_model, X_train[, final_variables, drop = FALSE], type = "prob")
test_predictions_prob <- predict(final_rf_model, X_test[, final_variables, drop = FALSE], type = "prob")

# 获取训练集和测试集的预测类别
train_predictions_class <- predict(final_rf_model, X_train[, final_variables, drop = FALSE])
test_predictions_class <- predict(final_rf_model, X_test[, final_variables, drop = FALSE])

# 创建训练集预测结果数据框
train_results <- data.frame(
  行号 = as.numeric(rownames(X_train)),
  预测类别 = train_predictions_class,
  预测概率_类别0 = train_predictions_prob[, 1],
  预测概率_类别1 = train_predictions_prob[, 2],
  数据集 = "训练集"
)

# 创建测试集预测结果数据框
test_results <- data.frame(
  行号 = as.numeric(rownames(X_test)),
  预测类别 = test_predictions_class,
  预测概率_类别0 = test_predictions_prob[, 1],
  预测概率_类别1 = test_predictions_prob[, 2],
  数据集 = "测试集"
)

# 合并训练集和测试集的预测结果
all_predictions <- rbind(train_results, test_results)

# 获取原始数据的行号
original_data <- cbind(
  行号 = as.numeric(rownames(mydata)),
  mydata
)

# 将预测结果与原始数据合并
final_combined_data <- merge(original_data, all_predictions, by = "行号", all.x = TRUE)

# 重新排列列的顺序，使预测相关列在最后
col_order <- c("行号", "数据集", names(mydata), 
               "预测类别", "预测概率_类别0", "预测概率_类别1")
final_combined_data <- final_combined_data[, col_order]

# 导出合并后的数据
write.csv(final_combined_data, output_table_file("新再传播3.26样本预测结果.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("\n样本预测结果已保存到 '样本预测结果.csv'\n")

# 显示前几行数据作为预览
cat("\n合并数据预览 (前10行):\n")
print(head(final_combined_data, 10))

# ==================== 修正绘图代码 ====================

# 创建更美观的重要性图
library(ggplot2)
library(tidyr)

# 使用正确的重要性数据对象 - 应该是 importance_data 而不是 importance_df
# 将数据转换为长格式以便绘图
if ("MeanDecreaseAccuracy" %in% colnames(importance_data) && "MeanDecreaseGini" %in% colnames(importance_data)) {
  importance_long <- importance_data %>%
    pivot_longer(cols = c(MeanDecreaseAccuracy, MeanDecreaseGini), 
                 names_to = "Importance_Type", 
                 values_to = "Value")
  
  # 绘制变量重要性图
  ggplot(importance_long, aes(x = reorder(Variable, Value), y = Value, fill = Importance_Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    labs(title = "Final Model Variable Importance",
         x = "Variables",
         y = "Importance Value",
         fill = "Importance Type") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
} else {
  cat("无法绘制重要性图：缺少MeanDecreaseAccuracy或MeanDecreaseGini列\n")
}

# 分别绘制MeanDecreaseAccuracy和MeanDecreaseGini
if ("MeanDecreaseAccuracy" %in% colnames(importance_data)) {
  # MeanDecreaseAccuracy图
  ggplot(importance_data, aes(x = reorder(Variable, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "Variable Importance - Mean Decrease Accuracy",
         x = "Variables",
         y = "Mean Decrease Accuracy") +
    theme_minimal()
}

if ("MeanDecreaseGini" %in% colnames(importance_data)) {
  # MeanDecreaseGini图
  ggplot(importance_data, aes(x = reorder(Variable, MeanDecreaseGini), y = MeanDecreaseGini)) +
    geom_bar(stat = "identity", fill = "darkorange") +
    coord_flip() +
    labs(title = "Variable Importance - Mean Decrease Gini",
         x = "Variables",
         y = "Mean Decrease Gini") +
    theme_minimal()
}

# 创建组合图
library(gridExtra)

p1 <- ggplot(importance_data, aes(x = reorder(Variable, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Mean Decrease Accuracy",
       x = "",
       y = "Mean Decrease Accuracy") +
  theme_minimal()

p2 <- ggplot(importance_data, aes(x = reorder(Variable, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  coord_flip() +
  labs(title = "Mean Decrease Gini",
       x = "",
       y = "Mean Decrease Gini") +
  theme_minimal()

# 并排显示两个图
grid.arrange(p1, p2, ncol = 2)



###3相交曲线绘图代码
# -------------------------------
# [FINAL ENGLISH] Sensitivity = Specificity Plot (0~1 Cutoff)
# -------------------------------
library(ggplot2)
library(reshape2)

# 1. Calculate accuracy for all thresholds
thresholds <- roc_coords$threshold
accuracy_list <- c()

for (c in thresholds) {
  pred_class <- ifelse(train_probs >= c, levels(y_train)[2], levels(y_train)[1])
  pred_class <- factor(pred_class, levels = levels(y_train))
  cm <- confusionMatrix(pred_class, y_train, positive = levels(y_train)[2])
  acc <- as.numeric(cm$overall["Accuracy"])
  accuracy_list <- c(accuracy_list, acc)
}

# 2. Build plot data
plot_data <- data.frame(
  cutoff = rep(roc_coords$threshold, 3),
  value = c(roc_coords$sensitivity, roc_coords$specificity, accuracy_list),
  metric = c(rep("Sensitivity", nrow(roc_coords)),
             rep("Specificity", nrow(roc_coords)),
             rep("Accuracy", nrow(roc_coords)))
)

# Keep only 0 ≤ cutoff ≤ 1
plot_data <- plot_data[plot_data$cutoff >= 0 & plot_data$cutoff <= 1, ]

# 3. Find TRUE intersection (Sensitivity = Specificity)
plot_data_wide <- dcast(plot_data, cutoff ~ metric, value.var = "value")
plot_data_wide$diff <- abs(plot_data_wide$Sensitivity - plot_data_wide$Specificity)
best_idx <- which.min(plot_data_wide$diff)

se_sp_point <- data.frame(
  cutoff = plot_data_wide$cutoff[best_idx],
  value = plot_data_wide$Sensitivity[best_idx]
)

# 4. Plot (PERFECT ALIGNMENT)
p <- ggplot(plot_data) +
  geom_line(aes(x = cutoff, y = value, color = metric), linewidth = 0.7) +
  
  geom_point(data = se_sp_point, 
             aes(x = cutoff, y = value), 
             color = "red", size = 4) +
  
  geom_segment(data = se_sp_point,
               aes(x = cutoff, y = value, xend = cutoff, yend = 0),
               linetype = "dashed", color = "black", linewidth = 0.5) +
  
  annotate("text",
           x = se_sp_point$cutoff,
           y = se_sp_point$value + 0.04,
           label = paste0("Se=Sp=", round(se_sp_point$value, 3)),
           color = "red", fontface = "bold", size = 4.5) +
  
  scale_color_manual(values = c(
    "Sensitivity" = "#2E86AB",
    "Specificity" = "#F18F01",
    "Accuracy" = "#00A600"
  )) +
  
  labs(
    x = "Cutoff Value",
    y = "Value",
    color = "Metric"
  ) +
  
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

# Show in plot pane
print(p)

# Save PDF
ggsave(output_figure_file("Train_Sensitivity_Specificity_Accuracy.pdf"),
       plot = p, device = "pdf", width = 10, height = 6)

cat("✅ PDF saved: Train_Sensitivity_Specificity_Accuracy.pdf\n")


# -------------------------------
# [FINAL ENGLISH] PPV + NPV + Accuracy Plot (0~1 Cutoff)
# -------------------------------
library(ggplot2)

# 1. Calculate accuracy for all thresholds
thresholds <- roc_coords$threshold
accuracy_list <- c()

for (c in thresholds) {
  pred_class <- ifelse(train_probs >= c, levels(y_train)[2], levels(y_train)[1])
  pred_class <- factor(pred_class, levels = levels(y_train))
  cm <- confusionMatrix(pred_class, y_train, positive = levels(y_train)[2])
  acc <- as.numeric(cm$overall["Accuracy"])
  accuracy_list <- c(accuracy_list, acc)
}

# 2. Create long format data
plot_data <- data.frame(
  cutoff = rep(roc_coords$threshold, 3),
  value = c(roc_coords$ppv, roc_coords$npv, accuracy_list),
  metric = c(rep("PPV", nrow(roc_coords)),
             rep("NPV", nrow(roc_coords)),
             rep("Accuracy", nrow(roc_coords)))
)

# Keep only 0 ≤ cutoff ≤ 1
plot_data <- plot_data[plot_data$cutoff >= 0 & plot_data$cutoff <= 1, ]

# 3. PPV = NPV optimal point
best_point <- data.frame(
  cutoff = best_cutoff_row$threshold,
  value = best_cutoff_row$ppv
)

# 4. Plot (clean, English version)
p <- ggplot(plot_data) +
  geom_line(aes(x = cutoff, y = value, color = metric), linewidth = 0.7) +
  
  # Intersection point
  geom_point(data = best_point, aes(x = cutoff, y = value), 
             color = "red", size = 3.5) +
  
  # Vertical dashed line to x-axis
  geom_segment(data = best_point,
               aes(x = cutoff, y = value, xend = cutoff, yend = 0),
               linetype = "dashed", color = "black", linewidth = 0.5) +
  
  # Label
  annotate("text",
           x = best_point$cutoff,
           y = best_point$value + 0.05,
           label = paste0("PPV=NPV=", round(best_point$value,3)),
           color = "red", fontface = "bold", size = 4.5) +
  
  # Color
  scale_color_manual(values = c(
    "PPV" = "#2E86AB",
    "NPV" = "#F18F01",
    "Accuracy" = "#00A600"
  )) +
  
  # Labels (ALL ENGLISH)
  labs(x = "Cutoff Value",
       y = "Value",
       color = "Metric") +
  
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

# Show plot
print(p)

# Save PDF (English filename)
ggsave(output_figure_file("死亡Train_PPV_NPV_Accuracy.pdf"),
       plot = p, device = "pdf", width = 10, height = 6)

cat("✅ PDF saved: Train_PPV_NPV_Accuracy.pdf\n")


# -------------------------------
# [FINAL ENGLISH] Youden Index + Accuracy Plot (0~1 Cutoff)
# -------------------------------
library(ggplot2)

# 1. Calculate accuracy for all thresholds
thresholds <- roc_coords$threshold
accuracy_list <- c()

for (c in thresholds) {
  pred_class <- ifelse(train_probs >= c, levels(y_train)[2], levels(y_train)[1])
  pred_class <- factor(pred_class, levels = levels(y_train))
  cm <- confusionMatrix(pred_class, y_train, positive = levels(y_train)[2])
  acc <- as.numeric(cm$overall["Accuracy"])
  accuracy_list <- c(accuracy_list, acc)
}

# 2. Create plot data
plot_data <- data.frame(
  cutoff = rep(roc_coords$threshold, 2),
  value = c(roc_coords$youden_index, accuracy_list),
  metric = c(rep("Youden Index", nrow(roc_coords)),
             rep("Accuracy", nrow(roc_coords)))
)

# Keep only 0 ≤ cutoff ≤ 1
plot_data <- plot_data[plot_data$cutoff >= 0 & plot_data$cutoff <= 1, ]

# 3. Best Youden point
best_point <- data.frame(
  cutoff = best_cutoff_row$threshold,
  youden = best_cutoff_row$youden_index
)

# 4. Plot (Full English)
p <- ggplot(plot_data) +
  geom_line(aes(x = cutoff, y = value, color = metric), linewidth = 0.7) +
  geom_point(data = best_point, aes(x = cutoff, y = youden), color = "red", size = 3.5) +
  geom_segment(data = best_point,
               aes(x = cutoff, y = youden, xend = cutoff, yend = 0),
               linetype = "dashed", color = "black", linewidth = 0.5) +
  annotate("text",
           x = best_point$cutoff,
           y = best_point$youden + 0.05,
           label = paste0("Best Youden=", round(best_point$youden,3)),
           color = "red", fontface = "bold", size = 4.5) +
  scale_color_manual(values = c("Youden Index"="#2E86AB", "Accuracy"="#00A600")) +
  labs(x="Cutoff Value", 
       y="Value", 
       color="Metric") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw(base_size=13) +
  theme(legend.position="top")

# Show plot
print(p)

# Save PDF
ggsave(output_figure_file("死亡Train_Youden_Accuracy.pdf"), plot=p, device="pdf", width=10, height=6)
cat("✅ PDF saved: Train_Youden_Accuracy.pdf\n")

###灵敏度等于特异度
# -------------------------------
library(ggplot2)
library(reshape2)

# 1. Calculate accuracy for all thresholds
thresholds <- roc_coords$threshold
accuracy_list <- c()

for (c in thresholds) {
  pred_class <- ifelse(train_probs >= c, levels(y_train)[2], levels(y_train)[1])
  pred_class <- factor(pred_class, levels = levels(y_train))
  cm <- caret::confusionMatrix(pred_class, y_train, positive = levels(y_train)[2])
  acc <- as.numeric(cm$overall["Accuracy"])
  accuracy_list <- c(accuracy_list, acc)
}

# 2. Build plot data
plot_data <- data.frame(
  cutoff = rep(roc_coords$threshold, 3),
  value = c(roc_coords$sensitivity, roc_coords$specificity, accuracy_list),
  metric = c(rep("Sensitivity", nrow(roc_coords)),
             rep("Specificity", nrow(roc_coords)),
             rep("Accuracy", nrow(roc_coords)))
)

# Keep only 0 ≤ cutoff ≤ 1
plot_data <- plot_data[plot_data$cutoff >= 0 & plot_data$cutoff <= 1, ]

# -------------------------------
# 核心修改：生成平滑曲线所需的密集点
# -------------------------------
smooth_n <- 500  # 平滑密度，越大越平滑
new_cutoff <- seq(0, 1, length.out = smooth_n)

# 对三个指标分别做平滑拟合
smooth_functions <- lapply(split(plot_data, plot_data$metric), function(df) {
  loess(value ~ cutoff, data = df, span = 0.15)  # span越小越贴近原始数据，越大越平滑
})

# 生成平滑后的数据
plot_data_smooth <- expand.grid(
  cutoff = new_cutoff,
  metric = c("Sensitivity", "Specificity", "Accuracy")
)
plot_data_smooth$value <- mapply(function(m, x) {
  predict(smooth_functions[[m]], newdata = data.frame(cutoff = x))
}, plot_data_smooth$metric, plot_data_smooth$cutoff)

# 限制 0~1 范围
plot_data_smooth$value <- pmax(pmin(plot_data_smooth$value, 1), 0)

# 3. Find TRUE intersection (Sensitivity = Specificity)
plot_data_wide <- dcast(plot_data, cutoff ~ metric, value.var = "value")
plot_data_wide$diff <- abs(plot_data_wide$Sensitivity - plot_data_wide$Specificity)
best_idx <- which.min(plot_data_wide$diff)

se_sp_point <- data.frame(
  cutoff = plot_data_wide$cutoff[best_idx],
  value = plot_data_wide$Sensitivity[best_idx]
)

# 4. 平滑曲线绘图
p <- ggplot() +
  # 平滑曲线（关键修改）
  geom_line(data = plot_data_smooth,
            aes(x = cutoff, y = value, color = metric),
            linewidth = 1.2) +
  
  geom_point(data = se_sp_point, 
             aes(x = cutoff, y = value), 
             color = "red", size = 4) +
  
  geom_segment(data = se_sp_point,
               aes(x = cutoff, y = value, xend = cutoff, yend = 0),
               linetype = "dashed", color = "black", linewidth = 0.5) +
  
  annotate("text",
           x = se_sp_point$cutoff,
           y = se_sp_point$value + 0.04,
           label = paste0("Se=Sp=", round(se_sp_point$value, 3)),
           color = "red", fontface = "bold", size = 4.5) +
  
  scale_color_manual(values = c(
    "Sensitivity" = "#2E86AB",
    "Specificity" = "#F18F01",
    "Accuracy" = "#00A600"
  )) +
  
  labs(
    x = "Cutoff Value",
    y = "Value",
    color = "Metric"
  ) +
  
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

print(p)

# Save PDF
ggsave(output_figure_file("0414再传播Train_Sensitivity_Specificity_Accuracy.pdf"),
       plot = p, device = "pdf", width = 10, height = 6)

cat("✅ PDF saved: Train_Sensitivity_Specificity_Accuracy.pdf\n")

###平滑NPV=PPV
# -------------------------------
# [FINAL ENGLISH] PPV + NPV + Accuracy Plot (SMOOTH CURVE)
# -------------------------------
library(ggplot2)
library(reshape2)

# 1. Calculate accuracy for all thresholds
thresholds <- roc_coords$threshold
accuracy_list <- c()

for (c in thresholds) {
  pred_class <- ifelse(train_probs >= c, levels(y_train)[2], levels(y_train)[1])
  pred_class <- factor(pred_class, levels = levels(y_train))
  cm <- caret::confusionMatrix(pred_class, y_train, positive = levels(y_train)[2])
  acc <- as.numeric(cm$overall["Accuracy"])
  accuracy_list <- c(accuracy_list, acc)
}

# 2. Create long format data
plot_data <- data.frame(
  cutoff = rep(roc_coords$threshold, 3),
  value = c(roc_coords$ppv, roc_coords$npv, accuracy_list),
  metric = c(rep("PPV", nrow(roc_coords)),
             rep("NPV", nrow(roc_coords)),
             rep("Accuracy", nrow(roc_coords)))
)

# Keep only 0 ≤ cutoff ≤ 1
plot_data <- plot_data[plot_data$cutoff >= 0 & plot_data$cutoff <= 1, ]

# ==============================
# 核心：生成平滑曲线（关键修改）
# ==============================
smooth_n <- 500    # 平滑密度，越大越顺滑
new_cutoff <- seq(0, 1, length.out = smooth_n)

# 对 PPV / NPV / Accuracy 分别做平滑拟合
smooth_fit <- lapply(split(plot_data, plot_data$metric), function(df) {
  loess(value ~ cutoff, data = df, span = 0.12)  # span 越小越贴原始，越大越平滑
})

# 生成平滑后的数据
plot_data_smooth <- expand.grid(
  cutoff = new_cutoff,
  metric = c("PPV", "NPV", "Accuracy")
)
plot_data_smooth$value <- mapply(function(m, x) {
  predict(smooth_fit[[m]], newdata = data.frame(cutoff = x))
}, plot_data_smooth$metric, plot_data_smooth$cutoff)

# 把数值限制在 0~1 之间
plot_data_smooth$value <- pmax(pmin(plot_data_smooth$value, 1), 0)

# 3. PPV = NPV optimal point
best_point <- data.frame(
  cutoff = best_cutoff_row$threshold,
  value = best_cutoff_row$ppv
)

# 4. Plot (clean, English version + SMOOTH)
p <- ggplot() +
  # 平滑曲线（替换原来的 geom_line）
  geom_line(data = plot_data_smooth,
            aes(x = cutoff, y = value, color = metric),
            linewidth = 1.0) +
  
  # Intersection point
  geom_point(data = best_point, aes(x = cutoff, y = value), 
             color = "red", size = 3.5) +
  
  # Vertical dashed line to x-axis
  geom_segment(data = best_point,
               aes(x = cutoff, y = value, xend = cutoff, yend = 0),
               linetype = "dashed", color = "black", linewidth = 0.5) +
  
  # Label
  annotate("text",
           x = best_point$cutoff,
           y = best_point$value + 0.05,
           label = paste0("PPV=NPV=", round(best_point$value,3)),
           color = "red", fontface = "bold", size = 4.5) +
  
  # Color
  scale_color_manual(values = c(
    "PPV" = "#2E86AB",
    "NPV" = "#F18F01",
    "Accuracy" = "#00A600"
  )) +
  
  # Labels (ALL ENGLISH)
  labs(x = "Cutoff Value",
       y = "Value",
       color = "Metric") +
  
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

# Show plot
print(p)

ggsave(output_figure_file("0414再传播Train_PPV_NPV_Accuracy.pdf"),
       plot = p, device = "pdf", width = 10, height = 6)

cat("✅ PDF saved: Train_PPV_NPV_Accuracy.pdf\n")

###Youden指数
# -------------------------------
# [FINAL ENGLISH] Youden Index + Accuracy Plot (SMOOTH CURVE)
# -------------------------------
library(ggplot2)
library(reshape2)

# 1. Calculate accuracy for all thresholds
thresholds <- roc_coords$threshold
accuracy_list <- c()

for (c in thresholds) {
  pred_class <- ifelse(train_probs >= c, levels(y_train)[2], levels(y_train)[1])
  pred_class <- factor(pred_class, levels = levels(y_train))
  cm <- caret::confusionMatrix(pred_class, y_train, positive = levels(y_train)[2])
  acc <- as.numeric(cm$overall["Accuracy"])
  accuracy_list <- c(accuracy_list, acc)
}

# 2. Create plot data
plot_data <- data.frame(
  cutoff = rep(roc_coords$threshold, 2),
  value = c(roc_coords$youden_index, accuracy_list),
  metric = c(rep("Youden Index", nrow(roc_coords)),
             rep("Accuracy", nrow(roc_coords)))
)

# Keep only 0 ≤ cutoff ≤ 1
plot_data <- plot_data[plot_data$cutoff >= 0 & plot_data$cutoff <= 1, ]

# ==============================
# 核心：平滑曲线处理（关键修改）
# ==============================
smooth_n <- 500    # 平滑密度，越大越顺滑
new_cutoff <- seq(0, 1, length.out = smooth_n)

# 对 Youden Index 和 Accuracy 分别做平滑
smooth_fit <- lapply(split(plot_data, plot_data$metric), function(df) {
  loess(value ~ cutoff, data = df, span = 0.12)  # 越小越贴原始，越大越平滑
})

# 生成平滑后的数据
plot_data_smooth <- expand.grid(
  cutoff = new_cutoff,
  metric = c("Youden Index", "Accuracy")
)
plot_data_smooth$value <- mapply(function(m, x) {
  predict(smooth_fit[[m]], newdata = data.frame(cutoff = x))
}, plot_data_smooth$metric, plot_data_smooth$cutoff)

# 数值限制在 0~1
plot_data_smooth$value <- pmax(pmin(plot_data_smooth$value, 1), 0)

# 3. Best Youden point
best_point <- data.frame(
  cutoff = best_cutoff_row$threshold,
  youden = best_cutoff_row$youden_index
)

# 4. Plot (Full English + SMOOTH)
p <- ggplot() +
  # 平滑曲线（替换原来的 geom_line）
  geom_line(data = plot_data_smooth,
            aes(x = cutoff, y = value, color = metric),
            linewidth = 1.0) +
  
  geom_point(data = best_point, aes(x = cutoff, y = youden), color = "red", size = 3.5) +
  geom_segment(data = best_point,
               aes(x = cutoff, y = youden, xend = cutoff, yend = 0),
               linetype = "dashed", color = "black", linewidth = 0.5) +
  annotate("text",
           x = best_point$cutoff,
           y = best_point$youden + 0.05,
           label = paste0("", round(best_point$youden,3)),
           color = "red", fontface = "bold", size = 4.5) +
  scale_color_manual(values = c("Youden Index"="#2E86AB", "Accuracy"="#00A600")) +
  labs(x="Cutoff Value", 
       y="Value", 
       color="Metric") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw(base_size=13) +
  theme(legend.position="top")

# Show plot
print(p)

# Save PDF
ggsave(output_figure_file("0414再传播Train_Youden_Accuracy.pdf"), plot=p, device="pdf", width=10, height=6)
cat("✅ PDF saved: Train_Youden_Accuracy.pdf\n")
