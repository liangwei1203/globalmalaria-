# ============================================================
# File: 04_onward_transmission_logistic_regression.R
# Original file: 3.18更新的logistic代码.R
# Purpose: Conventional logistic regression for onward transmission risk
# Study: Imported malaria analysis
#
# Inputs:
# - introduced260318_imputed.csv
#
# Outputs:
# - 0318基线表table1_output.csv
# - 3.19_再传播Multiv_logistic_rms格式.csv
#
# Notes:
# - Personal working-directory paths have been removed.
# - Run from the repository root using source("R/04_onward_transmission_logistic_regression.R").
# - Original unmodified script is retained in original_scripts/.
# ============================================================

source(file.path("R", "00_setup.R"))

# setwd() removed; project paths are defined in R/00_setup.R
mydata <- read.csv(input_file("introduced260318_imputed.csv"))
# View(mydata)  # Interactive inspection removed for reproducibility
library(dplyr)
library(tableone)         # 基线表
library(car)              # VIF
library(pROC)             # ROC / AUC
library(PRROC)            # PR 曲线（可选）
library(rms)              # 列线图
library(ResourceSelection)# Hosmer-Lemeshow 拟合优度
library(caret)            # preProcess 归一化

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


# 因子转换 - 只对分类变量进行因子转换
mydata$Outcome = as.factor(mydata$Outcome)
mydata$Region_of_acquisition = as.factor(mydata$Region_of_acquisition)
mydata$Country_of_destination = as.factor(mydata$Country_of_destination)
mydata$Plasmodium_species = as.factor(mydata$Plasmodium_species)
mydata$Occupation = as.factor(mydata$Occupation)
mydata$Groups_of_economies = as.factor(mydata$Groups_of_economies)

# 指定需要筛选的变量 - 明确区分分类变量和数值变量
categorical_vars <- c("Region_of_acquisition", "Country_of_destination", 
                      "Plasmodium_species", "Occupation", "Groups_of_economies")

continuous_vars <- c("Temp", "Precip", "Mosquito_borne_species")

# 合并所有变量
selected_columns <- c(categorical_vars, continuous_vars)

# 选择指定的变量
mydata <- mydata %>% select(all_of(selected_columns), Outcome)

# ----------------------------
# 基线表（Table 1）
# ----------------------------

# 设置分层变量（一般是 Outcome）
strata_var <- "Outcome"

# 创建 Table 1 对象
# 注意：factorVars 只包含分类变量，数值变量会自动按照均值±标准差展示
table1 <- CreateTableOne(vars = selected_columns,
                         strata = strata_var,
                         data = mydata,
                         factorVars = categorical_vars)  # 只将分类变量指定为因子

# 显示结果
print(table1, showAllLevels = TRUE, quote = TRUE, noSpaces = TRUE)

# 可选：以更美观的方式显示
print(table1, showAllLevels = TRUE, formatOptions = list(big.mark = ","))

# 可选：导出为CSV文件
table1_matrix <- print(table1, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table1_matrix, file = output_table_file("0318基线表table1_output.csv"))

# ----------------------------
# 筛选 p < 0.05 的显著变量
# ----------------------------
selected_significant_names <- c(
  "Region_of_acquisition",
  "Country_of_destination",
  "Plasmodium_species",
  "Occupation",
  "Temp"
)

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

####共线性检测
str(mydata[selected_significant_names])
vif_model <- glm(
  Outcome ~ Region_of_acquisition+Country_of_destination +
    Plasmodium_species +
    Occupation +
    Temp,
  data = mydata,
  family = binomial
)

vif_raw <- car::vif(vif_model)

# car 正常返回的是 matrix（变量层级）
vif_df <- data.frame(
  Variable = rownames(vif_raw),
  GVIF = vif_raw[, "GVIF"],
  Df = vif_raw[, "Df"],
  GVIF_adj = vif_raw[, "GVIF"]^(1 / (2 * vif_raw[, "Df"]))
)

print(vif_df)

## -------- 1. 明确要设定参考值的变量 --------
vars_ref_1 <- c("Country_of_destination")
vars_ref_2 <- c("Plasmodium_species")
vars_ref_3 <- c("Occupation")
vars_ref_6 <- c("Region_of_acquisition")

## -------- 2. 设置参考水平（自动判断是否为 factor） --------
# 设置参考 = "1"
for (var in vars_ref_1) {
  if (is.factor(mydata[[var]])) {
    mydata[[var]] <- relevel(mydata[[var]], ref = "1")
  } else {
    message(sprintf("变量 %s 不是因子，跳过。", var))
  }
}

# 设置参考 = "2"
for (var in vars_ref_2) {
  if (is.factor(mydata[[var]])) {
    mydata[[var]] <- relevel(mydata[[var]], ref = "2")
  } else {
    message(sprintf("变量 %s 不是因子，跳过。", var))
  }
}

# 设置参考 = "3"
for (var in vars_ref_3) {
  if (is.factor(mydata[[var]])) {
    mydata[[var]] <- relevel(mydata[[var]], ref = "3")
  } else {
    message(sprintf("变量 %s 不是因子，跳过。", var))
  }
}

# 设置参考 = "6"
for (var in vars_ref_6) {
  if (is.factor(mydata[[var]])) {
    mydata[[var]] <- relevel(mydata[[var]], ref = "6")
  } else {
    message(sprintf("变量 %s 不是因子，跳过。", var))
  }
}

# 3) 拟合多元 logistic 模型（用所有显著变量）
multiv_formula <- as.formula(paste("Outcome ~", paste(selected_significant_names, collapse = " + ")))
multiv_model <- glm(multiv_formula, data = mydata, family = binomial)


# ---- 从 glm(Logistic) 模型中提取结果，并按 rms 风格输出 ----

# 1) 提取系数与方差
coefs <- coef(multiv_model)                 # β
se    <- sqrt(diag(vcov(multiv_model)))    # 标准误（来自协方差矩阵）

# 2) 计算 Wald χ² 与 P 值
wald_chi2 <- (coefs / se)^2
p_val     <- 1 - pchisq(wald_chi2, df = 1)

# 3) 计算 OR 与 95% CI
OR      <- exp(coefs)
OR_low  <- exp(coefs - 1.96 * se)
OR_high <- exp(coefs + 1.96 * se)

OR_CI <- paste0(
  round(OR, 3), 
  " (", round(OR_low, 3), "-", round(OR_high, 3), ")"
)

# 4) 变量名称
vars <- names(coefs)

# 5) 构建结果表（与第二段格式一致）
result_table <- data.frame(
  自变量      = vars,
  β           = round(coefs, 3),
  Wald        = round(wald_chi2, 3),
  `OR (95% CI)` = OR_CI,
  P值         = signif(p_val, 3),
  stringsAsFactors = FALSE
)

# 6) 去掉截距
result_table <- subset(result_table, 自变量 != "(Intercept)")

# 7) 导出
write.csv(
  result_table,
  "3.19_再传播Multiv_logistic_rms格式.csv",
  row.names = FALSE
)

# 查看
print(result_table)

library(pROC)
library(caret)
library(ggplot2)

# ----------------------------
# 1. 得到预测概率（预测 Outcome = 1 的概率）
# ----------------------------
pred_prob <- predict(multiv_model, type = "response")

# ----------------------------
# 2. AUC 及 95% CI
# ----------------------------
roc_obj <- roc(mydata$Outcome, pred_prob, ci = TRUE)

auc_value <- auc(roc_obj)
auc_ci <- ci.auc(roc_obj)

cat("AUC =", auc_value, "\n")
cat("AUC 95% CI =", auc_ci[1], "-", auc_ci[3], "\n")

plot(roc_obj,
     col = "blue",
     lwd = 2,
     main = "ROC Curve for Multivariable Logistic Model")

text(0.6, 0.2,
     labels = paste0("AUC = ", round(auc_value, 3)),
     col = "red")

library(caret)


###默认阈值为0.5
# 1. 预测概率
pred_prob <- predict(multiv_model, type = "response")

# 2. 默认阈值 0.5 转为二分类
pred_class <- factor(ifelse(pred_prob >= 0.5, 1, 0), levels = c(0, 1))

# 3. 混淆矩阵
cm <- confusionMatrix(pred_class, mydata$Outcome, positive = "1")

# 4. 计算指标
sensitivity <- cm$byClass["Sensitivity"]
specificity <- cm$byClass["Specificity"]
ppv <- cm$byClass["Pos Pred Value"]
npv <- cm$byClass["Neg Pred Value"]
kappa <- cm$overall["Kappa"]    # <-- 新增的 Kappa 值

# 5. 输出
cat("Sensitivity =", sensitivity, "\n")
cat("Specificity =", specificity, "\n")
cat("PPV =", ppv, "\n")
cat("NPV =", npv, "\n")
cat("Kappa =", kappa, "\n")      # 输出 Kappa


##se=sp
library(pROC)
library(caret)

# 得到预测概率
pred_prob <- predict(multiv_model, type = "response")
true <- mydata$Outcome

# 使用 pROC 包计算 ROC
roc_obj <- roc(true, pred_prob, levels = c("0", "1"), direction = "<")

# 提取所有 cutoff、敏感性、特异性
roc_df <- data.frame(
  threshold = roc_obj$thresholds,
  sensitivity = roc_obj$sensitivities,
  specificity = roc_obj$specificities
)

# 找到 sensitivity 与 specificity 最接近的点
roc_df$diff <- abs(roc_df$sensitivity - roc_df$specificity)
best_row <- roc_df[which.min(roc_df$diff), ]

best_row

# 取最佳阈值
best_threshold <- best_row$threshold

# 使用该阈值进行分类
pred_class_best <- factor(ifelse(pred_prob >= best_threshold, "1", "0"), levels = c("0","1"))

# 生成混淆矩阵
cm_best <- confusionMatrix(pred_class_best, true, positive = "1")

# 提取指标
sensitivity_best  <- cm_best$byClass["Sensitivity"]
specificity_best  <- cm_best$byClass["Specificity"]
ppv_best          <- cm_best$byClass["Pos Pred Value"]
npv_best          <- cm_best$byClass["Neg Pred Value"]


cat("最佳阈值（敏感性≈特异性） =", best_threshold, "\n")
cat("Sensitivity =", sensitivity_best, "\n")
cat("Specificity =", specificity_best, "\n")
cat("PPV =", ppv_best, "\n")
cat("NPV =", npv_best, "\n")
kappa_best <- cm_best$overall["Kappa"]
cat("Kappa =", kappa_best, "\n")

###NPV=PPV
library(caret)

pred_prob <- predict(multiv_model, type = "response")
true <- mydata$Outcome

thresholds <- seq(0.01, 0.99, by = 0.001)   # 精细的阈值搜索
results <- data.frame()

for (t in thresholds) {
  
  pred_class <- factor(ifelse(pred_prob >= t, "1", "0"), levels = c("0","1"))
  
  cm <- confusionMatrix(pred_class, true, positive = "1")
  
  sens <- cm$byClass["Sensitivity"]
  spec <- cm$byClass["Specificity"]
  ppv  <- cm$byClass["Pos Pred Value"]
  npv  <- cm$byClass["Neg Pred Value"]
  kappa <- cm$overall["Kappa"]
  
  results <- rbind(results,
                   data.frame(
                     threshold = t,
                     sensitivity = sens,
                     specificity = spec,
                     PPV = ppv,
                     NPV = npv,
                     Kappa = kappa,
                     diff = abs(ppv - npv)   # <-- PPV 与 NPV 的差
                   ))
}

# 找到 PPV ≈ NPV 的阈值
best_row <- results[which.min(results$diff), ]
best_row

cat("最佳阈值（PPV = NPV） =", best_row$threshold, "\n")
cat("Sensitivity =", best_row$sensitivity, "\n")
cat("Specificity =", best_row$specificity, "\n")
cat("PPV =", best_row$PPV, "\n")
cat("NPV =", best_row$NPV, "\n")
cat("Kappa =", best_row$Kappa, "\n")




######绘制列线图
# 最终模型（用于列线图）
library(rms)

# 1. 指定变量（确保名称与你数据一致）
final_vars <- c("Plasmodium_species",
                "Occupation","Temp")

# 2. 构建公式
final_formula <- as.formula(
  paste("Outcome ~", paste(final_vars, collapse = " + "))
)

# 3. 设置 rms 所需的数据结构
dd <- datadist(mydata)
options(datadist = "dd")

# 4. 建立最终 Logistic 回归模型（rms::lrm）
final_model <- lrm(final_formula, data = mydata, x = TRUE, y = TRUE)
print(final_model)

nom <- nomogram(
  final_model,
  fun = list(function(x) 1 / (1 + exp(-x))),   # logistic 转换
  fun.at = list(c(0.1, 0.5, 0.9, 0.999)),      # ⭐ 自定义刻度
  funlabel = "Predicted Value",
  lp = TRUE
)

plot(nom, xfrac = 0.3)

