# ============================================================
# File: 06_mortality_logistic_regression.R
# Original file: 2026.5.9- 死亡逻辑回归代码.R
# Purpose: Conventional logistic regression for imported-malaria mortality risk
# Study: Imported malaria analysis
#
# Inputs:
# - death260318.csv
#
# Outputs:
# - 0325死亡基线表table1_output.csv
# - 59测_敏感死亡播Multiv_logistic_rms格式调整.csv
#
# Notes:
# - Personal working-directory paths have been removed.
# - Run from the repository root using source("R/06_mortality_logistic_regression.R").
# - Original unmodified script is retained in original_scripts/.
# ============================================================

source(file.path("R", "00_setup.R"))

# setwd() removed; project paths are defined in R/00_setup.R
mydata <- read.csv(input_file("death260318.csv"))
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
write.csv(table1_matrix, file = output_table_file("0325死亡基线表table1_output.csv"))

# ----------------------------
# 筛选 p < 0.05 的显著变量
# ----------------------------
selected_significant_names <- c(
  "Region_of_acquisition",
  "Country_of_destination",
  "Plasmodium_species",
  "Occupation",
  "Temp","Precip"
)

## -------- 1. 明确要设定参考值的变量 --------
vars_ref_1 <- c("Country_of_destination")
vars_ref_4 <- c("Plasmodium_species")
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

# 设置参考 = "4"
for (var in vars_ref_4) {
  if (is.factor(mydata[[var]])) {
    mydata[[var]] <- relevel(mydata[[var]], ref = "4")
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


# 7) 导出
write.csv(
  result_table,
  "59测_敏感死亡播Multiv_logistic_rms格式调整.csv",
  row.names = FALSE
)

# 查看
print(result_table)