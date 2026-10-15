# ============================================================
# File: 05_onward_transmission_firth_logistic.R
# Original file: 2026.5.9 再传播logistic的Firth's Logistic.R
# Purpose: Firth logistic regression for onward transmission risk
# Study: Imported malaria analysis
#
# Inputs:
# - introduced260318_imputed.csv
#
# Outputs:
# - 0407再传播基线表table1_output.csv
# - Firth回归_完整结果表.csv
#
# Notes:
# - Personal working-directory paths have been removed.
# - Run from the repository root using source("R/05_onward_transmission_firth_logistic.R").
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
# Package installation removed for reproducibility. Install "logistf" manually if needed.
library(logistf)  # Firth 回归核心包
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
write.csv(table1_matrix, file = output_table_file("0407再传播基线表table1_output.csv"))

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

## -------- 1. 明确要设定参考值的变量 --------
vars_ref_1 <- c("Country_of_destination")
vars_ref_4 <- c("Plasmodium_species")
vars_ref_2 <- c("Occupation")
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

# 设置参考 = "4"
for (var in vars_ref_4) {
  if (is.factor(mydata[[var]])) {
    mydata[[var]] <- relevel(mydata[[var]], ref = "4")
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
# 构建模型
firth_formula <- as.formula(paste("Outcome ~", paste(selected_significant_names, collapse = " + ")))
firth_model <- logistf(firth_formula, data = mydata)

# 查看结果
summary(firth_model)
# ==============================================
# 提取 Firth 回归：β、Wald、OR(95%CI)、P值
# ==============================================

# 提取系数、标准误、Wald、P值
beta <- firth_model$coefficients
se <- sqrt(diag(firth_model$var))
wald <- (beta/se)^2
p_val <- firth_model$prob

# 提取 OR 和 95%CI
or <- exp(beta)
ci_lower <- exp(beta - 1.96 * se)
ci_upper <- exp(beta + 1.96 * se)
or_ci <- paste0(round(or,3), " (", round(ci_lower,3), "-", round(ci_upper,3), ")")

# 合并成最终表格（你要的4列）
result_table <- data.frame(
  β = round(beta, 4),
  Wald = round(wald, 4),
  OR..95..CI. = or_ci,
  P值 = round(p_val, 4)
)

# 查看表格
cat("\n======= Firth 回归完整结果表 =======\n")
print(result_table)

# 导出 CSV（论文直接用）
write.csv(result_table, output_table_file("Firth回归_完整结果表.csv"), fileEncoding = "UTF-8")

