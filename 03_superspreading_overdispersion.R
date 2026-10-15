# ============================================================
# File: 03_superspreading_overdispersion.R
# Original file: 2025.11.4 超级传播的判定.R
# Purpose: Superspreading and overdispersion analysis using secondary-case counts
# Study: Imported malaria analysis
#
# Inputs:
# - second2.csv
#
# Outputs:
# - Console output only unless user adds write commands
#
# Notes:
# - Personal working-directory paths have been removed.
# - Run from the repository root using source("R/03_superspreading_overdispersion.R").
# - Original unmodified script is retained in original_scripts/.
# ============================================================

source(file.path("R", "00_setup.R"))

# setwd() removed; project paths are defined in R/00_setup.R
# 安装与加载必要包
library(fitdistrplus)
library(MASS)

mydata <- read.csv(input_file("second2.csv"))

# 检查数据
head(mydata)
mydata$secondary_cases <- as.numeric(mydata$secondary_cases)
# 计算平均传播数 R
R_mean <- mean(mydata$secondary_cases)
R_mean

# 平均传播数
R_mean <- mean(mydata$secondary_cases)

# Poisson 95% CI
R_ci <- poisson.test(
  sum(mydata$secondary_cases),
  T = length(mydata$secondary_cases)
)$conf.int

cat("R =", round(R_mean, 4), "\n")
cat("95% CI =", 
    paste0("(", round(R_ci[1], 4), ", ", round(R_ci[2], 4), ")"), "\n")


# -----------------------------
# 估计离散参数 k （负二项分布拟合）
# -----------------------------
# 使用MASS包中的fitdistr函数进行最大似然估计
fit_nb <- fitdistr(mydata$secondary_cases, "Negative Binomial")

fit_nb

# 提取参数
size_k <- fit_nb$estimate["size"]   # 即离散参数k
mu_R   <- fit_nb$estimate["mu"]     # 即平均传播数R的拟合值

cat("估计结果：\n")
cat("平均传播数 R =", round(mu_R, 3), "\n")
cat("离散参数 k =", round(size_k, 3), "\n")

# 拟合负二项
fit_nb <- fitdistr(mydata$secondary_cases, "Negative Binomial")

# 提取 k 和 SE
k_hat <- fit_nb$estimate["size"]
k_se  <- fit_nb$sd["size"]

# Wald 95% CI
k_lower <- k_hat - 1.96 * k_se
k_upper <- k_hat + 1.96 * k_se

cat("k =", round(k_hat, 4), "\n")
cat("95% CI for k (Wald): (",
    round(k_lower, 4), ",",
    round(k_upper, 4), ")\n")


# -----------------------------
# 与泊松分布进行似然比检验（检验是否存在聚集性）
# -----------------------------
fit_pois <- fitdistr(data$secondary_cases, "Poisson")

# 计算似然比统计量
LR_stat <- 2 * (fit_nb$loglik - fit_pois$loglik)
p_value <- pchisq(LR_stat, df = 1, lower.tail = FALSE)

cat("\n似然比检验结果：\n")
cat("LR统计量 =", round(LR_stat, 3), "\n")
cat("p值 =", round(p_value, 5), "\n")

if(p_value < 0.05){
  cat("结果：拒绝泊松假设，传播呈现超分散（聚集）模式。\n")
} else {
  cat("结果：未发现显著偏离泊松分布（随机传播）。\n")
}
