# ============================================================
# File: 01_impute_missing_data.R
# Original file: 2026.1.29随机森林插补死亡数据库（已去除高缺失）.R
# Purpose: Missing-data imputation with missForest for imported malaria and mortality datasets
# Study: Imported malaria analysis
#
# Inputs:
# - introduced260318.csv
# - introduced260316.csv
#
# Outputs:
# - death2rf_raw.csv
# - death2rf_imputed.csv
# - introduced260318rf_raw.csv
# - introduced260318_imputed.csv
#
# Notes:
# - Personal working-directory paths have been removed.
# - Run from the repository root using source("R/01_impute_missing_data.R").
# - Original unmodified script is retained in original_scripts/.
# ============================================================

source(file.path("R", "00_setup.R"))

# setwd() removed; project paths are defined in R/00_setup.R
mydata <- read.csv(input_file("introduced260318.csv"))
# View(mydata)  # Interactive inspection removed for reproducibility
# 每个变量的缺失值数量
na_count <- colSums(is.na(mydata))
na_count

mydata$Outcome <- as.factor(mydata$Outcome)
mydata$Country_of_destination <- as.factor(mydata$Country_of_destination)
mydata$Region_of_acquisition <- as.factor(mydata$Region_of_acquisition)
mydata$Groups_of_economies <- as.factor(mydata$Groups_of_economies)
mydata$Plasmodium_species <- as.factor(mydata$Plasmodium_species)
mydata$Occupation = as.factor(mydata$Occupation)
mydata$Temp <- as.numeric(mydata$Temp)
mydata$Precip <- as.numeric(mydata$Precip)
mydata$Mosquito_borne_species <- as.numeric(mydata$Mosquito_borne_species)

colSums(is.na(mydata))

library(missForest)

set.seed(123)

imp <- missForest(
  mydata,
  maxiter = 10,
  ntree = 500
)

mydata_imp <- mydata

# 只替换需要插补的 6 个变量
vars_to_impute <- c(
  "Region_of_acquisition","Country_of_destination","Groups_of_economies",
  "Plasmodium_species",
  "Temp",
  "Precip","Occupation"
)

mydata_imp[vars_to_impute] <- imp$ximp[vars_to_impute]

colSums(is.na(mydata_imp))

prop.table(table(mydata$Sex))
prop.table(table(mydata_imp$Sex))
prop.table(table(mydata$Age))
prop.table(table(mydata_imp$Age))
prop.table(table(mydata$Region_of_acquisition))
prop.table(table(mydata_imp$Region_of_acquisition))
prop.table(table(mydata$Plasmodium_species))
prop.table(table(mydata_imp$Plasmodium_species))

imp$OOBerror

colSums(is.na(mydata_imp))

lapply(mydata_imp, levels)

table(mydata$Outcome, mydata_imp$Outcome)
write.csv(
  mydata,
  file = processed_file("death2rf_raw.csv"),
  row.names = FALSE
)
write.csv(
  mydata_imp,
  file = processed_file("death2rf_imputed.csv"),
  row.names = FALSE
)
# setwd() removed; project paths are defined in R/00_setup.R
mydata <- read.csv(input_file("introduced260316.csv"))
# View(mydata)  # Interactive inspection removed for reproducibility
# 每个变量的缺失值数量
na_count <- colSums(is.na(mydata))
na_count

mydata$Outcome <- as.factor(mydata$Outcome)
mydata$Region_of_acquisition <- as.factor(mydata$Region_of_acquisition)
mydata$Country_of_destination <- as.factor(mydata$Country_of_destination)
mydata$Groups_of_economies <- as.factor(mydata$Groups_of_economies)
mydata$Plasmodium_species <- as.factor(mydata$Plasmodium_species)
mydata$Mosquito_borne_species <- as.numeric(mydata$Mosquito_borne_species)
mydata$Temp <- as.numeric(mydata$Temp)
mydata$Precip <- as.numeric(mydata$Precip)
mydata$Occupation <- as.factor(mydata$Occupation)
mydata$Secondary <- as.numeric(mydata$Secondary)

colSums(is.na(mydata))

library(missForest)

set.seed(123)

imp <- missForest(
  mydata,
  maxiter = 10,
  ntree = 100
)

mydata_imp <- mydata

# 只替换需要插补的 3 个变量
vars_to_impute <- c(
  "Region_of_acquisition","Country_of_destination","Groups_of_economies",
  "Plasmodium_species","Temp","Precip",
  "Occupation"
)

mydata_imp[vars_to_impute] <- imp$ximp[vars_to_impute]

colSums(is.na(mydata_imp))


prop.table(table(mydata$Region_of_acquisition))
prop.table(table(mydata_imp$Region_of_acquisition))
prop.table(table(mydata$Country_of_destination))
prop.table(table(mydata_imp$Country_of_destination))
prop.table(table(mydata$Groups_of_economies))
prop.table(table(mydata_imp$Groups_of_economies))
prop.table(table(mydata$Plasmodium_species))
prop.table(table(mydata_imp$Plasmodium_species))
prop.table(table(mydata$Temp))
prop.table(table(mydata_imp$Temp))
prop.table(table(mydata$Precip))
prop.table(table(mydata_imp$Precip))
prop.table(table(mydata$Occupation))
prop.table(table(mydata_imp$Occupation))


imp$OOBerror

colSums(is.na(mydata_imp))

lapply(mydata_imp, levels)

table(mydata$Outcome, mydata_imp$Outcome)
write.csv(
  mydata,
  file = processed_file("introduced260318rf_raw.csv"),
  row.names = FALSE
)
write.csv(
  mydata_imp,
  file = processed_file("introduced260318_imputed.csv"),
  row.names = FALSE
)
