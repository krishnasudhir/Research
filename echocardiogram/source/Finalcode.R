rm(list=ls())
setwd("~/NCSSM/Intro to Computational Sciences/Research Project")
library(Metrics)
library(InformationValue)
library(caTools)
library(ggplot2)
library(ROCR)
library(boot)

#---------------------------------------------------------------------------------------#
# Load the data from CSV file
#---------------------------------------------------------------------------------------#

echocardiogram <- read.csv("echocardiogram.csv")

importantecho <- echocardiogram[, c("age", "aliveat1", "pericardialeffusion", "wallmotion.index", "wallmotion.score", "fractionalshortening", "lvdd", "epss")]

importantecho$age[is.na(importantecho$age)] <- median(importantecho$age,na.rm=T)

#---------------------------------------------------------------------------------------#
# There are missing values that need to be cleaned up
#---------------------------------------------------------------------------------------#

cleanEcho <- na.omit(importantecho)
nrow(cleanEcho)

#---------------------------------------------------------------------------------------#
# Set seed for the random generators
#---------------------------------------------------------------------------------------#

set.seed(134)


#---------------------------------------------------------------------------------------#
# Get 75% of the dataset as training data
#---------------------------------------------------------------------------------------#

sample_size <- floor(0.65*(nrow(cleanEcho)))
echo_flag <- sample(seq_len(nrow(cleanEcho)), size = sample_size)

#------------------------------ ---------------------------------------------------------#
# Split into training data and test data (75/25 split)
#---------------------------------------------------------------------------------------#

echotrain <- cleanEcho[echo_flag, ]
echotest <- cleanEcho[-echo_flag, ]

nrow(echotrain)
nrow(echotest)

#---------------------------------------------------------------------------------------#
# Building the model
#---------------------------------------------------------------------------------------#
#
# Univariate models
#
formula1 <- aliveat1 ~ age
formula2 <- aliveat1 ~ pericardialeffusion
formula3 <- aliveat1 ~ fractionalshortening
formula4 <- aliveat1 ~ wallmotion.index
formula5 <- aliveat1 ~ wallmotion.score
formula6 <- aliveat1 ~ lvdd
formula7 <- aliveat1 ~ epss
formula8 <- aliveat1 ~ age + fractionalshortening + wallmotion.index
formula9 <- aliveat1 ~ age + fractionalshortening + wallmotion.score
formula10 <- aliveat1 ~ age + fractionalshortening + wallmotion.index + lvdd
formula11 <- aliveat1 ~ age + fractionalshortening + wallmotion.index + lvdd + epss
formula12 <- aliveat1 ~ age + fractionalshortening + wallmotion.index + epss
formula13 <- aliveat1 ~ fractionalshortening + wallmotion.index + epss + lvdd
formula14 <- aliveat1 ~ fractionalshortening + wallmotion.index + epss + lvdd + pericardialeffusion

formula_list_1 <- c(formula1, formula2, formula3, formula4, formula5, formula6, 
                    formula7, formula8, formula9, formula10, formula11, formula12, formula13, formula14)

regression_enmasse <- function(formula_list) {
  for (fmla in formula_list) {
    set.seed(1234)
    model1 <- glm(fmla, family=binomial(), data = echotrain)
    
    prediction1 <- predict(model1, newdata = echotest, type = 'response')
    fitted.results <- ifelse(prediction1 > 0.5, 1, 0)
    
    predictionerror <- mean(fitted.results != echotest$aliveat1)
    
    score <- rmse(prediction1, echotest$aliveat1)
    print("")
    print("-------------------------------------------------------------------------------------")
    print(fmla)
    print("-------------------------------------------------------------------------------------")
    print(summary(model1))
    print(paste("ACCURACY: ", 1-predictionerror))
    print(paste("RMSE: ", score))
    print("Concordance Analysis:")
    print(Concordance(echotest$aliveat1, fitted.results))
    print("Misclassification Error Analysis:")
    print(misClassError(echotest$aliveat1, fitted.results))
    print("Specificity Analysis:")
    print(specificity(echotest$aliveat1, fitted.results))
    print("Sensitivity Analysis:")
    print(sensitivity(echotest$aliveat1, fitted.results))
    print("Confusion Matrix:")
    print(confusionMatrix(echotest$aliveat1, fitted.results))
    sensMatrix <- plotROC(echotest$aliveat1, prediction1, Show.labels = F, returnSensitivityMat = T)
    #print(sensMatrix)
    ks_plot(echotest$aliveat1, prediction1)
  }
}

regression_enmasse(formula_list_1)


# This is to figure out multicollinearity using the VIF score. All VIF scores are greater than or equal to 1 however any values greater than 10 is concerning. Anything above 2.5 can be slightly concerning but all values below 2 so no problemo :) 

vifformula1 <- aliveat1 ~ age + pericardialeffusion + fractionalshortening + wallmotion.index + lvdd + epss 

model1 <- glm(vifformula1, family=binomial(), data = echotrain)
summary(model1)
vif(model1)

#---------------------------------------------------------------------------------------#
# Perform a cross validation 
#---------------------------------------------------------------------------------------#
 
final_model <- glm(formula8, family = binomial(), data = cleanEcho)
summary(final_model)

#---------------------------------------------------------------------------------------#
# Leave one out Cross Validation
#---------------------------------------------------------------------------------------#

cost <- function(r, pi = 0) mean(abs(r-pi) > 0.5)
cv.loocv <- cv.glm(cleanEcho, final_model, cost)
cv.loocv$delta

#---------------------------------------------------------------------------------------#
# K-fold Cross Validation
#---------------------------------------------------------------------------------------#

data2 <- data.frame()

for(idx in 2:12){
  set.seed(134)
  model2 <- glm(formula8, family = binomial(), data = cleanEcho)
  cv.out <- cv.glm(cleanEcho, model2, cost, K=idx)
  
  print(paste(idx, cv.out$K, cv.out$delta))
  
  data2 <- rbind(data2, c(cv.out$K, cv.out$delta[1], cv.out$delta[2]))
}
data2

#This is to show hyper parameter tuning
plot(data2$K, data2$`Adj Error`, type="o")

