# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Sat May 25 12:27:39 EDT 2019

################################################################
# LIBRARIES #
library(Biobase)
library(GEOquery)
library(limma)
library(Metrics)
library(InformationValue)
library(car)

#------CLEAN THINGS UP -------#
rm(list=ls())

#--------------------------------------------------------------#
# load series and platform data from GEO
#--------------------------------------------------------------#
#gset <- gset_backup

gset <- getGEO("GSE73002", GSEMatrix =TRUE, AnnotGPL=FALSE)

if (length(gset) > 1) idx <- grep("GPL18941", attr(gset, "names")) 
else idx <- 1
gset <- gset[[idx]]

#----------------------------------------------------------------------#
# make proper column names to match toptable 
#----------------------------------------------------------------------#
fvarLabels(gset) <- make.names(fvarLabels(gset))

#----------------------------------------------------------------------#
#Seperating all diagnosises as breast cancer and control
#----------------------------------------------------------------------#
YY <- ifelse(gset$diagnosis == "breast cancer", "1", 
             ifelse(gset$diagnosis == "non-cancer", "0", "X"))

sel <- which(YY != "X")
gset <- gset[ ,sel]

#----------------------------------------------------------------------#
# group names for all samples
#----------------------------------------------------------------------#
masks <- ifelse(gset$diagnosis == "breast cancer", "G1", "G0")


#----------------------------------------------------------------------#
# log2 transform
#----------------------------------------------------------------------#
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

#----------------------------------------------------------------------#
# set up the data and proceed with analysis
#----------------------------------------------------------------------#

fl <- as.factor(masks)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

#----------------------------------------------------------------------#
# FINAL Top Ten Table
#----------------------------------------------------------------------#
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","miRNA_ID","miRNA_ID_LIST"))
write.table(tT, file=stdout(), row.names=F, sep="\t")


#--------------------------------------------------------------------------------------------------

dataframe1 <- data.frame(gset) 
Y <- ifelse(gset$`diagnosis:ch1` == "breast cancer", 1, 0)
dataframe1$Y <- Y
miRNAIDs <- tT$ID[1:100]
final_dataframe <- subset(dataframe1, select=c(miRNAIDs, "Y"))
final_dataframe <- na.omit(final_dataframe)

#----------------------------------------------------------------------#
#-------------------------Create model formula-------------------------#
#----------------------------------------------------------------------#
idx <- 20
fmla <- paste("Y~")
for (var1 in miRNAIDs[1:idx]) {
  
  fmla <- paste(fmla, var1, sep = "") 
  if (var1 != miRNAIDs[idx])
    fmla <- paste(fmla, "+", sep = "")
}

#----------------------------------------------------------------------#
# Create the final data frame with all data 
#----------------------------------------------------------------------#

final_dataframe <- final_dataframe[sample(nrow(final_dataframe)),]

#----------------------------------------------------------------------#
# Set seed to ensure the model may be replicated. 
#----------------------------------------------------------------------#
set.seed(143)

#----------------------------------------------------------------------#
# Split the data frame into training and testing data sets. 
# 65% of the existing data frame will be used for training the model
# and the remaning 35% of the data will be used to test the model 
# for accuracy. 
#----------------------------------------------------------------------#
sample_size <- floor(0.65*(nrow(final_dataframe)))
my_flag <- sample(seq_len(nrow(final_dataframe)), size = sample_size)
df_train <- final_dataframe[my_flag, ]
df_test <- final_dataframe[-my_flag, ]

nrow(df_train)
nrow(df_test)

#----------------------------------------------------------------------#
#-------------------------The Model------------------------------------#
#----------------------------------------------------------------------#
model1 <- glm(fmla, family=binomial(), data = df_train, maxit = 100)
summary(model1)

# Uncomment to compare the models
# anova(model1, model2, test = "LRT")

#----------------------------------------------------------------------#
#-----------------------Prediction and Error Analysis------------------#
#----------------------------------------------------------------------#
prediction1 <- predict(model1, newdata = df_test, type = "response")
fitted_results <- ifelse(prediction1 > 0.5, 1, 0)
prediction_error <- mean(fitted_results != df_test$Y)

#----------------------------------------------------------------------#
#-------------------------Other Analyses-------------------------------#
#----------------------------------------------------------------------#
rmsescore <- rmse(prediction1, df_test$Y)
print(" ")
print("----------------------------------------------------------------------------")
print(fmla)
print("----------------------------------------------------------------------------")
print(summary(model1))
print(paste("ACCURACY:", 1-prediction_error))
print(paste("RMSE:", rmsescore) )
print("Concordance Analysis:")
print(Concordance(df_test$Y, fitted_results))
print("Misclassification Error Analysis:")
print(misClassError(df_test$Y, fitted_results))
print("Specificity Analysis:")
print(specificity(df_test$Y, fitted_results))
print("Confusion Matrix:")
print(confusionMatrix(df_test$Y, fitted_results))
sensMatrix <- plotROC(df_test$Y, prediction1, Show.labels = F, returnSensitivityMat = T)
ks_plot(df_test$Y, prediction1)

#---------------------------------------------------------------------------------------------------#
# This is a VIF model. A VIF score is generated to figure out multicollinearity using the VIF score. 
# All VIF scores are greater than or equal to 1 however any values greater than 10 is concerning. 
# Anything above 2.5 can be slightly concerning but all values below 2 are not concerning.
#---------------------------------------------------------------------------------------------------#
vif(model1)
