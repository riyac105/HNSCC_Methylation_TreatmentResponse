library(tidyverse)
library(dplyr)
library(randomForest)
library(class)
library(caret)
library(pROC)

#Part A: All Samples
#STEP 1: Read in data from input file with predictor variables and response
RFInput <- read.csv("RFInput.csv")  #This was from previous steps commented out
RFInput <- RFInput[,-c(1, 2)] #Removing unecessary column 1 (row number) and column 2 (patient ID)
RFInput$Treatment.Response <- as.factor(RFInput$Treatment.Response)
RFInput <- RFInput[, colSums(is.na(RFInput)) == 0] #removing any column with NA values

#STEP 2: 70-30 Train-Test Split
set.seed(0)
trainIndex <- sample(1:nrow(RFInput), 0.7 * nrow(RFInput))
trainData_Rev <- RFInput[trainIndex, ]
testData_Rev <- RFInput[-trainIndex, ]

#STEP 3: Train RF Model using 70% train split cohort
set.seed(0)
rf <- randomForest(Treatment.Response ~ ., data = trainData_Rev, mtry=45, maxnodes = 25, ntrees=500,  sampsize = c("No Response" = 11, "Response" = 11))
print(rf)

#Step 4: Test RF Model using 30% test split cohort
predictions_Rev<-predict(rf, testData_Rev)
confusionMatrix(predictions_Rev, testData_Rev$Treatment.Response)

#Step 5: Display ROC Curves
predictions_ROC <- predict(rf, newdata = testData_Rev, type = "prob")[, "No Response"]
roc_Rev <- roc(testData_Rev$Treatment.Response, predictions_ROC)
#Display ROC Curve
plot(roc_Rev, col = "red", lwd = 2, main = "ROC (TEST): Clinical vs Clinical+Methylation\nAUROC: Clin=0.6588, Clin+Meth=0.8986", xlab = "False Positive Rate (1 - Specificity)",
     ylab = "True Positive Rate (Sensitivity)",
     legacy.axes = TRUE, grid=TRUE)
#Display AUC
pROC::auc(roc_Rev)

#Step 6: Repeating steps 2-5 for clinical covariates alone
RFInput_clin <- RFInput[, 41:47]

set.seed(0)
trainIndex_clin <- sample(1:nrow(RFInput_clin), 0.7 * nrow(RFInput_clin))
trainData_clinRev <- RFInput_clin[trainIndex, ]
testData_clinRev <- RFInput_clin[-trainIndex, ]

set.seed(0)
rf_clin <- randomForest(Treatment.Response ~ ., data = trainData_clinRev, mtry=45, maxnodes = 25, ntrees=500,  sampsize = c("No Response" = 11, "Response" = 11))
print(rf_clin)

predictions_Rev_clin<-predict(rf_clin, testData_clinRev)
confusionMatrix(predictions_Rev_clin, testData_clinRev$Treatment.Response)

predictions_clin_ROC <- predict(rf_clin, newdata = testData_clinRev, type = "prob")[, "No Response"]
roc_clin_Rev <- roc(testData_clinRev$Treatment.Response, predictions_clin_ROC)
pROC::auc(roc_clin_Rev)

plot(roc_clin_Rev,
     col = "steelblue3",
     lwd = 2,
     add = TRUE,
     print.auc = FALSE)

#Step 7: delong's test for statistical significance in difference between ROC curves
test <- roc.test(roc_Rev, roc_clin_Rev, method="delong")
print(test)

#Step 8: concordant partial AUC calculation
#following framework proposed by Carrington et al. (doi: 10.1186/s12911-019-1014-6. PMID: 31906931. PMCID: PMC6945414)
#Thresholds: TRP>70%, FPR<30%
pAUCSens <- pROC::auc(roc_Rev, partial.auc=c(0.7,1), partial.auc.focus="sens")
pAUCSpec <- pROC::auc(roc_Rev, partial.auc=c(0.7,1), partial.auc.focus="spec")
cpAUC <- ((0.5)*pAUCSens)+((0.5)*pAUCSpec)
#Normalization of concordant partial AUC
cpAUCn <- cpAUC/(0.3)
print(cpAUCn)
