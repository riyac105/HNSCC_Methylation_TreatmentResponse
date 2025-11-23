##############################################
# 03_predictive_models_glmnet_RF_ROC_pAUC.R
##############################################

library(tidyverse)
library(glmnet)
library(randomForest)
library(pROC)

set.seed(123)

# Load merged data (includes risk_score)
dat <- readRDS("outputs/dat_TCGA_MRS_clinical.rds")

# Binary response: 1 = Nonresponder, 0 = Responder
dat <- dat %>%
  mutate(
    response_bin = ifelse(Treatment.Response %in% c("Progressive Disease",
                                                    "Persistent Disease"),
                          1, 0)
  )

#-------------------------
# Split train / test
#-------------------------

train_idx <- sample(seq_len(nrow(dat)), size = 0.7 * nrow(dat))
train <- dat[train_idx, ]
test  <- dat[-train_idx, ]

#-------------------------
# Model 1: glmnet with clinical +/- MRS
#-------------------------

# Design matrices
x_train_clin <- model.matrix(
  response_bin ~ Age + stage_simple + HPV + Alcohol.History + Smoking.History,
  data = train
)[, -1]

x_train_clin_mrs <- model.matrix(
  response_bin ~ Age + stage_simple + HPV + Alcohol.History +
    Smoking.History + risk_score,
  data = train
)[, -1]

y_train <- train$response_bin

x_test_clin <- model.matrix(
  response_bin ~ Age + stage_simple + HPV + Alcohol.History + Smoking.History,
  data = test
)[, -1]

x_test_clin_mrs <- model.matrix(
  response_bin ~ Age + stage_simple + HPV + Alcohol.History +
    Smoking.History + risk_score,
  data = test
)[, -1]

y_test <- test$response_bin

# Fit Elastic Net (alpha can be tuned; use alpha=0.5 as example)
fit_clin <- cv.glmnet(x_train_clin, y_train, family = "binomial", alpha = 0.5)
fit_clin_mrs <- cv.glmnet(x_train_clin_mrs, y_train, family = "binomial", alpha = 0.5)

# Predict probabilities on test set
pred_clin     <- predict(fit_clin,     newx = x_test_clin,     s = "lambda.min", type = "response")
pred_clin_mrs <- predict(fit_clin_mrs, newx = x_test_clin_mrs, s = "lambda.min", type = "response")

# ROC / AUC and DeLong test
roc_clin     <- roc(y_test, as.numeric(pred_clin),     quiet = TRUE)
roc_clin_mrs <- roc(y_test, as.numeric(pred_clin_mrs), quiet = TRUE)

auc(roc_clin)
auc(roc_clin_mrs)

# DeLong test
roc.test(roc_clin, roc_clin_mrs, method = "delong")

# Export glmnet coefficients for supplement
coef_clin_mrs <- coef(fit_clin_mrs, s = "lambda.min")
coef_df <- data.frame(
  Variable    = rownames(coef_clin_mrs),
  Coefficient = as.numeric(coef_clin_mrs)
)
write.csv(coef_df, "outputs/logistic_glmnet_coefficients_MRS.csv", row.names = FALSE)

#-------------------------
# Random Forest with M-values
#-------------------------

# Load M-value matrix (rows = CpGs, cols = samples); top 60 CpGs
mval_mat <- readRDS("data/Mvalue_top60_TCGA.rds")

# Ensure order of columns matches dat$patient_ID
mval_mat <- mval_mat[, dat$patient_ID, drop = FALSE]

# Build RF train / test data
mval_train <- t(mval_mat[, train_idx, drop = FALSE])
mval_test  <- t(mval_mat[, test_idx <- setdiff(seq_len(nrow(dat)), train_idx), drop = FALSE])

rf_train <- train %>%
  select(response_bin, Age, Stage, HPV, Alcohol.History,
         Smoking.History, sub_site_grouped) %>%
  bind_cols(as.data.frame(mval_train))

rf_test <- test %>%
  select(response_bin, Age, Stage, HPV, Alcohol.History,
         Smoking.History, sub_site_grouped) %>%
  bind_cols(as.data.frame(mval_test))

rf_train$response_bin <- factor(rf_train$response_bin, levels = c(0,1))

# Example RF (you can tune mtry, ntree, sampsize as in your analysis)
fit_rf <- randomForest(
  response_bin ~ .,
  data = rf_train,
  ntree = 1000,
  importance = TRUE
)

# Predict on test
pred_rf <- predict(fit_rf, newdata = rf_test, type = "prob")[, "1"]

roc_rf <- roc(y_test, as.numeric(pred_rf), quiet = TRUE)
auc(roc_rf)

# Concordant partial AUC: FPR ≤ 0.30, TPR ≥ 0.70
# This corresponds roughly to specificity ≥ 0.70
pAUC <- auc(
  roc_rf,
  partial.auc = c(0.7, 1),  # specificity range
  partial.auc.focus = "specificity",
  partial.auc.correct = TRUE
)

# Normalize to [0,1]
pAUC_norm <- as.numeric(pAUC) / (1 - 0.7)  # divide by maximum area in that region
pAUC_norm

# Save AUC/pAUC results
perf_df <- tibble(
  model = c("Clin_only", "Clin+MRS", "RF_Mvalues"),
  AUC   = c(auc(roc_clin), auc(roc_clin_mrs), auc(roc_rf))
)

write.csv(perf_df, "outputs/predictive_model_AUCs.csv", row.names = FALSE)
