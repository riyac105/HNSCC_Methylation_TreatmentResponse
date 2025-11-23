##############################################
# 03_predictive_models_glmnet_RF_ROC_pAUC.R
##############################################

library(tidyverse)
library(glmnet)
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
