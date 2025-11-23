##############################################
# 01_compute_MRS_and_merge_clinical.R
# - Computes methylation risk score (MRS)
# - Merges with clinical data
##############################################

# Load libraries
library(tidyverse)

#-------------------------
# 1. Load input data
#-------------------------

# Beta matrix: rows = CpGs, columns = samples
# Example: beta_mat[cpg, sample]
beta_mat <- readRDS("data/beta_matrix_TCGA_HNSC.rds")

# Limma results with top 300 CpGs used for MRS
# Required columns: CpG_ID, logFC
limma_res <- read.csv("data/limma_top300_for_MRS.csv")

# Clinical data (one row per patient)
# Required columns: patient_ID, Treatment.Response, Age, Stage, HPV,
# Alcohol.History, Smoking.History, OS.time, OS.event, sub_site_grouped
clin <- read.csv("data/clinical_TCGA_HNSC.csv")

# Ensure CpG IDs match rownames in beta_mat
stopifnot(all(limma_res$CpG_ID %in% rownames(beta_mat)))

# Subset beta matrix to MRS CpGs
beta_mrs <- beta_mat[limma_res$CpG_ID, , drop = FALSE]

#-------------------------
# 2. Z-score beta values across patients
#-------------------------
beta_mrs_z <- t(scale(t(beta_mrs)))  # z-score each CpG across samples

# Check dimensions
dim(beta_mrs_z)

#-------------------------
# 3. Compute MRS for each patient
#-------------------------

# Make sure weights vector is aligned to rows of beta_mrs_z
limma_res <- limma_res %>%
  arrange(match(CpG_ID, rownames(beta_mrs_z)))

stopifnot(all(limma_res$CpG_ID == rownames(beta_mrs_z)))

weights <- limma_res$logFC  # vector of length = number of CpGs

# risk_score for each sample = sum_j (z_ij * logFC_j)
mrs_vec <- as.numeric(t(weights) %*% beta_mrs_z)  # 1 x N_samples

# Put into data frame
mrs_df <- tibble(
  sample_ID = colnames(beta_mrs_z),
  risk_score = mrs_vec
)

#-------------------------
# 4. Merge with clinical data
#-------------------------

# Ensure sample_ID matches patient_ID
dat_merged <- clin %>%
  inner_join(mrs_df, by = c("patient_ID" = "sample_ID"))

# Save merged data for downstream analyses
write.csv(dat_merged, "outputs/dat_TCGA_MRS_clinical.csv", row.names = FALSE)
saveRDS(dat_merged, "outputs/dat_TCGA_MRS_clinical.rds")
