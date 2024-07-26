# Load necessary libraries
library(gap)

# Load data
load("prepared_believe_data.RData")
rm(X)

# Read data files
phenotype <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/BELIEVE.pheno", header = FALSE, col.names = c("FID", "IID", "T2D_status"))
albumin_discrete <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/Albumin_discrete_file.env", header = FALSE, col.names = c("FID", "IID", "Albumin"))
vigorous_activity <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/Vigorous10MinsActivity.env", header = FALSE, col.names = c("FID", "IID", "Vigorous10MinsActivity"))
categorical_covariates <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/C.covar", header = FALSE, col.names = c("FID", "IID", "Sex"))
quantitative_covariates <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/Q.covar", header = FALSE, col.names = c("FID", "IID", "Age", "BMI"))
study_discrete <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/study_discrete.env", header = FALSE, col.names = c("FID", "IID", "Study"))

# Combine all data for filtering
combined_data <- cbind(phenotype$T2D_status, albumin_discrete$Albumin, vigorous_activity$Vigorous10MinsActivity, categorical_covariates$Sex, quantitative_covariates$Age, quantitative_covariates$BMI, study_discrete$Study)

# Identify valid rows (excluding those with -9)
valid_indices <- apply(combined_data, 1, function(row) all(row != -9))

# Filter each dataset individually
phenotype <- phenotype[valid_indices, ]
albumin_discrete <- albumin_discrete[valid_indices, ]
vigorous_activity <- vigorous_activity[valid_indices, ]
categorical_covariates <- categorical_covariates[valid_indices, ]
quantitative_covariates <- quantitative_covariates[valid_indices, ]
study_discrete <- study_discrete[valid_indices, ]

# Use the valid indices to filter the GRM
GRM <- GRM[valid_indices, valid_indices]

# Number of batches
n_batches <- 6

# Create indices for the 6 batches
batch_indices <- split(1:nrow(GRM), sort(1:nrow(GRM) %% n_batches))

# Extract indices for Batch 1
indices <- batch_indices[[1]]

# Filter GRM and data for Batch 1
GRM_batch_1 <- GRM[indices, indices]
phenotype_batch_1 <- phenotype[indices, ]
albumin_discrete_batch_1 <- albumin_discrete[indices, ]
vigorous_activity_batch_1 <- vigorous_activity[indices, ]
categorical_covariates_batch_1 <- categorical_covariates[indices, ]
quantitative_covariates_batch_1 <- quantitative_covariates[indices, ]
study_discrete_batch_1 <- study_discrete[indices, ]

# Create IDs for Batch 1
ids_batch_1 <- phenotype_batch_1[, c("FID", "IID")]

# Write GRM files for Batch 1
write.table(ids_batch_1, file = "GRM_batch_1.grm.id", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Manually write GRM to binary file
write_grm_bin <- function(grm, prefix) {
  binfile <- paste0(prefix, ".grm.bin")
  n <- nrow(grm)
  grm_vec <- as.vector(grm)
  writeBin(grm_vec, binfile, size = 4)
  nfile <- paste0(prefix, ".grm.N.bin")
  writeBin(rep(1L, length(grm_vec)), nfile, size = 4)
}

write_grm_bin(GRM_batch_1, "GRM_batch_1")

# Write data files for Batch 1
write.table(phenotype_batch_1, file = "phenotype_batch_1.pheno", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(albumin_discrete_batch_1, file = "albumin_batch_1.env", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(vigorous_activity_batch_1, file = "vigorous_batch_1.env", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(categorical_covariates_batch_1, file = "categorical_covariates_batch_1.covar", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(quantitative_covariates_batch_1, file = "quantitative_covariates_batch_1.qcovar", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(study_discrete_batch_1, file = "study_batch_1.env", row.names = FALSE, col.names = FALSE, quote = FALSE)

