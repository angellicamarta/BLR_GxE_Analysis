# Load data
load("prepared_believe_data.RData")
rm(X)

library(gap)

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

# Define batch sizes
batch_size <- 10000
last_batch_size <- 13717

# Batch indices
indices_batch_1 <- 1:batch_size
indices_batch_2 <- (batch_size + 1):(2 * batch_size)
indices_batch_3 <- (2 * batch_size + 1):(3 * batch_size)
indices_batch_4 <- (3 * batch_size + 1):(4 * batch_size)
indices_batch_5 <- (4 * batch_size + 1):(5 * batch_size)
indices_batch_6 <- (5 * batch_size + 1):(5 * batch_size + last_batch_size)

# Function to create and save batch data
create_save_batch <- function(batch_indices, batch_number) {
  GRM_batch <- GRM[batch_indices, batch_indices]
  ids_batch <- phenotype[batch_indices, c("FID", "IID")]
  GRM_batch_vec <- as.vector(GRM_batch)
  
  write.table(ids_batch, file = paste0("GRM_batch_", batch_number, ".grm.id"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  writeBin(GRM_batch_vec, paste0("GRM_batch_", batch_number, ".grm.bin"), size = 4)
  writeBin(rep(1, length(GRM_batch_vec)), paste0("GRM_batch_", batch_number, ".grm.N.bin"), size = 4)
  
  phenotype_batch <- phenotype[batch_indices, ]
  albumin_discrete_batch <- albumin_discrete[batch_indices, ]
  vigorous_activity_batch <- vigorous_activity[batch_indices, ]
  categorical_covariates_batch <- categorical_covariates[batch_indices, ]
  quantitative_covariates_batch <- quantitative_covariates[batch_indices, ]
  study_discrete_batch <- study_discrete[batch_indices, ]
  
  write.table(phenotype_batch, file = paste0("phenotype_batch_", batch_number, ".pheno"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(albumin_discrete_batch, file = paste0("albumin_batch_", batch_number, ".env"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(vigorous_activity_batch, file = paste0("vigorous_batch_", batch_number, ".env"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(categorical_covariates_batch, file = paste0("categorical_covariates_batch_", batch_number, ".covar"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(quantitative_covariates_batch, file = paste0("quantitative_covariates_batch_", batch_number, ".qcovar"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(study_discrete_batch, file = paste0("study_batch_", batch_number, ".env"), row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Create and save batches
create_save_batch(indices_batch_1, 1)
create_save_batch(indices_batch_2, 2)
create_save_batch(indices_batch_3, 3)
create_save_batch(indices_batch_4, 4)
create_save_batch(indices_batch_5, 5)
create_save_batch(indices_batch_6, 6)
