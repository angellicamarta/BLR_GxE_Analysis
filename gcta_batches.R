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

# Create function to process each batch
process_batch <- function(indices, batch_number) {
  GRM_batch <- GRM[indices, indices]
  ids_batch <- phenotype[indices, c("FID", "IID")]
  
  # Convert GRM matrix into vector form
  GRM_batch_vec <- as.vector(GRM_batch)
  
  # Write GRM files for the batch using WriteGRMBin
  WriteGRMBin(prefix = paste0("GRM_batch_", batch_number), grm = GRM_batch_vec, N = nrow(GRM_batch), id = ids_batch, size = 4)
}

# Process each batch
process_batch(indices_batch_1, 1)
process_batch(indices_batch_2, 2)
process_batch(indices_batch_3, 3)
process_batch(indices_batch_4, 4)
process_batch(indices_batch_5, 5)
process_batch(indices_batch_6, 6)

# Ensure environmental variables and covariates are split correctly for each batch
write.table(albumin_discrete[indices_batch_1, ], file = "albumin_batch_1.env", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(albumin_discrete[indices_batch_2, ], file = "albumin_batch_2.env", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(albumin_discrete[indices_batch_3, ], file = "albumin_batch_3.env", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(albumin_discrete[indices_batch_4, ], file = "albumin_batch_4.env", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(albumin_discrete[indices_batch_5, ], file = "albumin_batch_5.env", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(albumin_discrete[indices_batch_6, ], file = "albumin_batch_6.env", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(vigorous_activity[indices_batch_1, ], file = "vigorous_batch_1.env", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(vigorous_activity[indices_batch_2, ], file = "vigorous_batch_2.env", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(vigorous_activity[indices_batch_3, ], file = "vigorous_batch_3.env", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(vigorous_activity[indices_batch_4, ], file = "vigorous_batch_4.env", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(vigorous_activity[indices_batch_5, ], file = "vigorous_batch_5.env", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(vigorous_activity[indices_batch_6, ], file = "vigorous_batch_6.env", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(categorical_covariates[indices_batch_1, ], file = "categorical_covariates_batch_1.covar", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(categorical_covariates[indices_batch_2, ], file = "categorical_covariates_batch_2.covar", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(categorical_covariates[indices_batch_3, ], file = "categorical_covariates_batch_3.covar", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(categorical_covariates[indices_batch_4, ], file = "categorical_covariates_batch_4.covar", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(categorical_covariates[indices_batch_5, ], file = "categorical_covariates_batch_5.covar", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(categorical_covariates[indices_batch_6, ], file = "categorical_covariates_batch_6.covar", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(quantitative_covariates[indices_batch_1, ], file = "quantitative_covariates_batch_1.qcovar", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(quantitative_covariates[indices_batch_2, ], file = "quantitative_covariates_batch_2.qcovar", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(quantitative_covariates[indices_batch_3, ], file = "quantitative_covariates_batch_3.qcovar", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(quantitative_covariates[indices_batch_4, ], file = "quantitative_covariates_batch_4.qcovar", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(quantitative_covariates[indices_batch_5, ], file = "quantitative_covariates_batch_5.qcovar", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(quantitative_covariates[indices_batch_6, ], file = "quantitative_covariates_batch_6.qcovar", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(study_discrete[indices_batch_1, ], file = "study_batch_1.env", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(study_discrete[indices_batch_2, ], file = "study_batch_2.env", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(study_discrete[indices_batch_3, ], file = "study_batch_3.env", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(study_discrete[indices_batch_4, ], file = "study_batch_4.env", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(study_discrete[indices_batch_5, ], file = "study_batch_5.env", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(study_discrete[indices_batch_6, ], file = "study_batch_6.env", quote = FALSE, row.names = FALSE, col.names = FALSE)

