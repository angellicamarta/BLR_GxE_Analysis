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

# Write GRM files for Batch 1
GRM_batch_1 <- GRM[indices_batch_1, indices_batch_1]
ids_batch_1 <- phenotype[indices_batch_1, c("FID", "IID")]
GRM_batch_1_vec <- as.vector(GRM_batch_1)
write.table(ids_batch_1, file = "GRM_batch_1.grm.id", quote = FALSE, row.names = FALSE, col.names = FALSE)
writeBin(GRM_batch_1_vec, "GRM_batch_1.grm.bin", size = 4)
writeBin(rep(1, length(GRM_batch_1_vec)), "GRM_batch_1.grm.N.bin", size = 4)

# Write GRM files for Batch 2
GRM_batch_2 <- GRM[indices_batch_2, indices_batch_2]
ids_batch_2 <- phenotype[indices_batch_2, c("FID", "IID")]
GRM_batch_2_vec <- as.vector(GRM_batch_2)
write.table(ids_batch_2, file = "GRM_batch_2.grm.id", quote = FALSE, row.names = FALSE, col.names = FALSE)
writeBin(GRM_batch_2_vec, "GRM_batch_2.grm.bin", size = 4)
writeBin(rep(1, length(GRM_batch_2_vec)), "GRM_batch_2.grm.N.bin", size = 4)

# Write GRM files for Batch 3
GRM_batch_3 <- GRM[indices_batch_3, indices_batch_3]
ids_batch_3 <- phenotype[indices_batch_3, c("FID", "IID")]
GRM_batch_3_vec <- as.vector(GRM_batch_3)
write.table(ids_batch_3, file = "GRM_batch_3.grm.id", quote = FALSE, row.names = FALSE, col.names = FALSE)
writeBin(GRM_batch_3_vec, "GRM_batch_3.grm.bin", size = 4)
writeBin(rep(1, length(GRM_batch_3_vec)), "GRM_batch_3.grm.N.bin", size = 4)

# Write GRM files for Batch 4
GRM_batch_4 <- GRM[indices_batch_4, indices_batch_4]
ids_batch_4 <- phenotype[indices_batch_4, c("FID", "IID")]
GRM_batch_4_vec <- as.vector(GRM_batch_4)
write.table(ids_batch_4, file = "GRM_batch_4.grm.id", quote = FALSE, row.names = FALSE, col.names = FALSE)
writeBin(GRM_batch_4_vec, "GRM_batch_4.grm.bin", size = 4)
writeBin(rep(1, length(GRM_batch_4_vec)), "GRM_batch_4.grm.N.bin", size = 4)

# Write GRM files for Batch 5
GRM_batch_5 <- GRM[indices_batch_5, indices_batch_5]
ids_batch_5 <- phenotype[indices_batch_5, c("FID", "IID")]
GRM_batch_5_vec <- as.vector(GRM_batch_5)
write.table(ids_batch_5, file = "GRM_batch_5.grm.id", quote = FALSE, row.names = FALSE, col.names = FALSE)
writeBin(GRM_batch_5_vec, "GRM_batch_5.grm.bin", size = 4)
writeBin(rep(1, length(GRM_batch_5_vec)), "GRM_batch_5.grm.N.bin", size = 4)

# Write GRM files for Batch 6
GRM_batch_6 <- GRM[indices_batch_6, indices_batch_6]
ids_batch_6 <- phenotype[indices_batch_6, c("FID", "IID")]
GRM_batch_6_vec <- as.vector(GRM_batch_6)
write.table(ids_batch_6, file = "GRM_batch_6.grm.id", quote = FALSE, row.names = FALSE, col.names = FALSE)
writeBin(GRM_batch_6_vec, "GRM_batch_6.grm.bin", size = 4)
writeBin(rep(1, length(GRM_batch_6_vec)), "GRM_batch_6.grm.N.bin", size = 4)


                       
