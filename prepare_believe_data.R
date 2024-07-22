# Load necessary libraries
library(gap)
library(BGLR)

# Read GRM data using the gap package
GRM <- with(ReadGRMBin("/home/am3194/rds/hpc-work/gcta_believe/merged_grm2"), GRM)

# Read phenotype and covariate data without headers and assign column names
phenotype <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/BELIEVE.pheno", header = FALSE, col.names = c("FID", "IID", "T2D_status"))
albumin <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/Albumin.env", header = FALSE, col.names = c("FID", "IID", "Albumin"))
vigorous_activity <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/Vigorous10MinsActivity.env", header = FALSE, col.names = c("FID", "IID", "Vigorous10MinsActivity"))
categorical_covariates <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/C.covar", header = FALSE, col.names = c("FID", "IID", "Sex"))
quantitative_covariates <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/Q.covar", header = FALSE, col.names = c("FID", "IID", "Age", "BMI"))
study_env <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/study.env", header = FALSE, col.names = c("FID", "IID", "study_BELIEVE", "study_BELURBAN", "study_BELSLUM", "study_BELRURAL"))

# Prepare phenotype and covariate matrices
y <- phenotype$T2D_status
PA <- as.numeric(vigorous_activity$Vigorous10MinsActivity)
BMI <- as.numeric(quantitative_covariates$BMI)
Age <- as.numeric(quantitative_covariates$Age)
Sex <- as.numeric(categorical_covariates$Sex)
Albumin <- as.numeric(albumin$Albumin)
Study <- as.matrix(study_env[, 3:6]) # study columns

# Combine all covariates including PA into the matrix X
X <- cbind(Sex, Age, BMI, Albumin, PA, Study)

# Save data for further analysis
save(GRM, y, X, file = "prepared_believe_data.RData")
