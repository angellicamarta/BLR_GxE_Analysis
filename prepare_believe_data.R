# Load necessary libraries
library(gap)
library(BGLR)

# Read GRM data using the gap package
GRM <- with(ReadGRMBin("/home/am3194/rds/hpc-work/gcta_believe/merged_grm2"), GRM)

# Read phenotype and covariate data without headers and assign column names
phenotype <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/BELIEVE.pheno", header = FALSE, col.names = c("FID", "IID", "T2D_status"))
albumin_discrete <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/Albumin_discrete_file.env", header = FALSE, col.names = c("FID", "IID", "Albumin"))
vigorous_activity <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/Vigorous10MinsActivity.env", header = FALSE, col.names = c("FID", "IID", "Vigorous10MinsActivity"))
categorical_covariates <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/C.covar", header = FALSE, col.names = c("FID", "IID", "Sex"))
quantitative_covariates <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/Q.covar", header = FALSE, col.names = c("FID", "IID", "Age", "BMI"))
study_discrete <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/study_discrete.env", header = FALSE, col.names = c("FID", "IID", "Study"))

# Prepare phenotype and covariate matrices
y <- phenotype$T2D_status

# Remove missing data coded as -9 and corresponding rows in covariates
valid_indices <- y != "-9"
y <- y[valid_indices]
vigorous_activity <- vigorous_activity[valid_indices, ]
quantitative_covariates <- quantitative_covariates[valid_indices, ]
categorical_covariates <- categorical_covariates[valid_indices, ]
albumin_discrete <- albumin_discrete[valid_indices, ]
study_discrete <- study_discrete[valid_indices, ]

# Convert necessary variables to factors
y <- as.factor(y)
PA <- as.numeric(vigorous_activity$Vigorous10MinsActivity)
BMI <- as.numeric(quantitative_covariates$BMI)
Age <- as.numeric(quantitative_covariates$Age)
Sex <- as.factor(categorical_covariates$Sex)
Albumin <- as.factor(albumin_discrete$Albumin)
Study <- as.factor(study_discrete$Study)

# Combine all covariates including PA into the matrix X
X <- cbind(Sex, Age, BMI, Albumin, PA, Study)

# Save data for further analysis
save(GRM, y, X, file = "prepared_believe_data.RData")
