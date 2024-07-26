load("prepared_believe_data.RData")
rm(X)

# Install required packages if not already installed
required_packages <- c("BGLR", "data.table")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

library(BGLR)
library(data.table)

# Prepare phenotype and covariate data
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
                       
# Ensure Age and BMI are numeric
quantitative_covariates$Age <- as.numeric(quantitative_covariates$Age)
quantitative_covariates$BMI <- as.numeric(quantitative_covariates$BMI)

# Ensure other covariates are factors
categorical_covariates$Sex <- as.factor(categorical_covariates$Sex)
albumin_discrete$Albumin <- as.factor(albumin_discrete$Albumin)
vigorous_activity$Vigorous10MinsActivity <- as.factor(vigorous_activity$Vigorous10MinsActivity)
study_discrete$Study <- as.factor(study_discrete$Study)

# Prepare covariates by converting factors to dummy variables
X_cov <- model.matrix(~ Sex + Age + BMI - 1, data = data.frame(categorical_covariates, quantitative_covariates))

# Prepare y after filtering
y <- phenotype$T2D_status

# Study
env_var_study <- study_discrete$Study
env_var_matrix_study <- model.matrix(~ env_var_study - 1)

# Check dimensions
cat("Dimensions of env_var_matrix_study:", dim(env_var_matrix_study), "\n")

if (nrow(env_var_matrix_study) != nrow(GRM)) {
  stop("The number of rows in env_var_matrix_study does not match the number of rows in GRM.")
}

# Interaction term (gene x environment) for Study
X_gxe_study <- as.matrix(GRM) %*% as.matrix(env_var_matrix_study)

# Define the BLR model for Study
fm_study <- BGLR(y = y, ETA = list(
  main = list(X = as.matrix(GRM), model = 'BRR'),
  env = list(X = as.matrix(env_var_matrix_study), model = 'BRR'),
  gxe = list(X = X_gxe_study, model = 'BRR'),
  cov = list(X = as.matrix(X_cov), model = 'FIXED')
), nIter = 6000, burnIn = 1000, saveAt = 'GxE_Study_', groups = rep(1, length(y)))

# Post-processing for Study
if (file.exists('GxE_Study_varE.dat')) {
  varU_main_study <- scan('GxE_Study_ETA_main_varB.dat')[-c(1:1000)]
  varU_env_study <- scan('GxE_Study_ETA_env_varB.dat')[-c(1:1000)]
  varU_gxe_study <- scan('GxE_Study_ETA_gxe_varB.dat')[-c(1:1000)]
  varE_study <- scan('GxE_Study_varE.dat', what = numeric(), skip = 1000)
  
  # Check lengths to ensure consistency
  cat("Length of varE_Study:", length(varE_study), "\n")
  cat("Length of varU1_Study:", length(varU_env_study), "\n")
  cat("Length of varU2_Study:", length(varU_gxe_study), "\n")
  
  # Calculating the variances
  varU1_study <- varU_main_study + varU_env_study
  varU2_study <- varU_main_study + varU_gxe_study
  
  # Ensure varE has the correct length
  if (length(varE_study) == length(varU1_study)) {
    # Calculating heritabilities
    h2_1_study <- varU1_study / (varU1_study + varE_study)
    h2_2_study <- varU2_study / (varU2_study + varE_study)
    
    # Calculating genetic correlation
    COR_study <- varU_main_study / sqrt(varU1_study * varU2_study)
    
    # Save results to a CSV file
    results_study <- data.frame(h2_1 = h2_1_study, h2_2 = h2_2_study, COR = COR_study)
    fwrite(results_study, 'GxE_Study_results.csv')
    
    # Output results
    print(list(h2_1 = h2_1_study, h2_2 = h2_2_study, COR = COR_study))
  } else {
    cat("Error: Length of varE_Study does not match the lengths of varU1 and varU2.\n")
  }
} else {
  cat("Error: The BGLR analysis for Study did not complete successfully.\n")
}
                       
