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

# Albumin
env_var_albumin <- albumin_discrete$Albumin
env_var_matrix_albumin <- model.matrix(~ env_var_albumin - 1)
                       
# Interaction term (gene x environment) for Albumin
X_gxe_albumin <- as.matrix(GRM) %*% as.matrix(env_var_matrix_albumin)
                       
# Define the BLR model for Albumin
fm_albumin <- BGLR(y = y, ETA = list(
  main = list(X = as.matrix(GRM), model = 'BRR'),
  env = list(X = as.matrix(env_var_matrix_albumin), model = 'BRR'),
  gxe = list(X = X_gxe_albumin, model = 'BRR'),
  cov = list(X = as.matrix(X_cov), model = 'FIXED')
), nIter = 6000, burnIn = 1000, saveAt = 'GxE_Albumin_', groups = rep(1, length(y)))

# Post-processing for Albumin
if (file.exists('GxE_Albumin_varE.dat')) {
  varU_main_albumin <- scan('GxE_Albumin_ETA_main_varB.dat')[-c(1:1000)]
  varU_env_albumin <- scan('GxE_Albumin_ETA_env_varB.dat')[-c(1:1000)]
  varU_gxe_albumin <- scan('GxE_Albumin_ETA_gxe_varB.dat')[-c(1:1000)]
  varE_albumin <- scan('GxE_Albumin_varE.dat', what = numeric(), skip = 1000)
  
  # Check lengths to ensure consistency
  cat("Length of varE_Albumin:", length(varE_albumin), "\n")
  cat("Length of varU1_Albumin:", length(varU_env_albumin), "\n")
  cat("Length of varU2_Albumin:", length(varU_gxe_albumin), "\n")
  
  # Calculating the variances
  varU1_albumin <- varU_main_albumin + varU_env_albumin
  varU2_albumin <- varU_main_albumin + varU_gxe_albumin
  
  # Ensure varE has the correct length
  if (length(varE_albumin) == length(varU1_albumin)) {
    # Calculating heritabilities
    h2_1_albumin <- varU1_albumin / (varU1_albumin + varE_albumin)
    h2_2_albumin <- varU2_albumin / (varU2_albumin + varE_albumin)
    
    # Calculating genetic correlation
    COR_albumin <- varU_main_albumin / sqrt(varU1_albumin * varU2_albumin)
    
    # Save results to a CSV file
    results_albumin <- data.frame(h2_1 = h2_1_albumin, h2_2 = h2_2_albumin, COR = COR_albumin)
    fwrite(results_albumin, 'GxE_Albumin_results.csv')
    
    # Output results
    print(list(h2_1 = h2_1_albumin, h2_2 = h2_2_albumin, COR = COR_albumin))
  } else {
    cat("Error: Length of varE_Albumin does not match the lengths of varU1 and varU2.\n")
  }
} else {
  cat("Error: The BGLR analysis for Albumin did not complete successfully.\n")
}
                
