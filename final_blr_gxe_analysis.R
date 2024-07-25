library(BGLR)

# Read GRM data using the gap package
GRM <- with(ReadGRMBin("/home/am3194/rds/hpc-work/gcta_believe/merged_grm2"), GRM)

# Function to filter missing data
filter_missing_data <- function(...) {
  datasets <- list(...)
  complete_cases <- complete.cases(do.call(cbind, datasets))
  return(lapply(datasets, function(x) x[complete_cases, , drop = FALSE]))
}

# Prepare phenotype and covariate data
phenotype <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/BELIEVE.pheno", header = FALSE, col.names = c("FID", "IID", "T2D_status"))
albumin_discrete <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/Albumin_discrete_file.env", header = FALSE, col.names = c("FID", "IID", "Albumin"))
vigorous_activity <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/Vigorous10MinsActivity.env", header = FALSE, col.names = c("FID", "IID", "Vigorous10MinsActivity"))
categorical_covariates <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/C.covar", header = FALSE, col.names = c("FID", "IID", "Sex"))
quantitative_covariates <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/Q.covar", header = FALSE, col.names = c("FID", "IID", "Age", "BMI"))
study_discrete <- read.table("/home/am3194/rds/hpc-work/gcta_believe/pheno_data/study_discrete.env", header = FALSE, col.names = c("FID", "IID", "Study"))

# Filter missing data from all datasets
filtered_data <- filter_missing_data(phenotype, albumin_discrete, vigorous_activity, categorical_covariates, quantitative_covariates, study_discrete)
phenotype <- filtered_data[[1]]
albumin_discrete <- filtered_data[[2]]
vigorous_activity <- filtered_data[[3]]
categorical_covariates <- filtered_data[[4]]
quantitative_covariates <- filtered_data[[5]]
study_discrete <- filtered_data[[6]]

# Ensure GRM matches the number of rows in filtered datasets
GRM <- GRM[complete.cases(cbind(phenotype$T2D_status, albumin_discrete$Albumin, vigorous_activity$Vigorous10MinsActivity, categorical_covariates$Sex, quantitative_covariates$Age, quantitative_covariates$BMI, study_discrete$Study)), ]

# Prepare covariates by converting factors to dummy variables
X_cov <- model.matrix(~ Sex + Age + BMI - 1, data = cbind(categorical_covariates, quantitative_covariates))

# Function to run the BLR model for a given environmental variable
run_blr <- function(env_name, env_var, save_prefix) {
  # Create model matrix for the environmental variable
  env_var_matrix <- model.matrix(~ env_var - 1)
  
  # Ensure the dimensions of env_var_matrix match the number of individuals in GRM
  if (nrow(env_var_matrix) != nrow(GRM)) {
    stop("The number of rows in env_var_matrix does not match the number of rows in GRM.")
  }
  
  # Interaction term (gene x environment)
  X_gxe <- as.matrix(GRM) %*% env_var_matrix
  
  # Define the BLR model
  fm <- BGLR(y = y, ETA = list(
    main = list(X = as.matrix(GRM), model = 'BRR'),
    env = list(X = env_var_matrix, model = 'BRR'),
    gxe = list(X = X_gxe, model = 'BRR'),
    cov = list(X = as.matrix(X_cov), model = 'FIXED')
  ), nIter = 6000, burnIn = 1000, saveAt = save_prefix, groups = rep(1, length(y)))
  
  # Post-processing
  if (file.exists(paste0(save_prefix, 'varE.dat'))) {
    varU_main <- scan(paste0(save_prefix, 'ETA_main_varB.dat'))[-c(1:1000)]
    varU_env <- scan(paste0(save_prefix, 'ETA_env_varB.dat'))[-c(1:1000)]
    varU_gxe <- scan(paste0(save_prefix, 'ETA_gxe_varB.dat'))[-c(1:1000)]
    varE <- scan(paste0(save_prefix, 'varE.dat'), what = numeric(), skip = 1000)
    
    # Check lengths to ensure consistency
    cat("Length of varE_", env_name, ":", length(varE), "\n")
    cat("Length of varU1_", env_name, ":", length(varU_env), "\n")
    cat("Length of varU2_", env_name, ":", length(varU_gxe), "\n")
    
    # Calculating the variances
    varU1 <- varU_main + varU_env
    varU2 <- varU_main + varU_gxe
    
    # Ensure varE has the correct length
    if (length(varE) == length(varU1)) {
      # Calculating heritabilities
      h2_1 <- varU1 / (varU1 + varE)
      h2_2 <- varU2 / (varU2 + varE)
      
      # Calculating genetic correlation
      COR <- varU_main / sqrt(varU1 * varU2)
      
      # Output results
      print(list(h2_1 = h2_1, h2_2 = h2_2, COR = COR))
    } else {
      cat("Error: Length of varE_", env_name, " does not match the lengths of varU1 and varU2.\n")
    }
  } else {
    cat("Error: The BGLR analysis for ", env_name, " did not complete successfully.\n")
  }
}

# Prepare y after filtering
y <- phenotype$T2D_status

# Run the model for each environmental variable

# Albumin
run_blr(env_name = "Albumin", env_var = albumin_discrete$Albumin, save_prefix = 'GxE_Albumin_')

# Vigorous10MinsActivity
run_blr(env_name = "Vigorous10MinsActivity", env_var = vigorous_activity$Vigorous10MinsActivity, save_prefix = 'GxE_Vigorous_')

# Study
run_blr(env_name = "Study", env_var = study_discrete$Study, save_prefix = 'GxE_Study_')
