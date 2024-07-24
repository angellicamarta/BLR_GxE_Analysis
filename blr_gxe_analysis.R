library(BGLR)

# Define interaction matrices for Albumin
X0_albumin <- matrix(0, nrow=nrow(X), ncol=ncol(X))
X_main_albumin <- X
X_1_albumin <- X
X_2_albumin <- X

# Check dimensions
cat("Dimensions of X_main_albumin:", dim(X_main_albumin), "\n")
cat("Dimensions of X_1_albumin:", dim(X_1_albumin), "\n")
cat("Dimensions of X_2_albumin:", dim(X_2_albumin), "\n")

# Correct the groups parameter
groups_albumin <- rep(1, length(y))  # Each individual is their own group for now
cat("Length of groups_albumin:", length(groups_albumin), "\n")
cat("Length of y:", length(y), "\n")

# Run BGLR analysis for Albumin.env
fm_albumin <- BGLR(y = y, ETA = list(
  main = list(X = X_main_albumin, model = 'BRR'),
  int1 = list(X = X_1_albumin, model = 'BRR'),
  int2 = list(X = X_2_albumin, model = 'BRR')
), nIter = 6000, burnIn = 1000, saveAt = 'GxE_Albumin_', groups = groups_albumin)

# Read the file again, ensuring correct format
varE_albumin <- read.table('GxE_Albumin_varE.dat', header = FALSE)
# Exclude the first 1000 iterations (burn-in period)
varE_albumin <- varE_albumin[-c(1:1000), ]

# Ensure varE_albumin is a numeric vector
varE_albumin <- as.numeric(varE_albumin$V1)

# Check lengths to ensure consistency
cat("Length of varE_albumin:", length(varE_albumin), "\n")
cat("Length of varU1_albumin:", length(varU1_albumin), "\n")
cat("Length of varU2_albumin:", length(varU2_albumin), "\n")

# Calculating the variances
varU1_albumin <- varU_main_albumin + varU_int1_albumin
varU2_albumin <- varU_main_albumin + varU_int2_albumin

# Ensure varE_albumin has the correct length
if (length(varE_albumin) == length(varU1_albumin)) {
  # Calculating heritabilities
  h2_1_albumin <- varU1_albumin / (varU1_albumin + varE_albumin)
  h2_2_albumin <- varU2_albumin / (varU2_albumin + varE_albumin)

  # Calculating genetic correlation
  COR_albumin <- varU_main_albumin / sqrt(varU1_albumin * varU2_albumin)
  
  # Output results
  print(list(h2_1_albumin = h2_1_albumin, h2_2_albumin = h2_2_albumin, COR_albumin = COR_albumin))
} else {
  cat("Error: Length of varE_albumin does not match the lengths of varU1_albumin and varU2_albumin.\n")
}




# Define interaction matrices for Vigorous10MinsActivity
X0_vigorous <- matrix(0, nrow=nrow(X), ncol=ncol(X))
X_main_vigorous <- X
X_1_vigorous <- X
X_2_vigorous <- X

# Check dimensions
cat("Dimensions of X_main_vigorous:", dim(X_main_vigorous), "\n")
cat("Dimensions of X_1_vigorous:", dim(X_1_vigorous), "\n")
cat("Dimensions of X_2_vigorous:", dim(X_2_vigorous), "\n")

# Correct the groups parameter
groups_vigorous <- rep(1, length(y))  # Each individual is their own group for now
cat("Length of groups_vigorous:", length(groups_vigorous), "\n")
cat("Length of y:", length(y), "\n")

# Run BGLR analysis for Vigorous10MinsActivity.env
fm_vigorous <- BGLR(y = y, ETA = list(
  main = list(X = X_main_vigorous, model = 'BRR'),
  int1 = list(X = X_1_vigorous, model = 'BRR'),
  int2 = list(X = X_2_vigorous, model = 'BRR')
), nIter = 6000, burnIn = 1000, saveAt = 'GxE_Vigorous_', groups = groups_vigorous)

# Check if the analysis ran successfully
if (file.exists('GxE_Vigorous_varE.dat')) {
  # Post-processing for Vigorous10MinsActivity.env
  varU_main_vigorous <- scan('GxE_Vigorous_ETA_main_varB.dat')[-c(1:1000)]
  varU_int1_vigorous <- scan('GxE_Vigorous_ETA_int1_varB.dat')[-c(1:1000)]
  varU_int2_vigorous <- scan('GxE_Vigorous_ETA_int2_varB.dat')[-c(1:1000)]
  varE_vigorous <- scan('GxE_Vigorous_varE.dat', what = numeric(), skip = 1000)
  
  # Check lengths to ensure consistency
  cat("Length of varE_vigorous:", length(varE_vigorous), "\n")
  cat("Length of varU1_vigorous:", length(varU_main_vigorous), "\n")
  cat("Length of varU2_vigorous:", length(varU_int1_vigorous), "\n")
  
  # Calculating the variances
  varU1_vigorous <- varU_main_vigorous + varU_int1_vigorous
  varU2_vigorous <- varU_main_vigorous + varU_int2_vigorous
  
  # Ensure varE_vigorous has the correct length
  if (length(varE_vigorous) == length(varU1_vigorous)) {
    # Calculating heritabilities
    h2_1_vigorous <- varU1_vigorous / (varU1_vigorous + varE_vigorous)
    h2_2_vigorous <- varU2_vigorous / (varU2_vigorous + varE_vigorous)
    
    # Calculating genetic correlation
    COR_vigorous <- varU_main_vigorous / sqrt(varU1_vigorous * varU2_vigorous)
    
    # Output results
    print(list(h2_1_vigorous = h2_1_vigorous, h2_2_vigorous = h2_2_vigorous, COR_vigorous = COR_vigorous))
  } else {
    cat("Error: Length of varE_vigorous does not match the lengths of varU_main_vigorous.\n")
  }
} else {
  cat("Error: The BGLR analysis for Vigorous10MinsActivity.env did not complete successfully.\n")
}




# Define interaction matrices for Study Site
# Each study site is treated as a separate environmental variable
X0_study <- matrix(0, nrow=nrow(X), ncol=ncol(X))

# Interaction matrices for each study site
X_main_study <- X
X_1_study <- X0_study
X_2_study <- X0_study

# Assuming the study site columns are the last 4 columns in X
for (i in (ncol(X) - 3):ncol(X)) {
  X_1_study[, i] <- X[, i]
  X_2_study[, i] <- X[, i]
}

# Check dimensions
cat("Dimensions of X_main_study:", dim(X_main_study), "\n")
cat("Dimensions of X_1_study:", dim(X_1_study), "\n")
cat("Dimensions of X_2_study:", dim(X_2_study), "\n")

# Correct the groups parameter
groups_study <- rep(1, length(y))  # Each individual is their own group for now
cat("Length of groups_study:", length(groups_study), "\n")
cat("Length of y:", length(y), "\n")

# Run BGLR analysis for Study Site
fm_study <- BGLR(y = y, ETA = list(
  main = list(X = X_main_study, model = 'BRR'),
  int1 = list(X = X_1_study, model = 'BRR'),
  int2 = list(X = X_2_study, model = 'BRR')
), nIter = 6000, burnIn = 1000, saveAt = 'GxE_Study_', groups = groups_study)

# Check if the analysis ran successfully
if (file.exists('GxE_Study_varE.dat')) {
  # Post-processing for Study Site
  varU_main_study <- scan('GxE_Study_ETA_main_varB.dat')[-c(1:1000)]
  varU_int1_study <- scan('GxE_Study_ETA_int1_varB.dat')[-c(1:1000)]
  varU_int2_study <- scan('GxE_Study_ETA_int2_varB.dat')[-c(1:1000)]
  varE_study <- read.table('GxE_Study_varE.dat', header = FALSE)[-c(1:1000), ]

  # Verify the structure of varE_study
  str(varE_study)
  
  # Ensure varE_study is a numeric vector
  varE_study <- as.numeric(varE_study)
  
  # Check lengths to ensure consistency
  cat("Length of varE_study:", length(varE_study), "\n")
  cat("Length of varU1_study:", length(varU_main_study), "\n")
  cat("Length of varU2_study:", length(varU_int1_study), "\n")
  
  # Calculating the variances
  varU1_study <- varU_main_study + varU_int1_study
  varU2_study <- varU_main_study + varU_int2_study
  
  # Ensure varE_study has the correct length
  if (length(varE_study) == length(varU1_study)) {
    # Calculating heritabilities
    h2_1_study <- varU1_study / (varU1_study + varE_study)
    h2_2_study <- varU2_study / (varU2_study + varE_study)
    
    # Calculating genetic correlation
    COR_study <- varU_main_study / sqrt(varU1_study * varU2_study)
    
    # Output results
    print(list(h2_1_study = h2_1_study, h2_2_study = h2_2_study, COR_study = COR_study))
  } else {
    cat("Error: Length of varE_study does not match the lengths of varU_main_study.\n")
  }
} else {
  cat("Error: The BGLR analysis for Study Site did not complete successfully.\n")
}

