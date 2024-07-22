# Load prepared data
load("prepared_believe_data.RData")

# Define interaction matrices for Albumin
X0_albumin <- matrix(0, nrow=nrow(X), ncol=ncol(X))
X_main_albumin <- rbind(X, X)
X_1_albumin <- rbind(X, X0_albumin)
X_2_albumin <- rbind(X0_albumin, X)

# Check dimensions
cat("Dimensions of X_main_albumin:", dim(X_main_albumin), "\n")
cat("Dimensions of X_1_albumin:", dim(X_1_albumin), "\n")
cat("Dimensions of X_2_albumin:", dim(X_2_albumin), "\n")

# Correct the groups parameter
groups_albumin <- rep(1:2, each = length(y) / 2)
cat("Length of groups_albumin:", length(groups_albumin), "\n")
cat("Length of y:", length(y), "\n")

# Run BGLR analysis for Albumin.env
fm_albumin <- BGLR(y = y, ETA = list(
  main = list(X = X_main_albumin, model = 'BRR'),
  int1 = list(X = X_1_albumin, model = 'BRR'),
  int2 = list(X = X_2_albumin, model = 'BRR')
), nIter = 6000, burnIn = 1000, saveAt = 'GxE_Albumin_', groups = groups_albumin)

# Post-processing for Albumin.env
varU_main_albumin <- scan('GxE_Albumin_ETA_main_varB.dat')[-c(1:1000)]
varU_int1_albumin <- scan('GxE_Albumin_ETA_int1_varB.dat')[-c(1:1000)]
varU_int2_albumin <- scan('GxE_Albumin_ETA_int2_varB.dat')[-c(1:1000)]
varE_albumin <- read.table('GxE_Albumin_varE.dat', header = FALSE)[-c(1:1000),]
varU1_albumin <- varU_main_albumin + varU_int1_albumin
varU2_albumin <- varU_main_albumin + varU_int2_albumin
h2_1_albumin <- varU1_albumin / (varU1_albumin + varE_albumin[, 1])
h2_2_albumin <- varU2_albumin / (varU2_albumin + varE_albumin[, 2])
COR_albumin <- varU_main_albumin / sqrt(varU1_albumin * varU2_albumin)

# Define interaction matrices for Vigorous10MinsActivity
X0_vigorous <- matrix(0, nrow=nrow(X), ncol=ncol(X))
X_main_vigorous <- rbind(X, X)
X_1_vigorous <- rbind(X, X0_vigorous)
X_2_vigorous <- rbind(X0_vigorous, X)

# Check dimensions
cat("Dimensions of X_main_vigorous:", dim(X_main_vigorous), "\n")
cat("Dimensions of X_1_vigorous:", dim(X_1_vigorous), "\n")
cat("Dimensions of X_2_vigorous:", dim(X_2_vigorous), "\n")

# Correct the groups parameter
groups_vigorous <- rep(1:2, each = length(y) / 2)
cat("Length of groups_vigorous:", length(groups_vigorous), "\n")
cat("Length of y:", length(y), "\n")

# Run BGLR analysis for Vigorous10MinsActivity.env
fm_vigorous <- BGLR(y = y, ETA = list(
  main = list(X = X_main_vigorous, model = 'BRR'),
  int1 = list(X = X_1_vigorous, model = 'BRR'),
  int2 = list(X = X_2_vigorous, model = 'BRR')
), nIter = 6000, burnIn = 1000, saveAt = 'GxE_Vigorous_', groups = groups_vigorous)

# Post-processing for Vigorous10MinsActivity.env
varU_main_vigorous <- scan('GxE_Vigorous_ETA_main_varB.dat')[-c(1:1000)]
varU_int1_vigorous <- scan('GxE_Vigorous_ETA_int1_varB.dat')[-c(1:1000)]
varU_int2_vigorous <- scan('GxE_Vigorous_ETA_int2_varB.dat')[-c(1:1000)]
varE_vigorous <- read.table('GxE_Vigorous_varE.dat', header = FALSE)[-c(1:1000),]
varU1_vigorous <- varU_main_vigorous + varU_int1_vigorous
varU2_vigorous <- varU_main_vigorous + varU_int2_vigorous
h2_1_vigorous <- varU1_vigorous / (varU1_vigorous + varE_vigorous[, 1])
h2_2_vigorous <- varU2_vigorous / (varU2_vigorous + varE_vigorous[, 2])
COR_vigorous <- varU_main_vigorous / sqrt(varU1_vigorous * varU2_vigorous)

# Define interaction matrices for Study environment
X0_study <- matrix(0, nrow=nrow(X), ncol=ncol(X))
X_main_study <- rbind(X, X)
X_1_study <- rbind(X, X0_study)
X_2_study <- rbind(X0_study, X)

# Check dimensions
cat("Dimensions of X_main_study:", dim(X_main_study), "\n")
cat("Dimensions of X_1_study:", dim(X_1_study), "\n")
cat("Dimensions of X_2_study:", dim(X_2_study), "\n")

# Correct the groups parameter
groups_study <- rep(1:2, each = length(y) / 2)
cat("Length of groups_study:", length(groups_study), "\n")
cat("Length of y:", length(y), "\n")

# Run BGLR analysis for study.env
fm_study <- BGLR(y = y, ETA = list(
  main = list(X = X_main_study, model = 'BRR'),
  int1 = list(X = X_1_study, model = 'BRR'),
  int2 = list(X = X_2_study, model = 'BRR')
), nIter = 6000, burnIn = 1000, saveAt = 'GxE_Study_', groups = groups_study)

# Post-processing for study.env
varU_main_study <- scan('GxE_Study_ETA_main_varB.dat')[-c(1:1000)]
varU_int1_study <- scan('GxE_Study_ETA_int1_varB.dat')[-c(1:1000)]
varU_int2_study <- scan('GxE_Study_ETA_int2_varB.dat')[-c(1:1000)]
varE_study <- read.table('GxE_Study_varE.dat', header = FALSE)[-c(1:1000),]
varU1_study <- varU_main_study + varU_int1_study
varU2_study <- varU_main_study + varU_int2_study
h2_1_study <- varU1_study / (varU1_study + varE_study[, 1])
h2_2_study <- varU2_study / (varU2_study + varE_study[, 2])
COR_study <- varU_main_study / sqrt(varU1_study * varU2_study)

# Save results
save(fm_albumin, h2_1_albumin, h2_2_albumin, COR_albumin, 
     fm_vigorous, h2_1_vigorous, h2_2_vigorous, COR_vigorous,
     fm_study, h2_1_study, h2_2_study, COR_study,
     file = "blr_gxe_results.RData")
