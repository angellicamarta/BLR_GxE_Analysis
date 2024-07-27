# Load necessary libraries
library(gap)

# Function to write the GRM matrix in binary format using writeBin
write_grm_bin <- function(grm, prefix, ids) {
  n <- nrow(grm)
  # Ensure the GRM matrix is in vector form (upper triangular matrix)
  grm_vec <- as.vector(grm[upper.tri(grm, diag = TRUE)])
  grm_n <- rep(n, length(grm_vec))
  
  # Write GRM id file
  write.table(ids, file = paste0(prefix, ".grm.id"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Write GRM binary file
  writeBin(grm_vec, paste0(prefix, ".grm.bin"), size = 4)
  
  # Write GRM N binary file
  writeBin(grm_n, paste0(prefix, ".grm.N.bin"), size = 4)
}

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

# Define the function to process each batch
process_batch <- function(indices, batch_number) {
  grm_batch <- GRM[indices, indices]
  ids_batch <- phenotype[indices, c("FID", "IID")]
  
  # Write GRM files for current batch
  write_grm_bin(grm_batch, paste0("GRM_batch_", batch_number), ids_batch)
}

# Process each batch
process_batch(indices_batch_1, 1)
process_batch(indices_batch_2, 2)
process_batch(indices_batch_3, 3)
process_batch(indices_batch_4, 4)
process_batch(indices_batch_5, 5)
process_batch(indices_batch_6, 6)

# Verification function with detailed diagnostics
check_grm_files <- function(prefix, original_grm, n_individuals) {
  ids <- read.table(paste0(prefix, ".grm.id"), header = FALSE)
  cat("Number of IDs read:", nrow(ids), "\n")
  
  if (nrow(ids) != n_individuals) {
    stop("Number of IDs does not match the number of individuals in the original GRM.")
  }
  
  grm_vec <- readBin(paste0(prefix, ".grm.bin"), what = numeric(), n = n_individuals * (n_individuals + 1) / 2, size = 4)
  cat("Length of GRM vector read:", length(grm_vec), "\n")
  
  grm_n <- readBin(paste0(prefix, ".grm.N.bin"), what = integer(), n = n_individuals * (n_individuals + 1) / 2, size = 4)
  cat("Length of GRM N vector read:", length(grm_n), "\n")
  
  cat("First few values from the original GRM vector:\n")
  print(as.vector(original_grm[upper.tri(original_grm, diag = TRUE)])[1:10])
  
  cat("First few values from the GRM vector read from binary file:\n")
  print(grm_vec[1:10])
  
  cat("Comparing with the original GRM vector...\n")
  comparison_result <- identical(as.vector(original_grm[upper.tri(original_grm, diag = TRUE)]), grm_vec)
  
  if (!comparison_result) {
    cat("Difference between original GRM vector and read GRM vector:\n")
    print(as.vector(original_grm[upper.tri(original_grm, diag = TRUE)])[1:10] - grm_vec[1:10])
  }
  
  return(comparison_result)
}

# Example usage for batch 1
prefix <- "GRM_batch_1"
original_grm <- GRM[indices_batch_1, indices_batch_1]
n_individuals <- nrow(original_grm)
check_grm_files(prefix, original_grm, n_individuals)
