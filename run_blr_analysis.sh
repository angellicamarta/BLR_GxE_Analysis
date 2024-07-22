#!/bin/bash
#SBATCH --job-name=blr_gxe_analysis  # Job name
#SBATCH --output=blr_gxe_%j.out      # Standard output and error log
#SBATCH --error=blr_gxe_%j.err       # Error log
#SBATCH --time=11:59:00              # Time limit hrs:min:sec
#SBATCH --partition=icelake-himem    # Partition to submit to
#SBATCH --ntasks=1                   # Number of tasks (processes)
#SBATCH --cpus-per-task=30           # Number of CPU cores per task (30 cores to match memory request)
#SBATCH --mem=200G                   # Total memory (200 GB)
#SBATCH --mail-type=START,END,FAIL   # Notifications for job done & fail
#SBATCH --mail-user=am3194@cam.ac.uk # Email to which notifications are sent
#SBATCH --account=BUTTERWORTH-SL3-CPU# Account name

module load R

# Step 1: Prepare data
Rscript prepare_believe_data.R

# Step 2: Run BLR analysis
Rscript blr_gxe_analysis.R
