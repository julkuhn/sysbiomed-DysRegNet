#!/bin/bash
#SBATCH --job-name=genie3_job          # Job name
#SBATCH --output=/nfs/home/students/f.gillhuber/preprocessing/R_GENIE3/output_GENIE3/genie3_output.%j.log # Output file (%j will be replaced by the job ID)
#SBATCH --error=/nfs/home/students/f.gillhuber/preprocessing/R_GENIE3/output_GENIE3/genie3_error.%j.log   # Error file
#SBATCH --ntasks=1                    # Number of tasks (set to 1 for serial jobs)
#SBATCH --cpus-per-task=50           # Number of cores for this job
#SBATCH --mem=128G                     # Memory allocation
#SBATCH --time=24:00:00               # Time limit (HH:MM:SS)
#SBATCH --mail-type=ALL               # Send email on job events
#SBATCH --mail-user=f.gillhuber@tum.de  # Replace with your email address

# Set variables
#EXPR_MATRIX="/nfs/data2/dysregnet_gtex/gtex/gene_tpm_v10_breast_mammary_tissu_filtered.csv" # Path to expression matrix
#EXPR_MATRIX="/nfs/data2/dysregnet_gtex/gtex/gene_tpm_v10_lung_filtered.csv" # Path to expression matrix
EXPR_MATRIX="/nfs/data2/dysregnet_gtex/gtex/tcga_lung_controls.csv"
#EXPR_MATRIX="/nfs/data2/dysregnet_gtex/GTEx_by_tissue_filtered/tissue_Brain_filtered.csv"
REGULATOR_FILE="/nfs/home/students/f.gillhuber/preprocessing/R_GENIE3/input_GENIE3/restriction.txt" # Path to regulators (optional)
NORMALIZE="TRUE"                               # TRUE or FALSE (optional)

# Run the R script
Rscript genie3_script.R $EXPR_MATRIX $REGULATOR_FILE $NORMALIZE
