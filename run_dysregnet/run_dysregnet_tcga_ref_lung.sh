#!/bin/bash
#SBATCH --job-name=dysregnet
#SBATCH --output=dysregnet_%j.out
#SBATCH --error=dysregnet_%j.log
#SBATCH --time=24:00:00         # Set the time limit for the job
#SBATCH --ntasks=1              # Number of tasks
#SBATCH --cpus-per-task=4       # Number of CPU cores per task
#SBATCH --mem=32G               # Memory allocation
#SBATCH --mail-type=ALL               # Send email on job events
#SBATCH --mail-user=lara.lechner@tum.de  # Replace with your email address


# Activate the conda environment
   

# Assign paths to parameters
EXPR="/nfs/data/patients_networks/DysRegNet_workflow/results/expression_processed/LUAD/tpm.csv"
META="/nfs/data/patients_networks/DysRegNet_workflow/results/expression_processed/LUAD/tpm_meta.csv"
GRN="/nfs/data/patients_networks/DysRegNet_workflow/results/reference_networks/genie3/LUAD/genie3_LUAD_tpm.top_100k.csv"
OUT="/nfs/data2/dysregnet_gtex/results/tcga_lung_no_confounders_binary.fea"
OUT_STATS="/nfs/data2/dysregnet_gtex/results/tcga_lung_no_confounders.csv"

# Run the Python script with the parameters
python run_dysregnet.py --expr $EXPR --meta $META --grn $GRN --no_direction --output $OUT --output_stats $OUT_STATS