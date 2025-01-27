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
EXPR="/nfs/data/patients_networks/DysRegNet_workflow/results/expression_processed/BRCA/tpm.csv"
META="/nfs/data/patients_networks/DysRegNet_workflow/results/expression_processed/BRCA/tpm_meta.csv"
#GRN="/nfs/data/patients_networks/DysRegNet_workflow/results/reference_networks/genie3/BRCA/genie3_BRCA_tpm.top_100k.csv"
#Felix network
GRN="/nfs/data2/dysregnet_gtex/GENIE3_output/linkedList_output_slurm_tcga_breast_controls.csv"
#OUT="/nfs/data2/dysregnet_gtex/results/tcga_breast_no_confounders_binary.fea"
#OUT_STATS="/nfs/data2/dysregnet_gtex/results/tcga_breast_no_confounders.csv"
#out felix
OUT="/nfs/data2/dysregnet_gtex/results/tcga_refNets_Felix/tcga_breast_felix.fea"
OUT_STATS="/nfs/data2/dysregnet_gtex/results/tcga_refNets_Felix/tcga_breast_Felix.csv"

# Run the Python script with the parameters
python run_dysregnet.py --expr $EXPR --meta $META --grn $GRN --no_direction --output $OUT --output_stats $OUT_STATS