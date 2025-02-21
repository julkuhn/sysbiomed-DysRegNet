#!/bin/bash
#SBATCH --job-name=genie3_pipeline_job         
#SBATCH --output=/nfs/home/students/f.gillhuber/preprocessing/R_GENIE3/output_GENIE3/genie3_pipeline_output.%j.log 
#SBATCH --error=/nfs/home/students/f.gillhuber/preprocessing/R_GENIE3/output_GENIE3/genie3_pipeline_error.%j.log  
#SBATCH --ntasks=1                    
#SBATCH --cpus-per-task=50           
#SBATCH --mem=64G                     
#SBATCH --time=95:00:00               
#SBATCH --mail-type=ALL               
#SBATCH --mail-user=f.gillhuber@tum.de 
#SBATCH --array=0-30

# Set variables for Python filtering
PYTHON_SCRIPT="/nfs/home/students/f.gillhuber/preprocessing/R_GENIE3/filtered_expressed_genes.py"
INPUT_DIR="/nfs/data2/dysregnet_gtex/GTEx_by_tissue/"
OUTPUT_DIR="/nfs/data2/dysregnet_gtex/GTEx_by_tissue_filtered/"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run Python filtering on all files in the directory
echo "Starting Python filtering..."
for input_file in "$INPUT_DIR"/*.csv; do
    base_name=$(basename "$input_file" .csv)
    output_file="$OUTPUT_DIR/${base_name}_filtered.csv"

    python "$PYTHON_SCRIPT" --input "$input_file" --output "$output_file"

    # Check if the Python script ran successfully for each file
    if [ $? -ne 0 ]; then
        echo "Python filtering failed for $input_file. Skipping this file."
        continue
    fi

    echo "Filtered file saved to $output_file"
done

INPUT_DIR="/nfs/data2/dysregnet_gtex/GTEx_by_tissue_filtered/"
OUTPUT_DIR="/nfs/data2/dysregnet_gtex/results/"
NORMALIZE="TRUE"
REGULATOR_FILE="/nfs/home/students/f.gillhuber/preprocessing/R_GENIE3/input_GENIE3/restriction.txt"

# Create a list of files and get the file corresponding to the job array index
file_list=($INPUT_DIR/*.csv)
input_file="${file_list[${SLURM_ARRAY_TASK_ID}]}"

# Check if the input file is valid
if [ ! -f "$input_file" ]; then
    echo "Error: Input file not found for task ID ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Starting GENIE3 for file: $input_file"

# Check if the regulatory file is defined and run the R script accordingly
if [ -n "$REGULATOR_FILE" ]; then
    echo "Using regulator file: $REGULATOR_FILE"
    Rscript genie3_script.R "$input_file" "$REGULATOR_FILE" "$NORMALIZE"
else
    echo "No regulator file provided. Running without it."
    Rscript genie3_script.R "$input_file" "" "$NORMALIZE"
fi   

# Check if the R script ran successfully
if [ $? -ne 0 ]; then
    echo "GENIE3 script failed for input file $input_file. Exiting job."
    exit 1
fi

echo "Finished processing $input_file"
echo "Pipeline completed successfully!"