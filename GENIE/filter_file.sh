#!/bin/bash


PYTHON_SCRIPT="/nfs/home/students/f.gillhuber/preprocessing/R_GENIE3/filtered_expressed_genes.py"

# Directory containing the split files
INPUT_DIR="/nfs/data2/dysregnet_gtex/GTEx_by_tissue/"
OUTPUT_DIR="/nfs/data2/dysregnet_gtex/GTEx_by_tissue_filtered/"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

for input_file in "$INPUT_DIR"/*.csv; do
    base_name=$(basename "$input_file" .csv)

    output_file="$OUTPUT_DIR/${base_name}_filtered.csv"

    python "$PYTHON_SCRIPT" --input "$input_file" --output "$output_file"

    echo "Filtered file saved to $output_file"
done