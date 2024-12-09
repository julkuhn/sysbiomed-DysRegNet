import argparse
import pandas as pd
import sys


parser = argparse.ArgumentParser()

parser.add_argument('--input', type=str, required=True, help='Path to the expression matrix from gtex(genes as rows, samples as columns)')
parser.add_argument('--output', type=str, help="Output file for the filtered matrix", required=True)

args = parser.parse_args()


input_file = args.input
output_file = args.output

#read in matrix
tpm_gtex = pd.read_csv(input_file, sep = "\t", skiprows=2) 

#filter out gene descriptions that occur more than once
print('number of genes before removing duplicates ' + str(tpm_gtex.shape[0]))
duplicated_names = tpm_gtex['Description'].duplicated(keep=False)
print(f'{duplicated_names.sum()} genes names occured multiple times in the Description column')
tpm_gtex_no_dupl = tpm_gtex[~duplicated_names]
print('number of genes left after removing duplicates ' + str(tpm_gtex_no_dupl.shape[0]))

#filter out genes that are not expressed in at least 80% of the samples

non_expressed_counts = (tpm_gtex_no_dupl == 0).sum(axis=1)

#Calculate the proportion of non-expressed samples per gene
non_expressed_proportion = non_expressed_counts / tpm_gtex_no_dupl.shape[1]

#remove those genes
low_expression_genes = tpm_gtex_no_dupl.index[non_expressed_proportion >= 0.2]
tpm_gtex_filtered = tpm_gtex_no_dupl.drop(low_expression_genes, axis=0)

# Output the list of genes
print('number of genes before removing the not expressed: ' + str(tpm_gtex_no_dupl.shape[0]))
print(f"Number of genes not expressed in ≥80% of samples: {len(low_expression_genes)}")
print('number of genes after removing the not expressed: ' + str(tpm_gtex_filtered.shape[0]))

tpm_gtex_filtered.iloc[:, 1:].to_csv(output_file, index=False)
print(f'filtered expression matrix saved to {output_file}')