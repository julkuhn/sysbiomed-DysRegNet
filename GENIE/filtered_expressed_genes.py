import argparse
import pandas as pd
import sys


parser = argparse.ArgumentParser()

parser.add_argument('--input', type=str, required=True, help='Path to the expression matrix from gtex(genes as rows, samples as columns)')
parser.add_argument('--output', type=str, help="Output file for the filtered matrix", required=True)

args = parser.parse_args()


input_file = args.input
output_file = args.output
#input_file = "/nfs/data2/dysregnet_gtex/GTEx_by_tissue//tissue_Blood.csv"
#output_file = "/nfs/home/students/f.gillhuber/preprocessing/R_GENIE3/test1.csv"

#read in matrix
tpm_gtex = pd.read_csv(input_file, sep = ",")#, skiprows=1) 

splitted = tpm_gtex['Combined'].str.split('_',n=1, expand=True)
splitted.columns = ['Names_{}'.format(x+1) for x in splitted.columns]

tpm_gtex = splitted.join(tpm_gtex)
tpm_gtex = tpm_gtex.rename(columns = {'Names_1':'Name', 'Names_2':'Description'})

tpm_gtex = tpm_gtex.drop(columns=['Combined', 'Unnamed: 0'])

#tpm_gtex.to_csv("/nfs/home/students/f.gillhuber/preprocessing/R_GENIE3/test1.csv")


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
tpm_gtex_no_ensemble = tpm_gtex_filtered[~(tpm_gtex_filtered['Description'].str.startswith('ENS'))]

# Output the list of genes
print('Number of genes before removing the not expressed: ' + str(tpm_gtex_no_dupl.shape[0]))
print(f"Number of genes not expressed in ≥80% of samples: {len(low_expression_genes)}")
print('Number of genes after removing the not expressed: ' + str(tpm_gtex_filtered.shape[0]))
print('Number of genes after removing the ones with ENSEMBLE IDs in the Description column ' + str(tpm_gtex_no_ensemble.shape[0]))

tpm_gtex_no_ensemble.iloc[:, 1:].to_csv(output_file, index=False)
print(f'filtered expression matrix saved to {output_file}')