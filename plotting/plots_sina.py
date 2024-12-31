import argparse
import os

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import seaborn as sns
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--input_gtex', type=str, required=True)
parser.add_argument('--input_tcga', type=str, required=True)
parser.add_argument('--output_dir', type=str, required=True)
args = parser.parse_args()

gtex_data = args.input_gtex
tcga_data = args.input_tcga
output_dir = args.output_dir

#read in
gtex_df = pd.read_csv(gtex_data, sep=",")
tcga_df = pd.read_csv(tcga_data, sep=",")

# filter the dfs for their common genes
unique_gtex_genes = set(gtex_df['Description'])  # use unique descriptions and samples
unique_tcga_samples = set(tcga_df['sample'])
print("Unique GTEx genes ('Description' column): ", len(unique_gtex_genes))
print("Unique TCGA genes ('sample' column): ", len(unique_tcga_samples))

common_genes = unique_gtex_genes.intersection(unique_tcga_samples)
print("Intersection of unique GTEx genes and unique TCGA genes: ", len(common_genes), " common genes")

gtex_intersection = gtex_df[gtex_df['Description'].isin(common_genes)]
tcga_intersection = tcga_df[tcga_df['sample'].isin(common_genes)]

print("Shape of GTEx intersection dataframe:")
print(gtex_intersection.shape)
print("Shape of TCGA intersection dataframe:")
print(tcga_intersection.shape)

# melt the dfs for distribution plot
gtex_melted = gtex_intersection.melt(id_vars=['Description'], var_name='GTEx_Sample', value_name='GTEx_Value')
tcga_melted = tcga_intersection.melt(id_vars=['sample'], var_name='TCGA_Sample', value_name='TCGA_Value')
# TODO maybe add .dropna() to the melting of the dfs

# Distribution plot with all values:
plt.figure(figsize=(10, 6))
sns.histplot(gtex_melted['GTEx_Value'], kde=True, color='red', label='GTEx', bins=30)
sns.histplot(tcga_melted['TCGA_Value'], kde=True, color='blue', label='TCGA', bins=30)
plt.title('Distribution of GTEx and TCGA values')
plt.xlabel('Expression value')
plt.ylabel('Density')
plt.legend()
# save file
os.makedirs(output_dir, exist_ok=True)  # ensure the specified output_dir exists
output_path = os.path.join(output_dir, "distribution_plot.png")
plt.savefig(output_path, dpi=300)
plt.close()

# Distribution plot without outliers:
def distribution_plot_without_outliers(gtex_dataframe, tcga_dataframe, outlier_threshold):
    gtex_filtered = gtex_dataframe.iloc[:, 1:].applymap(lambda x: x if x <= outlier_threshold else None)
    tcga_filtered = tcga_dataframe.iloc[:, 1:].applymap(lambda x: x if x <= outlier_threshold else None)
    melted_gtex = gtex_filtered.melt(value_name='Expression', var_name='Sample').dropna()
    melted_tcga = tcga_filtered.melt(value_name='Expression', var_name='Sample').dropna()

    plot, ax = plt.subplots(figsize=(10, 6))
    sns.histplot(melted_gtex['Expression'], bins=50, kde=True, color='red', label='GTEx', alpha=0.7, ax=ax)
    sns.histplot(melted_tcga['Expression'], bins=50, kde=True, color='blue', label='TCGA', alpha=0.7, ax=ax)
    ax.set_title(f'Distribution of GTEx and TCGA values with threshold at {outlier_threshold}')
    ax.set_xlabel('Expression value')
    ax.set_ylabel('Density')
    ax.legend()
    return plot

# create and save plots
for threshold in [1000, 200]:
    plot = distribution_plot_without_outliers(gtex_intersection, tcga_intersection, threshold)
    os.makedirs(output_dir, exist_ok=True)  # ensure the specified output_dir exists
    output_path = os.path.join(output_dir, f"distribution_plot_with_threshold_{threshold}.png")
    plot.savefig(output_path, dpi=300)
    plt.close(plot)

# Distribution plot with log values:
gtex_melted['GTEx_Value'] = np.log1p(gtex_melted['GTEx_Value'])  # log gtex (ln(1+x) to avoid log(0))
tcga_melted['TCGA_Value'] = np.log1p(tcga_melted['TCGA_Value'])  # log tcga
plt.figure(figsize=(10, 6))
sns.histplot(gtex_melted['GTEx_Value'], kde=True, color='red', label='GTEx log', bins=30)
sns.histplot(tcga_melted['TCGA_Value'], kde=True, color='blue', label='TCGA log', bins=30)
plt.title('Distribution of log-transformed GTEx and TCGA values')
plt.xlabel('Natural logarithmic expression value (log1p)')
plt.ylabel('Density')
plt.legend()
# save file
os.makedirs(output_dir, exist_ok=True)  # ensure the specified output_dir exists
output_path = os.path.join(output_dir, "distribution_plot_log_values.png")
plt.savefig(output_path, dpi=300)
plt.close()

# Distribution plot with log values and medians of the genes: -> use already logged gtex_melted and tcga_melted
gtex_medians = gtex_intersection.iloc[:, 1:].median(axis=1)  # gene-wise medians, skipping 'Description' column
tcga_medians = tcga_intersection.iloc[:, 1:].median(axis=1)  # gene-wise medians, skipping 'sample' column
gtex_medians = np.log1p(gtex_medians)  # log the medians
tcga_medians = np.log1p(tcga_medians)  # log the medians
print("GTEx medians:")
print(gtex_medians.head())
print("TCGA medians:")
print(tcga_medians.head())
plt.figure(figsize=(10, 6))
sns.histplot(gtex_medians, kde=True, color='red', label='GTEx logged gene medians', bins=30)
sns.histplot(tcga_medians, kde=True, color='blue', label='TCGA logged gene medians', bins=30)
plt.title('Distribution of Log-transformed Gene Medians')
plt.xlabel('Logarithmic Gene Medians (log1p)')
plt.ylabel('Density')
plt.legend()
# save file
os.makedirs(output_dir, exist_ok=True)  # ensure the specified output_dir exists
output_path = os.path.join(output_dir, "distribution_plot_log_gene_medians.png")
plt.savefig(output_path, dpi=300)
plt.close()

# Venn diagram: use original dfs (not intersected)
unique_gtex = set(gtex_df['Description'])  # take unique 'Description' values and unique 'sample' values
unique_tcga = set(tcga_df['sample'])
gtex_only = unique_gtex - unique_tcga  # genes in GTEx but not in TCGA
tcga_only = unique_tcga - unique_gtex  # genes in TCGA but not in GTEx
intersection = unique_gtex.intersection(unique_tcga)  # common genes

print("--VENN diagram--")
print("Exclusive GTEx genes, not in TCGA genes:", len(gtex_only))  # the ones in 'Description' and not in 'sample'
print("Intersection (common genes):", len(intersection))
print("Exclusive TCGA genes, not in GTEx genes:", len(tcga_only))  # the ones in 'sample' and not in 'Description'

plt.figure(figsize=(8, 8))
venn = venn2(
    subsets=(len(gtex_only), len(tcga_only), len(intersection)),  # attention to order: intersection as third value
    set_labels=('Exclusive GTEx genes', 'Exclusive TCGA genes', 'Common genes')
)
venn.get_patch_by_id('10').set_color('red')  # red for GTEx exclusive
venn.get_patch_by_id('01').set_color('blue')  # blue for TCGA exclusive
venn.get_patch_by_id('11').set_color('purple')  # purple for overlap
plt.title('Venn diagram of common genes and exclusive GTEx/TCGA genes')  # unique genes are used
os.makedirs(output_dir, exist_ok=True)  # ensure the output directory exists
output_path = os.path.join(output_dir, "venn_diagram.png")
plt.savefig(output_path, dpi=300)
plt.close()