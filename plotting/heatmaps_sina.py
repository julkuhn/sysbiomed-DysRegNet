import argparse
import os
import pandas as pd

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


parser = argparse.ArgumentParser()
parser.add_argument('--input_gtex_lung', type=str, required=True)
parser.add_argument('--input_tcga_lung', type=str, required=True)
parser.add_argument('--input_gtex_breast', type=str, required=True)
parser.add_argument('--input_tcga_breast', type=str, required=True)
parser.add_argument('--output_dir', type=str, required=True)
parser.add_argument('--top_n_variant', type=int, required=False)  # option to set number of top n variant genes to be used
args = parser.parse_args()

gtex_data_lung = args.input_gtex_lung
tcga_data_lung = args.input_tcga_lung
gtex_data_breast = args.input_gtex_breast
tcga_data_breast = args.input_tcga_breast

output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)  # ensure the specified output_dir exists

#read in
gtex_lung = pd.read_csv(gtex_data_lung, sep=",")
tcga_lung = pd.read_csv(tcga_data_lung, sep=",")
gtex_breast = pd.read_csv(gtex_data_breast, sep=",")
tcga_breast = pd.read_csv(tcga_data_breast, sep=",")

# filter the dfs for their common genes
def common_genes(gtex_df, tcga_df):
    unique_gtex_genes = set(gtex_df['Description'])  # use unique descriptions and samples
    unique_tcga_samples = set(tcga_df['sample'])
    common = unique_gtex_genes.intersection(unique_tcga_samples)
    return common

common_genes_lung = common_genes(gtex_lung, tcga_lung)
common_genes_breast = common_genes(gtex_breast, tcga_breast)

gtex_lung_intersection = gtex_lung[gtex_lung['Description'].isin(common_genes_lung)]
tcga_lung_intersection = tcga_lung[tcga_lung['sample'].isin(common_genes_lung)]

gtex_breast_intersection = gtex_breast[gtex_breast['Description'].isin(common_genes_breast)]
tcga_breast_intersection = tcga_breast[tcga_breast['sample'].isin(common_genes_breast)]

def plot_heatmap(gtex_df, tcga_df, title, output_path=None, top_n=500):
    if 'Tissue' in gtex_df.columns:
        gtex_df = gtex_df.set_index(['Description', 'Tissue'])
    else:
        gtex_df = gtex_df.set_index('Description')

    if 'Tissue' in tcga_df.columns:
        tcga_df = tcga_df.set_index(['sample', 'Tissue'])
    else:
        tcga_df = tcga_df.set_index('sample')

    combined = pd.concat([gtex_df, tcga_df], axis=1)  # concatenate
    combined_log = combined.apply(lambda x: np.log1p(x))  # log

    if args.top_n_variant is not None:
        top_n = args.top_n_variant  # if number of top n variant genes specified in arguments, use this

    variances = combined_log.var(axis=1)  # variance for each gene
    top_variant_genes = variances.nlargest(top_n).index  # pick n genes with the highest variance
    most_variant = combined_log.loc[top_variant_genes]  # filter the most variant

    plt.figure(figsize=(14, 12))
    sns.heatmap(
        most_variant,
        cmap="viridis",
        cbar=True,
        xticklabels=False,  # remove sample names to avoid overlapping labels
        yticklabels=False  # remove gene names to avoid overlapping labels
    )
    plt.title(f"{title} (Top {top_n} variant genes, log-transformed)")
    plt.xlabel("Samples")
    plt.ylabel("Genes/Tissues" if 'Tissue' in gtex_df.index.names else "Genes")

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')

# heatmap for lung data only:
plot_heatmap(
    gtex_lung_intersection,
    tcga_lung_intersection,
    "Heatmap for Lung (GTEx and TCGA Intersections)",
    output_path=os.path.join(output_dir, "heatmap_lung.png")
)

# heatmap for breast data only:
plot_heatmap(
    gtex_breast_intersection,
    tcga_breast_intersection,
    "Heatmap for Breast (GTEx and TCGA Intersections)",
    output_path=os.path.join(output_dir, "heatmap_breast.png")
)

# filter the four dfs to the common genes for all
common_genes_all = common_genes_lung.intersection(common_genes_breast)
gtex_lung_common = gtex_lung_intersection[gtex_lung_intersection['Description'].isin(common_genes_all)]
tcga_lung_common = tcga_lung_intersection[tcga_lung_intersection['sample'].isin(common_genes_all)]
gtex_breast_common = gtex_breast_intersection[gtex_breast_intersection['Description'].isin(common_genes_all)]
tcga_breast_common = tcga_breast_intersection[tcga_breast_intersection['sample'].isin(common_genes_all)]

# combine  gtex lung and gtex breast & label the tissue
gtex_combined = pd.concat(
    [gtex_lung_common.assign(Tissue='Lung'), gtex_breast_common.assign(Tissue='Breast')],
    axis=0
)
# combine tcga lung and tcga breast & label the tissue
tcga_combined = pd.concat(
    [tcga_lung_common.assign(Tissue='Lung'), tcga_breast_common.assign(Tissue='Breast')],
    axis=0
)

# heatmap for lung and breast data together:
plot_heatmap(
    gtex_combined,
    tcga_combined,
    "Combined Heatmap (Lung and Breast Intersections)",
    output_path=os.path.join(output_dir, "heatmap_combined.png")
)
