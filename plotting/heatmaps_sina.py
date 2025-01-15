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
parser.add_argument('--top_n_variant', type=int, default=500, required=False)  # option to set number of top n variant genes to be used
parser.add_argument('--show_gene_names', action='store_true', required=False)  # if set, args.show_gene_names will be True, otherwise it will be False
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

# Heatmaps without clustering:

def plot_heatmap(gtex_df, tcga_df, title, output_path=None, top_n=500, show_gene_names=False):
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

    variances = combined_log.var(axis=1)  # variance for each gene
    top_variant_genes = variances.nlargest(min(top_n, len(variances))).index  # pick n genes with the highest variance, without selecting more than available genes
    most_variant = combined_log.loc[top_variant_genes]  # filter the most variant

    plt.figure(figsize=(14, 12))
    sns.heatmap(
        most_variant,
        cmap="viridis",
        cbar=True,
        xticklabels=False,  # remove sample names
        yticklabels=show_gene_names  # show gene names if specified
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
    output_path=os.path.join(output_dir, "heatmap_lung.png"),
    top_n=args.top_n_variant,
    show_gene_names=args.show_gene_names
)

# heatmap for breast data only:
plot_heatmap(
    gtex_breast_intersection,
    tcga_breast_intersection,
    "Heatmap for Breast (GTEx and TCGA Intersections)",
    output_path=os.path.join(output_dir, "heatmap_breast.png"),
    top_n=args.top_n_variant,
    show_gene_names=args.show_gene_names
)

# filter the four dfs to the common genes for all
common_genes_all = common_genes_lung.intersection(common_genes_breast)
#print("Number of common genes for all 4 datasets: ", len(common_genes_all))
gtex_lung_common = gtex_lung_intersection[gtex_lung_intersection['Description'].isin(common_genes_all)]
tcga_lung_common = tcga_lung_intersection[tcga_lung_intersection['sample'].isin(common_genes_all)]
gtex_breast_common = gtex_breast_intersection[gtex_breast_intersection['Description'].isin(common_genes_all)]
tcga_breast_common = tcga_breast_intersection[tcga_breast_intersection['sample'].isin(common_genes_all)]

# combine  gtex lung and gtex breast
gtex_combined = pd.merge(gtex_lung_common, gtex_breast_common, on="Description", how="inner")

# combine tcga lung and tcga breast
tcga_combined = pd.merge(tcga_lung_common, tcga_breast_common, on="sample", how="inner")

#print("Shape of gtex lung and breast combined: ", gtex_combined.shape)
#print("Shape of tcga lung and breast combined: ", tcga_combined.shape)
#print("gtex_combined:")
#print(gtex_combined.head())
#print("tcga_combined:")
#print(tcga_combined.head())

# heatmap for lung and breast data together:
plot_heatmap(
    gtex_combined,
    tcga_combined,
    "Combined Heatmap (Lung and Breast Intersections)",
    output_path=os.path.join(output_dir, "heatmap_combined.png"),
    top_n=args.top_n_variant,
    show_gene_names=args.show_gene_names
)

# Seaborn clustermap:

def plot_clustermap(gtex_df, tcga_df, title, output_path=None, top_n=500, show_gene_names=False):
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

    variances = combined_log.var(axis=1)  # variance for each gene
    top_variant_genes = variances.nlargest(min(top_n, len(variances))).index  # pick n genes with the highest variance, without selecting more than available genes
    most_variant = combined_log.loc[top_variant_genes]  # filter the most variant

    # colors for gtex and tcga columns
    gtex_cols = ['red'] * gtex_df.shape[1]
    tcga_cols = ['blue'] * tcga_df.shape[1]
    col_colors = gtex_cols + tcga_cols

    sns.set(font_scale=1.0)
    cmap = sns.clustermap(
        most_variant,
        cmap="viridis",
        metric="euclidean",
        method="average",  # more balanced clusters that are not too sensitive to outliers
        figsize=(14, 12),
        xticklabels=False,  # remove sample names
        yticklabels=show_gene_names,  # show gene names if specified
        col_colors=col_colors  # assign colors for gtex and tcga columns
    )

    # legend for column colors
    legend_labels = ['GTEx', 'TCGA']
    legend_colors = ['red', 'blue']
    for label, color in zip(legend_labels, legend_colors):
        cmap.ax_col_dendrogram.bar(0, 0, color=color, label=label, linewidth=0)
    cmap.ax_col_dendrogram.legend(loc="center", ncol=2, bbox_to_anchor=(0.5, 1.1), title="Dataset")

    cmap.ax_heatmap.set_title(f"{title} (Top {top_n} variant genes, log-transformed)", fontsize=12)

    if output_path:
        cmap.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

# clustermap for lung data only:
plot_clustermap(
    gtex_lung_intersection,
    tcga_lung_intersection,
    "Clustermap for Lung (GTEx and TCGA Intersections)",
    output_path=os.path.join(output_dir, "clustermap_lung.png"),
    top_n=args.top_n_variant,
    show_gene_names=args.show_gene_names
)

# clustermap for breast data only:
plot_clustermap(
    gtex_breast_intersection,
    tcga_breast_intersection,
    "Clustermap for Breast (GTEx and TCGA Intersections)",
    output_path=os.path.join(output_dir, "clustermap_breast.png"),
    top_n=args.top_n_variant,
    show_gene_names=args.show_gene_names
)

# clustermap for lung and breast data together:
plot_clustermap(
    gtex_combined,
    tcga_combined,
    "Combined Clustermap (Lung and Breast Intersections)",
    output_path=os.path.join(output_dir, "clustermap_combined.png"),
    top_n=args.top_n_variant,
    show_gene_names=args.show_gene_names
)