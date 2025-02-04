import argparse
import os
import pandas as pd

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

plot_colors = {
    'GTEX': 'lightcoral',  # Light red
    'TCGA': 'lightblue',  # Light blue
    'lung': 'yellow',
    'breast': 'orange',
}

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

'''
print("gtex_lung_intersection contains NaN? ", gtex_lung_intersection.isna().any().any())
print("gtex_breast_intersection contains NaN? ", gtex_breast_intersection.isna().any().any())
print("tcga_lung_intersection contains NaN? ", tcga_lung_intersection.isna().any().any())
print("tcga_breast_intersection contains NaN? ", tcga_breast_intersection.isna().any().any())
'''

# Heatmaps without clustering:

def plot_heatmap(gtex_df, tcga_df, title, output_path=None, top_n=500, show_gene_names=False):
    gtex_df = gtex_df.set_index('Description')
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

'''
print("gtex_lung_common contains NaN? ", gtex_lung_common.isna().any().any())
print("tcga_lung_common contains NaN? ", tcga_lung_common.isna().any().any())
print("gtex_breast_common contains NaN? ", gtex_breast_common.isna().any().any())
print("tcga_breast_common contains NaN? ", tcga_breast_common.isna().any().any())
'''

# combine  gtex lung and gtex breast
gtex_combined = pd.merge(gtex_lung_common, gtex_breast_common, on="Description", how="inner")
# combine tcga lung and tcga breast
tcga_combined = pd.merge(tcga_lung_common, tcga_breast_common, on="sample", how="inner")

# heatmap for lung and breast data together:
plot_heatmap(
    gtex_combined,
    tcga_combined,
    "Combined Heatmap (Lung and Breast Intersections)",
    output_path=os.path.join(output_dir, "heatmap_combined.png"),
    top_n=args.top_n_variant,
    show_gene_names=args.show_gene_names
)

# -> Seaborn clustermap:

def plot_clustermap(gtex_df, tcga_df, title, output_path=None, top_n=500, show_gene_names=False):
    gtex_df = gtex_df.set_index('Description')
    tcga_df = tcga_df.set_index('sample')

    combined = pd.concat([gtex_df, tcga_df], axis=1)  # concatenate
    combined_log = combined.apply(lambda x: np.log1p(x))  # log

    variances = combined_log.var(axis=1)  # variance for each gene
    top_variant_genes = variances.nlargest(min(top_n, len(variances))).index  # pick n genes with the highest variance, without selecting more than available genes
    most_variant = combined_log.loc[top_variant_genes]  # filter the most variant

    # colors for gtex and tcga columns
    gtex_cols = [plot_colors['GTEX']] * gtex_df.shape[1]
    tcga_cols = [plot_colors['TCGA']] * tcga_df.shape[1]
    col_colors = gtex_cols + tcga_cols

    sns.set(font_scale=1.0)
    cmap = sns.clustermap(
        most_variant,
        cmap="viridis",
        metric="euclidean",
        method="average",  # more balanced clusters that are not too sensitive to outliers
        figsize=(20, 20),
        xticklabels=False,  # remove sample names
        yticklabels=show_gene_names,  # show gene names if specified
        col_colors=col_colors  # assign colors for gtex and tcga columns
    )

    # legend for column colors
    legend_labels = ['GTEx', 'TCGA']
    legend_colors = [plot_colors['GTEX'], plot_colors['TCGA']]
    for label, color in zip(legend_labels, legend_colors):
        plt.plot([], [], marker="o", label=label, color=color, ls="", markersize=18)

    cmap.fig.legend(
        title='Dataset',
        title_fontsize=28,
        fontsize=24,
        bbox_to_anchor=(1.0, 0.90),  # first number: 0 left, 1 right, second number: 0 bottom, 1 top
        loc='upper right',
        markerscale=1.5,
    )

    # title
    cmap.fig.suptitle(
        f"{title} \n(Top {top_n} variant genes, log-transformed)",
        fontsize=35,
        y=0.95,  # position at the top
        x=0.5,  # center horizontally
        ha='center'  # horizontal alignment
    )

    if output_path:
        plt.tight_layout(rect=[0, 0, 1, 0.95])  # provide space above the plot
        cmap.savefig(output_path, dpi=600, format="pdf", bbox_inches="tight")
    plt.close()


def plot_clustermap_combined(gtex_lung_df, gtex_breast_df, tcga_lung_df, tcga_breast_df, title, output_path=None, top_n=500, show_gene_names=False):
    '''
    print("Any NaN in gtex_lung_df:", gtex_lung_df.isna().any().any())
    print("Any NaN in gtex_breast_df:", gtex_breast_df.isna().any().any())
    print("Any NaN in tcga_lung_df:", tcga_lung_df.isna().any().any())
    print("Any NaN in tcga_breast_df:", tcga_breast_df.isna().any().any())
    print("nrow gtex_lung_df:", gtex_lung_df.shape[0])
    print("nrow gtex_breast_df:", gtex_breast_df.shape[0])
    print("nrow tcga_lung_df:", tcga_lung_df.shape[0])
    print("nrow tcga_breast_df:", tcga_breast_df.shape[0])
    '''

    # important: set row index to 'Description' or 'sample' for the case that dfs filtered for common genes are used
    gtex_lung_df = gtex_lung_df.set_index('Description')
    tcga_lung_df = tcga_lung_df.set_index('sample')
    gtex_breast_df = gtex_breast_df.set_index('Description')
    tcga_breast_df = tcga_breast_df.set_index('sample')

    # important: modify column names to avoid duplicates
    gtex_lung_df.columns = [f"gtex_lung_{col}" for col in gtex_lung_df.columns]
    tcga_lung_df.columns = [f"tcga_lung_{col}" for col in tcga_lung_df.columns]
    gtex_breast_df.columns = [f"gtex_breast_{col}" for col in gtex_breast_df.columns]
    tcga_breast_df.columns = [f"tcga_breast_{col}" for col in tcga_breast_df.columns]

    # important: take only the numeric columns
    gtex_lung_df = gtex_lung_df.select_dtypes(include=[np.number])
    gtex_breast_df = gtex_breast_df.select_dtypes(include=[np.number])
    tcga_lung_df = tcga_lung_df.select_dtypes(include=[np.number])
    tcga_breast_df = tcga_breast_df.select_dtypes(include=[np.number])

    # concatenate the dataframes column-wise
    combined_data = pd.concat(
        [gtex_lung_df, gtex_breast_df, tcga_lung_df, tcga_breast_df],
        axis=1
    )

    # metadata for the columns
    column_metadata = pd.DataFrame({
        'Source': ['GTEx'] * gtex_lung_df.shape[1] + ['GTEx'] * gtex_breast_df.shape[1] +
                  ['TCGA'] * tcga_lung_df.shape[1] + ['TCGA'] * tcga_breast_df.shape[1],
        'Tissue': ['Lung'] * gtex_lung_df.shape[1] + ['Breast'] * gtex_breast_df.shape[1] +
                  ['Lung'] * tcga_lung_df.shape[1] + ['Breast'] * tcga_breast_df.shape[1]
    }, index=combined_data.columns)

    '''
    if (combined_data.select_dtypes(include=[np.number]) <0).any().any():
        print("combined_data contains negative values")
    else: print("combined_data contains no negative values")
    print("combined_data contains NaN?: ", combined_data.isna().any().any())
    '''

    numeric_data = combined_data.select_dtypes(include=[np.number])  # filter only numeric columns

    '''
    if (numeric_data < 0).any().any():
        print("numeric_data contains negative values (*)")
    else: print("numeric_data contains no negative values (*)")
    print("numeric_data contains NaN? (*): ", numeric_data.isna().any().any())
    '''

    # before logarithm, replace small negative values with 0
    threshold = 0
    numeric_data[numeric_data < threshold] = numeric_data[numeric_data < threshold].clip(lower=0)

    '''
    if (numeric_data <0).any().any():
        print("numeric_data STILL contains negative values (**)")
    else: print("numeric_data NO LONGER contains negative values (**)")
    print("numeric_data contains zeros? (*): ", (numeric_data == 0).any().any())
    '''

    numeric_data += 1e-8  # add small constant to avoid log(0) in the following

    '''
    if (numeric_data < 0).any().any():
        print("numeric_data STILL contains negative values (***)")
    else: print("numeric_data NO LONGER contains negative values (***)")
    print("numeric_data contains zeros after pseudocount? (**): ", (numeric_data == 0).any().any())
    print("numeric_data contains NaN? (**): ", numeric_data.isna().any().any())
    '''

    scaled_data = numeric_data.apply(lambda x: np.log1p(x)) # log

    '''
    print("scaled_data contains Negatives?: ", (scaled_data <0).any().any())
    print("scaled_data contains zeros?: ", (scaled_data == 0).any().any())
    print("scaled_data contains NaN?: ", scaled_data.isna().any().any())
    '''

    variances = scaled_data.var(axis=1)  # variance for each gene
    top_variant_genes = variances.nlargest(min(top_n, len(variances))).index  # pick n genes with the highest variance, without selecting more than available genes
    most_variant = scaled_data.loc[top_variant_genes]  # filter the most variant

    '''
    print("Any NaN values in most_variant:", most_variant.isna().any().any())
    print("Any Inf values in most_variant:", np.isinf(most_variant).any().any())
    print("+++++")
    print("NaN values in variances: ", variances.isna().sum())
    print("NaN values in most_variant: ", most_variant.isna().sum())
    if np.isinf(most_variant.values).any():
        print("most_variant contains infinity values.")
    else: print("most_variant has no inf values")
    if most_variant.index.duplicated().any():
        print("most_variant contains duplicate rows.")
    else: print("most_variant has no duplicate rows")
    zero_variance_rows = most_variant.var(axis=1) == 0
    if zero_variance_rows.any():
        print(f"{zero_variance_rows.sum()} rows of most_variant have zero variance.")
    print("Any NaN values in most_variant:", most_variant.isna().any().any())
    print("Any inf values in most_variant:", np.isinf(most_variant.values).any())
    '''

    duplicates = combined_data.columns[combined_data.columns.duplicated()]
    if len(duplicates) > 0:
        # add suffix to duplicate column names to make them unique
        combined_data.columns = [f"{col}_{i}" if combined_data.columns.tolist().count(col) > 1 else col
                                 for i, col in enumerate(combined_data.columns)]

        # ensure metadata still corresponds to the correct columns
        column_metadata = pd.DataFrame({
            'Source': ['GTEx'] * gtex_lung_df.shape[1] + ['GTEx'] * gtex_breast_df.shape[1] +
                      ['TCGA'] * tcga_lung_df.shape[1] + ['TCGA'] * tcga_breast_df.shape[1],
            'Tissue': ['Lung'] * gtex_lung_df.shape[1] + ['Breast'] * gtex_breast_df.shape[1] +
                      ['Lung'] * tcga_lung_df.shape[1] + ['Breast'] * tcga_breast_df.shape[1]
        }, index=combined_data.columns)

    # map metadata to colors
    source_palette = {'GTEx': plot_colors['GTEX'], 'TCGA': plot_colors['TCGA']}
    tissue_palette = {'Lung': plot_colors['lung'], 'Breast': plot_colors['breast']}

    col_colors = pd.DataFrame({
        'Source': column_metadata['Source'].map(source_palette),
        'Tissue': column_metadata['Tissue'].map(tissue_palette)
    }, index=column_metadata.index)

    sns.set(font_scale=1.0)
    cmap = sns.clustermap(
        most_variant,
        cmap="viridis",
        metric="euclidean",
        method="average",  # more balanced clusters that are not too sensitive to outliers
        figsize=(20, 20),
        xticklabels=False,  # remove sample names
        yticklabels=show_gene_names,  # show gene names if specified
        col_colors=col_colors
    )

    # legend for the metadata (labels)
    for label, color in source_palette.items():
        plt.plot([], [], marker="o", label=label, color=color, ls="", markersize=18)

    for label, color in tissue_palette.items():
        plt.plot([], [], marker="o", label=label, color=color, ls="", markersize=18)

    cmap.fig.legend(
        title='Metadata',
        title_fontsize=28,
        fontsize=24,
        bbox_to_anchor=(1.0, 0.90),  # first number: 0 left, 1 right, second number: 0 bottom, 1 top
        loc='upper right',
        markerscale=1.5,
    )

    # title
    cmap.fig.suptitle(
        f"{title} \n(Top {top_n} variant genes, log-transformed)",
        fontsize=35,
        y=0.95,  # position at the top
        x=0.5,  # center horizontally
        ha='center'  # horizontal alignment
    )

    if output_path:
        plt.tight_layout(rect=[0, 0, 1, 0.95])  # provide space above the plot
        cmap.savefig(output_path, dpi=600, format="pdf", bbox_inches="tight")
    plt.close()


# clustermap for lung data only:
plot_clustermap(
    gtex_lung_intersection,
    tcga_lung_intersection,
    "Clustermap for Lung (GTEx and TCGA Intersections)",
    output_path=os.path.join(output_dir, "clustermap_lung.pdf"),
    top_n=args.top_n_variant,
    show_gene_names=args.show_gene_names
)

# clustermap for breast data only:
plot_clustermap(
    gtex_breast_intersection,
    tcga_breast_intersection,
    "Clustermap for Breast (GTEx and TCGA Intersections)",
    output_path=os.path.join(output_dir, "clustermap_breast.pdf"),
    top_n=args.top_n_variant,
    show_gene_names=args.show_gene_names
)

# clustermap for lung and breast data together: - old version without lung/breast labeling
plot_clustermap(
    gtex_combined,
    tcga_combined,
    "Combined Clustermap (Lung and Breast Intersections)",
    output_path=os.path.join(output_dir, "clustermap_combined.pdf"),
    top_n=args.top_n_variant,
    show_gene_names=args.show_gene_names
)

# -> clustermap for lung and breast data together with lung/breast labeling
plot_clustermap_combined(
    gtex_lung_df=gtex_lung_common,
    gtex_breast_df=gtex_breast_common,
    tcga_lung_df=tcga_lung_common,
    tcga_breast_df=tcga_breast_common,
    title="Combined Clustermap (Lung and Breast Intersections)",
    output_path=os.path.join(output_dir, "clustermap_combined_new.pdf"),
    top_n=args.top_n_variant,
    show_gene_names=args.show_gene_names
)

'''
print("--- gtex_lung_intersection:")
print("Any NaNs at all?: ", gtex_lung_intersection.isna().any().any())
print("Number of NaNs overall: ", gtex_lung_intersection.isna().sum().sum())
print("Negatives?: ", (gtex_lung_intersection.select_dtypes(include=['number']) < 0).any().any())
print(gtex_lung_intersection.head())
print(gtex_lung_intersection.select_dtypes(include=['number']).head())
print("--- gtex_breast_intersection:")
print("Any NaNs at all?: ", gtex_breast_intersection.isna().any().any())
print("Number of NaNs overall: ", gtex_breast_intersection.isna().sum().sum())
print("Negatives?: ", (gtex_breast_intersection.select_dtypes(include=['number']) < 0).any().any())
print("--- tcga_lung_intersection:")
print("Any NaNs at all?: ", tcga_lung_intersection.isna().any().any())
print("Number of NaNs overall: ", tcga_lung_intersection.isna().sum().sum())
print("Negatives?: ", (tcga_lung_intersection.select_dtypes(include=['number']) < 0).any().any())
print("--- tcga_breast_intersection:")
print("Any NaNs at all?: ", tcga_breast_intersection.isna().any().any())
print("Number of NaNs overall: ", tcga_breast_intersection.isna().sum().sum())
print("Negatives?: ", (tcga_breast_intersection.select_dtypes(include=['number']) < 0).any().any())
'''

# Reduce number of columns:

least_cols = min(len(gtex_lung.columns), len(tcga_lung.columns), len(gtex_breast.columns), len(tcga_breast.columns))
print("Smallest number of columns in one of the datasets: ", least_cols)

def reduce_columns(data, n_columns_to_sample, metadata_column_names):
    n_columns_to_sample = n_columns_to_sample - len(metadata_column_names)  # subtract number of meta columns to get the number of numeric columns
    numeric_columns = data.select_dtypes(include=[float, int]).columns
    columns_to_sample = data[numeric_columns].sample(
        n=n_columns_to_sample,
        axis=1,
        random_state=42  # for reproducibility
    ).columns

    sampled_df = data[metadata_column_names + list(columns_to_sample)]
    sampled_df = sampled_df.reset_index(drop=True)
    return sampled_df

# reduce all dfs to the smallest number of columns, attention: this number includes 'Description'/'sample' column
gtex_lung_intersection_reduced = reduce_columns(gtex_lung_intersection, least_cols, ['Description'])
tcga_lung_intersection_reduced = reduce_columns(tcga_lung_intersection, least_cols, ['sample'])
gtex_breast_intersection_reduced = reduce_columns(gtex_breast_intersection, least_cols, ['Description'])
tcga_breast_intersection_reduced = reduce_columns(tcga_breast_intersection, least_cols, ['sample'])

# for combined reduced: filter for common genes for both tissues (as above)
gtex_lung_reduced_common = gtex_lung_intersection_reduced[gtex_lung_intersection_reduced['Description'].isin(common_genes_all)]
tcga_lung_reduced_common = tcga_lung_intersection_reduced[tcga_lung_intersection_reduced['sample'].isin(common_genes_all)]
gtex_breast_reduced_common = gtex_breast_intersection_reduced[gtex_breast_intersection_reduced['Description'].isin(common_genes_all)]
tcga_breast_reduced_common = tcga_breast_intersection_reduced[tcga_breast_intersection_reduced['sample'].isin(common_genes_all)]

# clustermap for lung and breast data combined with reduced columns:
plot_clustermap_combined(gtex_lung_df=gtex_lung_reduced_common, gtex_breast_df=gtex_breast_reduced_common,
                         tcga_lung_df=tcga_lung_reduced_common, tcga_breast_df=tcga_breast_reduced_common,
                         title="Clustermap for Lung and Breast (GTEx and TCGA Intersections) \n- reduced number of samples",
                         output_path=os.path.join(output_dir, "clustermap_combined_reduced.pdf"),
                         top_n=args.top_n_variant,
                         show_gene_names=args.show_gene_names
                         )

# clustermap for lung data only, with reduced columns:
plot_clustermap(
    gtex_lung_intersection_reduced,
    tcga_lung_intersection_reduced,
    "Clustermap for Lung (GTEx and TCGA Intersections) \n- reduced number of samples",
    output_path=os.path.join(output_dir, "clustermap_lung_reduced.pdf"),
    top_n=args.top_n_variant,
    show_gene_names=args.show_gene_names
)

# clustermap for breast data only, with reduced columns:
plot_clustermap(
    gtex_breast_intersection_reduced,
    tcga_breast_intersection_reduced,
    "Clustermap for Breast (GTEx and TCGA Intersections) \n- reduced number of samples",
    output_path=os.path.join(output_dir, "clustermap_breast_reduced.pdf"),
    top_n=args.top_n_variant,
    show_gene_names=args.show_gene_names
)
