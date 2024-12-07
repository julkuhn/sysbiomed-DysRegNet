import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--input_meta', type=str, required=True)
parser.add_argument('--input_tcga_data', type=str, required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()

meta_input = args.input_meta
tpm_tcga_input = args.input_tcga_data
output_file = args.output

#read in
meta = pd.read_csv(meta_input, sep=",")
tcga_data = pd.read_csv(tpm_tcga_input, sep=",")

# Check how the control samples are called inside the sample_type column of the meta file
# print(meta['sample_type'].unique())  # -> "Solid Tissue Normal" represents controls
print("Assumes the controls are named 'Solid Tissue Normal' inside the sample_type column of the meta file.")

# filter the controls
controls = meta[meta['sample_type'] == 'Solid Tissue Normal']

# names of the control samples
control_samples = controls['sample'].tolist()

# filter the TCGA data for these control samples
tcga_data_filtered = tcga_data[control_samples]

# re-insert the "sample" column of the TCGA data into the filtered tcga
tcga_data_filtered.insert(0, 'sample', tcga_data['sample'])

print("Number of TCGA columns (samples) before filtering: ", tcga_data.shape[1] -1 )  # minus the sample column
print("Number of TCGA columns (samples) filtered for controls: ", tcga_data_filtered.shape[1] -1)

# write filtered tcga data to file
tcga_data_filtered.to_csv(args.output, index=False)

print(f"Output file saved to: {output_file}")
