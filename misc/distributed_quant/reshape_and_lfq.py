# Python script to reshape input precursor data and run directlfq quantification

import argparse
import pandas as pd
import os
import gc
import sys
import yaml
import directlfq.utils
import directlfq.lfq_manager

# parse input parameters
parser = argparse.ArgumentParser(
    prog="DistributedDirectLFQ",
    description="Reshape input precursor data and run DirectLFQ quantification"
)
parser.add_argument("--config", type=str, required=True)
args = parser.parse_args()

# read the config .yaml file into a dictionary
if not os.path.exists(args.config):
    print(f"Could not find {args.config}.")
    sys.exit(1)

with open(args.config) as file:
    config = yaml.safe_load(file)

# Parse input parameters explicitly
protein_list = config["proteins"]
parquet_file_path = config["precursor_file"]
results_dir = config["output_dir"]

# Functions to handle the input data
_, config_dict, _ = directlfq.utils.get_input_type_and_config_dict(parquet_file_path)

# Get the proteins from the .parquet file
filters = [("Protein.Group", "in", protein_list)]
df_sub = pd.read_parquet(parquet_file_path, filters=filters)

# Adapt and reshape in memory. In the normal directlfq.lfq_manager, this happens as the 
# file is loaded. Here, we bypass the loading step and directly convert the table into
# the format corresponding to the 'aq_reshaped' data.
df_sub_adapted = directlfq.utils.adapt_subtable(df_sub, config_dict)
df_sub_reshaped = directlfq.utils.reshape_input_df(df_sub_adapted, config_dict)

# Importantly, force all runs to be present in the output: since individual protein outputs
# will be merged on runs, it is crucial that all runs are retained in the chunk outputs
df_metadata = pd.read_parquet(parquet_file_path, columns=[config_dict.get("sample_ID")])
all_sample_ids = df_metadata[config_dict.get("sample_ID")].unique().tolist()
del df_metadata
df_sub_reshaped = df_sub_reshaped.reindex(columns=['protein', 'ion'] + all_sample_ids, fill_value=0)

# Run the lfq manager with dataframe input
directlfq.lfq_manager.run_lfq(
    input_df=df_sub_reshaped, 
    outfile_path=results_dir, # the current chunk output directory
    deactivate_normalization=True,
    compile_normalized_ion_table=False,
)

del df_sub, df_sub_adapted, df_sub_reshaped
gc.collect()

print(f"DirectLFQ finished in {results_dir}.")


