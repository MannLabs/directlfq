# Python script to prepare distributed DirectLFQ chunks

import argparse
import pandas as pd
import os
import sys
import yaml

# Add keys to config if they don't exist
def safe_add_key(config, parent_key, key, value):
    """Safely add keys and values to config file dictionary"""
    if parent_key is None:
        config[key] = value
    else:
        if parent_key not in config:
            config[parent_key] = {}
        config[parent_key][key] = value

# parse input parameters
parser = argparse.ArgumentParser(
    prog="DistributedDirectLFQ",
    description="Generate chunk folders to run distributed DirectLFQ quantifications"
)
parser.add_argument("--precursor_file", type=str, required=True)
parser.add_argument("--output_dir", type=str, required=True)
parser.add_argument("--chunk_size", type=int, required=True)
args = parser.parse_args()

# Check that the precursor file exists
if not os.path.exists(args.precursor_file):
    print(f"Could not find {args.precursor_file}.")
    sys.exit(1)

# Check that the output directory exists, otherwise create it
if not os.path.exists(args.output_dir):
    try:
        os.makedirs(args.output_dir)
    except Exception:
        print(f"Could not create {args.output_dir}.")    
        sys.exit(1)

# Get the proteins from the .parquet file
protein_df = pd.read_parquet(args.precursor_file, columns = ['Protein.Group'])
unique_proteins = protein_df['Protein.Group'].unique()
num_proteins = len(unique_proteins)
print(f"The precursor file contains {num_proteins} unique proteins.")

if num_proteins == 0:
    print("No proteins found in the precursor file, nothing to do.")
    sys.exit(1)

# Determine the number of chunks: shift num_proteins up so that even the smallest possible remainder (i.e. 1)
# causes an extra chunk to be created. But subtract the smallest possible remainder to account for the case 
# when num_proteins is in fact cleanly divisible by the chunk size.
n_chunks = (num_proteins + args.chunk_size - 1) // args.chunk_size

# Iterate the number of chunks and generate the respective folders in the output
for i in range(n_chunks):

    # get the corresponding slice of the protein list
    current_proteins = list(unique_proteins[i*args.chunk_size:(i+1)*args.chunk_size])
    
    # write the output config
    out_config = {
        'proteins' : current_proteins,
        'precursor_file' : args.precursor_file,
        'output_dir' : os.path.join(args.output_dir,f"chunk_{i}"),
    }

    # create the chunk folder
    os.makedirs(os.path.join(args.output_dir,f"chunk_{i}"), exist_ok=True)

    # overwrite the existing config file
    with open(os.path.join(args.output_dir,f"chunk_{i}/dist_lfq_config.yaml"), 'w') as file:
        yaml.dump(out_config, file)

# Return the number of chunks, i.e. the number of tasks for the job scheduler
print(n_chunks, flush=True)