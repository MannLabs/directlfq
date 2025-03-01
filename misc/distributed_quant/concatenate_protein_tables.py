import argparse
import pandas as pd
import os
import gc
import sys

# print to terminal
sys.stdout.flush()

# parse input parameters
parser = argparse.ArgumentParser(
    prog="DistributedDirectLFQ",
    description="Reshape input precursor data and run DirectLFQ quantification"
)
parser.add_argument("--input_dir", type=str, required=True)
parser.add_argument("--output_dir", type=str, required=False)
args = parser.parse_args()

if not args.output_dir:
    args.output_dir = args.input_dir

# Iterate over all files in the input directory and map the chunk_X folders
if not os.path.exists(args.input_dir):
    print(f"Could not find {args.input_dir}.")
    sys.exit(1) 

all_items = os.listdir(args.input_dir)
all_directories = [f for f in all_items if os.path.isdir(os.path.join(args.input_dir, f))]
chunk_directories = [f for f in all_directories if f.startswith("chunk_")]

# Sort the directories by their chunk number
chunk_directories = sorted(chunk_directories, key=lambda x: int(x.split("_")[1]))

# Iterate over all chunk directories and concatenate the protein tables
protein_dfs = []
for chunk_directory in chunk_directories:
    files = os.listdir(os.path.join(args.input_dir, chunk_directory))
    protein_table_path = [f for f in files if f.endswith("protein_intensities.tsv")]
    if len(protein_table_path) == 0:
        print(f"No protein table found in {chunk_directory}, skipping.")
        continue
    elif len(protein_table_path) > 1:
        print(f"Error: Multiple protein tables found in {chunk_directory}, skipping.")
        sys.exit(1)
    protein_table_path = protein_table_path[0]
    protein_df = pd.read_csv(os.path.join(args.input_dir, chunk_directory, protein_table_path), sep="\t")
    if len(protein_df) == 0:
        print(f"Empty protein table found in {chunk_directory}, skipping.")
        continue
    protein_dfs.append(protein_df)

# Concatenate all protein tables
if len(protein_dfs) == 0:
    print("No protein tables found, exiting.")
    sys.exit(1)

# Join on column indices
protein_df = pd.concat(protein_dfs, axis=0, ignore_index=True)

# Check if there are any duplicated protein names
if any(protein_df["protein"].duplicated()):
    print("Error: Duplicated protein names found in the concatenated protein tables.")
    sys.exit(1)

# Save the concatenated protein table
output_path = os.path.join(args.output_dir, "concatenated_protein_intensities.tsv")
protein_df.to_csv(output_path, sep="\t", index=False)

