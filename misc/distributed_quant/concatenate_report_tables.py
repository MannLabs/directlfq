import argparse
import pandas as pd
import os
import gc
import sys

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
    protein_table = [f for f in files if f.endswith("protein_intensities.tsv")]
    if len(protein_table) == 0:
        print(f"No protein table found in {chunk_directory}, skipping.")
        continue
    elif len(protein_table) > 1:
        print(f"Multiple protein tables found in {chunk_directory}, skipping.")
        continue
    protein_table = protein_table[0]