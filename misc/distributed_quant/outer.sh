#!/usr/bin/env bash

# Main script, set source and destination folders and creates chunk
# folders for each batch of proteins to search: 
# 1. Read in I/O parameters, chunk size (C) and number of nodes to occupy in parallel (n)
# 2. Read all protein groups from the parquet file (N)
# 3. Define the sbatch array string: N // C + 1 if N % C else N // C --> this sets the number of chunks
# 4. Create a chunk directory in the output folder for each chunk (i)
# 5. Subset the protein list into as many lists as there are chunks
# 5. Write the configuration yaml into each chunk directory: 
#     - the respective chunk's protein list
#     - the source parquet directory
#     - the destination directory
# 6. Submit the job array to the cluster

#SBATCH --job-name=dist_LFQ
#SBATCH --time=2-00:00:00
#SBATCH --output=./logs/%j-%x-slurm.out

# Set behavior when errors are encountered
set -u -x

# default parameters
nnodes=1
ntasks_per_node=1 # only change if you know that your HPCL supports more than 1 task per node
cpus=32
mem='250G'

# create logs directory for slurm jobs if it does not exist
mkdir -p ./logs

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --precursor_file) precursor_file="$2"; shift ;;
        --output_dir) output_dir="$2"; shift ;;
        --nnodes) nnodes="$2"; shift ;;
        --chunk_size) chunk_size="$2"; shift ;;
        *) echo "Unknown parameter $1"; exit 1 ;;
    esac
    shift
done

# Verify that the precursor file exists
if [[ ! -f "${precursor_file}" ]]; then
    echo "Precursor file not found: ${precursor_file}"
    exit 1
fi

# Verify that the output directory exists
mkdir -p "${output_dir}"

# Run the python script to handle chunk creation and array generation and get the number of tasks
num_tasks=$(python ./prepare_quant_chunks.py \
    --precursor_file "${precursor_file}" \
    --output_dir "${output_dir}" \
    --chunk_size "${chunk_size}" | tail -n 1)

# Check on num_tasks
if ! [[ "$num_tasks" =~ ^[0-9]+$ ]]; then
    echo "Failed to get the number of tasks"
    exit 1
fi

# Generate the sbatch array
slurm_array="0-$(($num_tasks-1))%${nnodes}"
echo "Submitting quantification jobs: ${slurm_array}"

# Submit the job array (Slurm passes the array index as the first argument to the script)
sbatch --array=${slurm_array} \
    --wait \
    --nodes=1 \
    --ntasks-per-node=${ntasks_per_node} \
    --cpus-per-task=${cpus} \
    --mem=${mem} \
    --export=ALL,output_dir=${output_dir} ./reshape_and_lfq.sh

# Check if the job submission was successful
if [[ $? -ne 0 ]]; then
    echo "Error: SLURM job submission failed."
    exit 1
fi

echo "All quantification chunks have been run"

# Run python script to traverse all chunk directories and concatenate the protein tables
python ./concatenate_protein_tables.py \
    --input_dir "${output_dir}"

echo "All quantification chunks have been concatenated, output table written to ${output_dir}/protein_table.tsv"