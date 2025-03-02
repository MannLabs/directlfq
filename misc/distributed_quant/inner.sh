#!/usr/bin/env bash

# Script to run distributed reshaping and quantification via DirectLFQ

#SBATCH --job-name=DLFQ
#SBATCH --time=2-00:00:00
#SBATCH --output=./logs/%A_%a-%x-slurm.out

# the script directory must be passed by the outer script
if [[ -z "${script_dir:-}" ]]; then
    echo "Error: script_dir is not set. Ensure it is passed via sbatch --export."
    exit 1
fi

# output_dir is passed by the outer script
if [[ -z "${output_dir:-}" ]]; then
    echo "Error: output_dir is not set. Ensure it is passed via sbatch --export."
    exit 1
fi

# Navigate to chunk directory
slurm_index=${SLURM_ARRAY_TASK_ID}
chunk_directory="${output_dir}/chunk_${slurm_index}/"

if ! cd "${chunk_directory}"; then
    echo "Failed to navigate to chunk directory: ${chunk_directory}"
    exit 1
fi

# Set the config file name
config_filename="dist_lfq_config.yaml"

if [[ ! -f "${config_filename}" ]]; then
    echo "Config file not found: ${config_filename}"
    exit 1
fi

# run directlfq wrapper script, which is in the same directory as this script
python "${script_dir}/reshape_and_lfq.py" --config "${config_filename}"
# python ./reshape_and_lfq.py --config "${config_filename}"

echo "DirectLFQ completed successfully"

