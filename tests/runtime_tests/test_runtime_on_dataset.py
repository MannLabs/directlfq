import sys
import directlfq.benchmarking as lfq_benchmark
import pandas as pd


args = sys.argv[1:]
template_file = args[0]
number_of_samples = int(args[1])

lfq_timer = lfq_benchmark.LFQTimer(template_file, list_of_samplenumbers_to_check=[number_of_samples], name="check_dataset")

timedrun = lfq_timer.timed_lfq_runs[0]


df_run = pd.Series({'time_normalization' : timedrun.runtime_info.runtime_samplenorm, 
                    'time_protein_norm' : timedrun.runtime_info.runtime_protein_norm,
                    'time_overall': timedrun.runtime_info.overall_runtime})

template_name = template_file.split("/")[-1]
outfile_name = f"runtimes_{template_name}_{number_of_samples}.tsv"
df_run.to_csv(outfile_name, sep = "\t", header=None)
print(f"Saved runtimes to {outfile_name}")







