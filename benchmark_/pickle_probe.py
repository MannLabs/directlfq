import time, pickle
import numpy as np
import directlfq.config as config
import directlfq.utils as lfqutils
import directlfq.normalization as lfqnorm
import directlfq.protein_intensity_estimation as lfqprot
INPUT="/home/magnus/speedup_directlfq/directlfq/input_data/HYE1-6_repeat/lfq.dia.features.csv.aq_reformat.tsv"
config.set_global_protein_and_ion_id("protein","ion"); config.set_log_processed_proteins(False)
config.set_compile_normalized_ion_table(True); config.check_wether_to_copy_numpy_arrays_derived_from_pandas()
inp=lfqutils.add_mq_protein_group_ids_if_applicable_and_obtain_annotated_file(INPUT,None,None,[])
df=lfqutils.import_data(input_file=inp,input_type_to_use=None,filter_dict=None)
df=lfqutils.sort_input_df_by_protein_and_quant_id(df); df=lfqutils.remove_potential_quant_id_duplicates(df)
df=lfqutils.index_and_log_transform_input_df(df); df=lfqutils.remove_allnan_rows_input_df(df)
normed=lfqnorm.NormalizationManagerSamplesOnSelectedProteins(df,num_samples_quadratic=50,selected_proteins_file=None).complete_dataframe
dfs=lfqprot.get_normed_dfs(normed)
print("n sub-dataframes:", len(dfs))
# pickle the list of per-protein DataFrames (what Pool ships to workers)
t=time.perf_counter(); blob=pickle.dumps(dfs, protocol=pickle.HIGHEST_PROTOCOL); t_df=time.perf_counter()-t
print(f"DataFrames: pickle {t_df:.2f}s  size {len(blob)/1e6:.1f} MB")
t=time.perf_counter(); _=pickle.loads(blob); t_dfl=time.perf_counter()-t
print(f"DataFrames: unpickle {t_dfl:.2f}s")
# equivalent as raw numpy arrays (values only)
arrs=[d.to_numpy() for d in dfs]
t=time.perf_counter(); blob2=pickle.dumps(arrs, protocol=pickle.HIGHEST_PROTOCOL); t_np=time.perf_counter()-t
print(f"numpy arrays: pickle {t_np:.2f}s  size {len(blob2)/1e6:.1f} MB")
