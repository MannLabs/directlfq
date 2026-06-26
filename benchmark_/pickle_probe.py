import time, pickle
import directlfq.config as config
import directlfq.utils as lfqutils
import directlfq.normalization as lfqnorm
import directlfq.protein_intensity_estimation as lfqprot

INPUT = "/home/magnus/speedup_directlfq/directlfq/input_data/HYE1-6_repeat/lfq.dia.features.csv.aq_reformat.tsv"
config.set_global_protein_and_ion_id("protein", "ion")
config.set_log_processed_proteins(False)
config.set_compile_normalized_ion_table(True)
config.check_wether_to_copy_numpy_arrays_derived_from_pandas()
inp = lfqutils.add_mq_protein_group_ids_if_applicable_and_obtain_annotated_file(
    INPUT, None, None, []
)
df = lfqutils.import_data(input_file=inp, input_type_to_use=None, filter_dict=None)
df = lfqutils.sort_input_df_by_protein_and_quant_id(df)
df = lfqutils.remove_potential_quant_id_duplicates(df)
df = lfqutils.index_and_log_transform_input_df(df)
df = lfqutils.remove_allnan_rows_input_df(df)
normed = lfqnorm.NormalizationManagerSamplesOnSelectedProteins(
    df, num_samples_quadratic=50, selected_proteins_file=None
).complete_dataframe
workitems = list(lfqprot.get_protein_workitems(normed, 10, 1))
print("n work items:", len(workitems))
# pickle the list of per-protein numpy work items (what Pool ships to workers)
t = time.perf_counter()
blob = pickle.dumps(workitems, protocol=pickle.HIGHEST_PROTOCOL)
t_wi = time.perf_counter() - t
print(f"work items: pickle {t_wi:.2f}s  size {len(blob) / 1e6:.1f} MB")
t = time.perf_counter()
_ = pickle.loads(blob)
t_wil = time.perf_counter() - t
print(f"work items: unpickle {t_wil:.2f}s")
