import time, sys
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
spec = list(
    lfqprot.get_input_specification_tuplelist_idx__df__num_samples_quadratic__min_nonan(
        normed, 10, 1
    )
)
# process just first 50 proteins single-core: time dominated by one-time JIT compile
t = time.perf_counter()
lfqprot.get_list_with_sequential_processing(iter(spec[:50]))
t1 = time.perf_counter() - t
# process next 50 (warm): pure compute
t = time.perf_counter()
lfqprot.get_list_with_sequential_processing(iter(spec[50:100]))
t2 = time.perf_counter() - t
print(
    f"first50_cold={t1:.2f}s  next50_warm={t2:.2f}s  -> one-time JIT+overhead ~= {t1 - t2:.2f}s per process"
)
