import cProfile, pstats, io, time
import directlfq.config as config, directlfq.utils as lfqutils
import directlfq.normalization as lfqnorm
INPUT="/home/magnus/speedup_directlfq/directlfq/input_data/HYE1-6_repeat/lfq.dia.features.csv.aq_reformat.tsv"
config.set_global_protein_and_ion_id("protein","ion"); config.set_log_processed_proteins(False)
config.set_compile_normalized_ion_table(True); config.check_wether_to_copy_numpy_arrays_derived_from_pandas()
inp=lfqutils.add_mq_protein_group_ids_if_applicable_and_obtain_annotated_file(INPUT,None,None,[])
df=lfqutils.import_data(input_file=inp,input_type_to_use=None,filter_dict=None)
df=lfqutils.sort_input_df_by_protein_and_quant_id(df); df=lfqutils.remove_potential_quant_id_duplicates(df)
df=lfqutils.index_and_log_transform_input_df(df); df=lfqutils.remove_allnan_rows_input_df(df)
# warm + time
t=time.perf_counter(); _=lfqnorm.NormalizationManagerSamplesOnSelectedProteins(df.copy(), num_samples_quadratic=50).complete_dataframe; print("norm time:", round(time.perf_counter()-t,3))
pr=cProfile.Profile(); pr.enable()
_=lfqnorm.NormalizationManagerSamplesOnSelectedProteins(df.copy(), num_samples_quadratic=50).complete_dataframe
pr.disable()
s=io.StringIO(); pstats.Stats(pr,stream=s).sort_stats('cumtime').print_stats(20); print(s.getvalue())
