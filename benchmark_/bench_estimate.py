"""Break the 'estimate' stage into its serial vs parallel sub-steps.

estimate_protein_intensities does, in order:
  1. build per-protein work items (get_protein_workitems) -> SERIAL
  2. map calculate_peptide_and_protein_intensities       -> PARALLEL (pool) / serial
  3. assemble protein dataframe                           -> SERIAL
  4. compile ion dataframe                                -> SERIAL

Usage: python bench_estimate.py <num_cores>
"""

import sys
import time

import directlfq.config as config
import directlfq.utils as lfqutils
import directlfq.normalization as lfqnorm
import directlfq.protein_intensity_estimation as lfqprot

INPUT = (
    "/home/magnus/speedup_directlfq/directlfq/input_data/"
    "HYE1-6_repeat/lfq.dia.features.csv.aq_reformat.tsv"
)


def main() -> None:
    num_cores = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    config.set_global_protein_and_ion_id(protein_id="protein", quant_id="ion")
    config.set_log_processed_proteins(log_processed_proteins=False)
    config.set_compile_normalized_ion_table(compile_normalized_ion_table=True)
    config.check_wether_to_copy_numpy_arrays_derived_from_pandas()

    inp = lfqutils.add_mq_protein_group_ids_if_applicable_and_obtain_annotated_file(
        INPUT, None, None, []
    )
    df = lfqutils.import_data(input_file=inp, input_type_to_use=None, filter_dict=None)
    df = lfqutils.sort_input_df_by_protein_and_quant_id(df)
    df = lfqutils.remove_potential_quant_id_duplicates(df)
    df = lfqutils.index_and_log_transform_input_df(df)
    df = lfqutils.remove_allnan_rows_input_df(df)
    normed_df = lfqnorm.NormalizationManagerSamplesOnSelectedProteins(
        df, num_samples_quadratic=50, selected_proteins_file=None
    ).complete_dataframe

    # 1. build per-protein work items (this is what gets pickled to workers)
    t = time.perf_counter()
    spec = list(lfqprot.get_protein_workitems(normed_df, 10, 1))
    t_build = time.perf_counter() - t
    n = len(spec)

    # 2. map (parallel or sequential)
    t = time.perf_counter()
    if num_cores is not None and num_cores <= 1:
        results = lfqprot.get_list_with_sequential_processing(iter(spec))
    else:
        results = lfqprot.get_list_with_multiprocessing(iter(spec), num_cores)
    t_map = time.perf_counter() - t

    # 3. assemble protein dataframe
    t = time.perf_counter()
    protein_df = lfqprot.get_protein_dataframe_from_list_of_protein_profiles(
        list_of_tuple_w_protein_profiles_and_shifted_peptides=results,
        normed_df=normed_df,
    )
    t_prot = time.perf_counter() - t

    # 4. compile ion dataframe
    t = time.perf_counter()
    _ = lfqprot.get_ion_intensity_dataframe_from_list_of_shifted_peptides(
        results, column_names=normed_df.columns
    )
    t_ion = time.perf_counter() - t

    print(
        f"RESULT cores={num_cores} nprot={n} build_subdfs={t_build:.2f} "
        f"map={t_map:.2f} assemble_prot={t_prot:.2f} compile_ion={t_ion:.2f} "
        f"serial_overhead={t_build + t_prot + t_ion:.2f}"
    )


if __name__ == "__main__":
    main()
