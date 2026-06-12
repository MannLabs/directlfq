"""Stage-resolved benchmark of the directLFQ pipeline.

Runs the same steps as lfq_manager.run_lfq but times the three big phases
separately so we can see which part actually scales with num_cores:

    import   : read TSV + sort + dedup + log-transform   (serial)
    norm     : between-sample normalization               (serial)
    estimate : per-protein LFQ estimation                 (parallel, uses num_cores)

Run as a fresh process per core count so numba JIT cost is paid every time
(this is what a real user experiences when invoking directlfq from the CLI).

Usage: python bench_dlfq.py <num_cores>
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
    config.set_log_processed_proteins(log_processed_proteins=False)  # cut per-protein log noise
    config.set_compile_normalized_ion_table(compile_normalized_ion_table=True)
    config.check_wether_to_copy_numpy_arrays_derived_from_pandas()

    wall0 = time.perf_counter()

    t = time.perf_counter()
    inp = lfqutils.add_mq_protein_group_ids_if_applicable_and_obtain_annotated_file(
        INPUT, None, None, []
    )
    df = lfqutils.import_data(input_file=inp, input_type_to_use=None, filter_dict=None)
    df = lfqutils.sort_input_df_by_protein_and_quant_id(df)
    df = lfqutils.remove_potential_quant_id_duplicates(df)
    df = lfqutils.index_and_log_transform_input_df(df)
    df = lfqutils.remove_allnan_rows_input_df(df)
    t_import = time.perf_counter() - t

    t = time.perf_counter()
    df = lfqnorm.NormalizationManagerSamplesOnSelectedProteins(
        df, num_samples_quadratic=50, selected_proteins_file=None
    ).complete_dataframe
    t_norm = time.perf_counter() - t

    t = time.perf_counter()
    protein_df, ion_df = lfqprot.estimate_protein_intensities(
        df, min_nonan=1, num_samples_quadratic=10, num_cores=num_cores
    )
    t_est = time.perf_counter() - t

    wall = time.perf_counter() - wall0
    print(
        f"RESULT cores={num_cores} import={t_import:.2f} norm={t_norm:.2f} "
        f"estimate={t_est:.2f} wall={wall:.2f} nprot={len(protein_df)}"
    )


if __name__ == "__main__":
    main()
