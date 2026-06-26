"""Single-core optimization benchmark + correctness check for directLFQ.

- Loads + normalizes the HYE dataset once (serial prep).
- Runs estimate_protein_intensities(num_cores=1), timed.
- Compares protein + ion outputs against a saved reference (exact, equal_nan).
  First run (or `save`) writes the reference.

Usage:
    python optbench.py save     # generate reference from current code
    python optbench.py          # run, time, compare to reference
"""

import os
import sys
import time

import numpy as np

import directlfq.config as config
import directlfq.utils as lfqutils
import directlfq.normalization as lfqnorm
import directlfq.protein_intensity_estimation as lfqprot

INPUT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "input_data/HYE1-6_repeat/lfq.dia.features.csv.aq_reformat.tsv",
)
# Reference data lives next to this script so the benchmarks dir is self-contained.
_HERE = os.path.dirname(os.path.abspath(__file__))
REF_PROT = os.path.join(_HERE, "opt_protein.reference.npz")
REF_ION = os.path.join(_HERE, "opt_ion.reference.npz")


def prep():
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
    return lfqnorm.NormalizationManagerSamplesOnSelectedProteins(
        df, num_samples_quadratic=50, selected_proteins_file=None
    ).complete_dataframe


def protvalcutter_reference_sorted_index(protvals_df):
    """Reference: the original ProtvalCutter index sort, kept for the probes to compare
    the numpy _cut_peptide_values against now that ProtvalCutter is removed.

    Primary key: number of NaNs (ascending); secondary: summed intensity (descending).
    """
    idxs = protvals_df.index
    return sorted(
        idxs,
        key=lambda idx: (
            sum(
                np.isnan(protvals_df.loc[idx].to_numpy())
            ),  # First by number of NaNs (ascending)
            -np.nansum(
                protvals_df.loc[idx].to_numpy()
            ),  # Then by sum of intensities (descending)
        ),
    )


def protvalcutter_reference_cut(protvals_df, maximum_df_length=100):
    """Reference: the original ProtvalCutter cut -- keep the top maximum_df_length rows."""
    if len(protvals_df.index) <= maximum_df_length:
        return protvals_df
    shortened_index = protvalcutter_reference_sorted_index(protvals_df)[
        :maximum_df_length
    ]
    return protvals_df.loc[shortened_index]


def df_to_arrays(df):
    """Stable representation: id strings (from protein/ion id cols) + sorted float matrix."""
    df = df.reset_index()
    id_cols = [c for c in df.columns if c in ("protein", "ion")]
    val_cols = [c for c in df.columns if c not in id_cols and df[c].dtype.kind in "fiu"]
    ids = df[id_cols].astype(str).agg("\t".join, axis=1).to_numpy()
    order = np.argsort(ids, kind="stable")
    ids = ids[order]
    vals = df[val_cols].to_numpy(dtype=float)[order]
    return ids, vals


def save_ref(path, df):
    ids, vals = df_to_arrays(df)
    np.savez(path, ids=ids, vals=vals)


def compare(path, df, name):
    ref = np.load(path, allow_pickle=True)
    ids, vals = df_to_arrays(df)
    rids, rvals = ref["ids"], ref["vals"]
    if not np.array_equal(ids, rids):
        print(f"  [{name}] INDEX MISMATCH new={len(ids)} ref={len(rids)}")
        return False
    if vals.shape != rvals.shape:
        print(f"  [{name}] SHAPE MISMATCH new={vals.shape} ref={rvals.shape}")
        return False
    exact = np.array_equal(vals, rvals, equal_nan=True)
    nan_ok = np.array_equal(np.isnan(vals), np.isnan(rvals))
    d = np.abs(vals - rvals)
    d[np.isnan(d)] = 0
    maxabs = float(np.max(d)) if d.size else 0.0
    print(
        f"  [{name}] exact={exact} nan_pattern_ok={nan_ok} max_abs_diff={maxabs:.3e} n={len(ids)}"
    )
    return exact


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"

    normed = prep()

    t = time.perf_counter()
    protein_df, ion_df = lfqprot.estimate_protein_intensities(
        normed, min_nonan=1, num_samples_quadratic=10, num_cores=1
    )
    t_est = time.perf_counter() - t

    if mode == "save":
        save_ref(REF_PROT, protein_df)
        save_ref(REF_ION, ion_df)
        print(f"SAVED reference (estimate={t_est:.2f}s, nprot={len(protein_df)})")
        return

    print(f"TIME estimate={t_est:.2f}s  nprot={len(protein_df)}")
    ok_p = compare(REF_PROT, protein_df, "protein")
    ok_i = compare(REF_ION, ion_df, "ion")
    print(f"VERDICT identical={ok_p and ok_i}")


if __name__ == "__main__":
    main()
