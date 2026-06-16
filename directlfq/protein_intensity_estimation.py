__all__ = [
    "estimate_protein_intensities",
    "get_list_of_tuple_w_protein_profiles_and_shifted_peptides",
    "get_list_with_sequential_processing",
    "get_list_with_multiprocessing",
    "get_configured_multiprocessing_pool",
    "get_input_specification_tuplelist_idx__df__num_samples_quadratic__min_nonan",
    "get_normed_dfs",
    "get_ion_intensity_dataframe_from_list_of_shifted_peptides",
    "add_protein_names_to_ion_ints",
    "add_protein_name_to_ion_df",
    "get_protein_dataframe_from_list_of_protein_profiles",
    "calculate_peptide_and_protein_intensities",
    "get_protein_profile_from_shifted_peptides",
    "get_list_with_protein_value_for_each_sample",
    "ProtvalCutter",
    "OrphanIonRemover",
    "OrphanIonsForDeletionSelector",
    "IonCheckedForOrphan",
]

import pandas as pd
import numpy as np
import directlfq.normalization as lfqnorm
import multiprocess
import itertools
import logging
import warnings
import directlfq.config as config
from typing import Iterator, Optional

config.setup_logging()

LOGGER = logging.getLogger(__name__)


def estimate_protein_intensities(
    normed_df, min_nonan, num_samples_quadratic, num_cores
):
    """ "Performs the LFQ protein intensity estimation. The input is a normalized ion intensity dataframe. The output is the lfq protein intensity dataframe and the lfq ion intensity dataframe.

    Args:
        normed_df (pd.DataFrame): A pandas dataframe that contains the multi-index of [config.PROTEIN_ID, config.QUANT_ID]. The columns are the sample names and the values are the ion intensities. Zero value are replaced with NaNs and the dataframe is subsequently log-transformed. The dataframe needs to be sorted by the protein IDs, as the subsequent functions assume this.
        min_nonan (int): minimum number of NaNs
        num_samples_quadratic (int): minimum number of samples to use for the quadratic normalization
        num_cores (int): number of cores to use for the multiprocessing. If set to 1, the processing is done sequentially.

    Returns:
        tuple[protein_intensity_df, ion_intensity_df]: protein intensity dataframe and an ion intensity dataframe. The ion intensity dataframe is only compiled if the config.COMPILE_NORMALIZED_ION_TABLE is set to True.
    """

    allprots = list(normed_df.index.get_level_values(0).unique())
    LOGGER.info(f"{len(allprots)} lfq-groups total")

    list_of_tuple_w_protein_profiles_and_shifted_peptides = (
        get_list_of_tuple_w_protein_profiles_and_shifted_peptides(
            normed_df, num_samples_quadratic, min_nonan, num_cores
        )
    )
    protein_df = get_protein_dataframe_from_list_of_protein_profiles(
        list_of_tuple_w_protein_profiles_and_shifted_peptides=list_of_tuple_w_protein_profiles_and_shifted_peptides,
        normed_df=normed_df,
    )
    if config.COMPILE_NORMALIZED_ION_TABLE:
        ion_df = get_ion_intensity_dataframe_from_list_of_shifted_peptides(
            list_of_tuple_w_protein_profiles_and_shifted_peptides,
            column_names=normed_df.columns,
        )
    else:
        ion_df = None

    return protein_df, ion_df


def get_list_of_tuple_w_protein_profiles_and_shifted_peptides(
    normed_df, num_samples_quadratic, min_nonan, num_cores
):
    input_specification_tuplelist_idx__df__num_samples_quadratic__min_nonan = (
        get_input_specification_tuplelist_idx__df__num_samples_quadratic__min_nonan(
            normed_df, num_samples_quadratic, min_nonan
        )
    )

    if num_cores is not None and num_cores <= 1:
        list_of_tuple_w_protein_profiles_and_shifted_peptides = (
            get_list_with_sequential_processing(
                input_specification_tuplelist_idx__df__num_samples_quadratic__min_nonan
            )
        )
    else:
        list_of_tuple_w_protein_profiles_and_shifted_peptides = (
            get_list_with_multiprocessing(
                input_specification_tuplelist_idx__df__num_samples_quadratic__min_nonan,
                num_cores,
            )
        )
    return list_of_tuple_w_protein_profiles_and_shifted_peptides


def get_input_specification_tuplelist_idx__df__num_samples_quadratic__min_nonan(
    normed_df: pd.DataFrame, num_samples_quadratic: int, min_nonan: int
) -> Iterator[tuple[int, str, np.ndarray, np.ndarray, int, int]]:
    """Yield one work item per protein for the intensity-estimation workers.

    The normalized frame is split into contiguous numpy slices (one per protein,
    detected via ``find_nameswitch_indices`` on the protein-level index) so the
    hot path ships raw arrays instead of per-protein DataFrames. Each item is
    ``(idx, protein_name, ion_names, peptide_values, num_samples_quadratic,
    min_nonan)``, where ``peptide_values`` has rows for ions and columns for
    samples.
    """
    protein_names = normed_df.index.get_level_values(0).to_numpy()
    ion_names = normed_df.index.get_level_values(1).to_numpy()
    normed_array = normed_df.to_numpy()
    switch = find_nameswitch_indices(protein_names)
    n_proteins = len(switch) - 1
    return zip(
        range(n_proteins),
        (protein_names[switch[i]] for i in range(n_proteins)),
        (ion_names[switch[i] : switch[i + 1]] for i in range(n_proteins)),
        (normed_array[switch[i] : switch[i + 1]] for i in range(n_proteins)),
        itertools.repeat(num_samples_quadratic),
        itertools.repeat(min_nonan),
    )


def get_normed_dfs(normed_df):
    protein_names = normed_df.index.get_level_values(0).to_numpy()
    ion_names = normed_df.index.get_level_values(1).to_numpy()
    normed_array = normed_df.to_numpy()
    indices_of_proteinname_switch = find_nameswitch_indices(protein_names)
    results_list = [
        get_subdf(
            normed_array, indices_of_proteinname_switch, idx, protein_names, ion_names
        )
        for idx in range(len(indices_of_proteinname_switch) - 1)
    ]

    return results_list


def find_nameswitch_indices(arr):
    change_indices = np.where(arr[:-1] != arr[1:])[0] + 1

    # Add the index 0 for the start of the first element
    start_indices = np.insert(change_indices, 0, 0)

    # Append the index of the last element
    start_indices = np.append(start_indices, len(arr))

    return start_indices


def get_subdf(
    normed_array, indices_of_proteinname_switch, idx, protein_names, ion_names
):
    start_switch = indices_of_proteinname_switch[idx]
    end_switch = indices_of_proteinname_switch[idx + 1]
    sub_array = normed_array[start_switch:end_switch]
    index_sub_array = pd.MultiIndex.from_arrays(
        [protein_names[start_switch:end_switch], ion_names[start_switch:end_switch]],
        names=[config.PROTEIN_ID, config.QUANT_ID],
    )
    return pd.DataFrame(sub_array, index=index_sub_array)


def get_list_with_sequential_processing(
    input_specification_tuplelist_idx__df__num_samples_quadratic__min_nonan,
):
    list_of_tuple_w_protein_profiles_and_shifted_peptides = list(
        map(
            lambda x: calculate_peptide_and_protein_intensities(*x),
            input_specification_tuplelist_idx__df__num_samples_quadratic__min_nonan,
        )
    )
    return list_of_tuple_w_protein_profiles_and_shifted_peptides


def get_list_with_multiprocessing(
    input_specification_tuplelist_idx__df__num_samples_quadratic__min_nonan, num_cores
):
    pool = get_configured_multiprocessing_pool(num_cores)
    list_of_tuple_w_protein_profiles_and_shifted_peptides = pool.starmap(
        calculate_peptide_and_protein_intensities,
        input_specification_tuplelist_idx__df__num_samples_quadratic__min_nonan,
    )
    pool.close()
    return list_of_tuple_w_protein_profiles_and_shifted_peptides


def get_configured_multiprocessing_pool(num_cores):
    multiprocess.freeze_support()
    if num_cores is None:
        num_cores = (
            multiprocess.cpu_count() if multiprocess.cpu_count() < 60 else 60
        )  # windows upper thread limit
    pool = multiprocess.Pool(num_cores)
    LOGGER.info(f"using {pool._processes} processes")
    return pool


def calculate_peptide_and_protein_intensities(
    idx: int,
    protein_name: str,
    ion_names: np.ndarray,
    peptide_values: np.ndarray,
    num_samples_quadratic: int,
    min_nonan: int,
) -> tuple[Optional[np.ndarray], str, np.ndarray, np.ndarray]:
    """Compute one protein's profile and its shifted peptide values from numpy slices.

    Rows of ``peptide_values`` are ions, columns are samples. Proteins with more
    than 100 ions are first reduced via ``_cut_peptide_values``. Returns
    ``(protein_profile, protein_name, ion_names, shifted_values)``; the protein
    profile is ``None`` when every sample collapses to NaN.
    """
    if peptide_values.shape[0] > 1:
        peptide_values, ion_names = _cut_peptide_values(
            peptide_values, ion_names, maximum=100
        )

    if config.LOG_PROCESSED_PROTEINS and (
        idx % config.LOG_PROCESSED_PROTEINS_INTERVAL == 0
    ):
        LOGGER.info(f"lfq-object {idx}")
    # asfortranarray reproduces the column-major summation order of the old
    # np.nansum(2**peptide_intensity_df) path (pandas stores a homogeneous float
    # frame column-major) -> keeps summed_pepint bit-identical.
    summed_pepint = np.nansum(np.asfortranarray(2**peptide_values))

    shifted_values = _normalize_protein_values(peptide_values, num_samples_quadratic)
    protein_profile = get_protein_profile_from_shifted_peptides(
        shifted_values, summed_pepint, min_nonan
    )
    return protein_profile, protein_name, ion_names, shifted_values


def _cut_peptide_values(
    peptide_values: np.ndarray, ion_names: np.ndarray, maximum: int = 100
) -> tuple[np.ndarray, np.ndarray]:
    """Numpy equivalent of ProtvalCutter: keep at most ``maximum`` ions, sorted by
    NaN count asc then summed intensity desc. Only reorders when > maximum ions."""
    if peptide_values.shape[0] <= maximum:
        return peptide_values, ion_names
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        neg_summed = -np.nansum(peptide_values, axis=1)
    nan_counts = np.isnan(peptide_values).sum(axis=1)
    # lexsort's last key is primary and the sort is stable, so full ties keep
    # original order -> reproduces ProtvalCutter's sorted() tie-breaking.
    order = np.lexsort((neg_summed, nan_counts))[:maximum]
    return peptide_values[order], ion_names[order]


def _normalize_protein_values(
    peptide_values: np.ndarray, num_samples_quadratic: int
) -> np.ndarray:
    """Numpy equivalent of NormalizationManagerProtein. Rows are ions, cols samples."""
    if peptide_values.shape[0] <= num_samples_quadratic:
        values = peptide_values.copy()
        sample2shift = lfqnorm.get_normfacts(values)
        return lfqnorm.apply_sampleshifts(values, sample2shift)

    arr = peptide_values.copy()
    k = num_samples_quadratic
    nan_counts = np.isnan(arr).sum(axis=1)
    order = np.argsort(nan_counts, kind="stable")
    q_pos = order[:k]
    linear_mask = np.ones(arr.shape[0], dtype=bool)
    linear_mask[q_pos] = False
    lin_pos = np.flatnonzero(linear_mask)

    q_vals = arr[q_pos].copy()
    sample2shift = lfqnorm.get_normfacts(q_vals)
    q_normed = lfqnorm.apply_sampleshifts(q_vals, sample2shift)
    arr[q_pos] = q_normed

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        reference = np.nanmedian(q_normed, axis=0)
        if lin_pos.size:
            shifts = np.nanmedian(reference - arr[lin_pos], axis=1)
            arr[lin_pos] = arr[lin_pos] + shifts[:, None]
    return arr


def get_protein_profile_from_shifted_peptides(
    shifted_values: np.ndarray, summed_pepints: float, min_nonan: int
) -> Optional[np.ndarray]:
    """Collapse shifted peptide values into a per-sample protein profile.

    The per-sample nanmedian profile is rescaled so its summed linear intensity
    matches ``summed_pepints`` (the protein's total peptide intensity). Returns
    ``None`` when every sample is NaN.
    """
    intens_vec = get_list_with_protein_value_for_each_sample(shifted_values, min_nonan)
    intens_vec = np.array(intens_vec)
    summed_intensity = np.nansum(2**intens_vec)
    if summed_intensity == 0:  # this means all elements in intens vec are nans
        return None
    intens_conversion_factor = summed_pepints / summed_intensity
    return intens_vec + np.log2(intens_conversion_factor)


def get_list_with_protein_value_for_each_sample(
    shifted_values: np.ndarray, min_nonan: int
) -> np.ndarray:
    """Collapse a protein's normalized peptide profiles into one value per sample.

    Each sample (column) is summarized by the nanmedian of its ion values. Samples
    with fewer than ``min_nonan`` finite ions are set to NaN.
    """
    nonan_counts = np.sum(~np.isnan(shifted_values), axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter(
            "ignore", category=RuntimeWarning
        )  # all-NaN columns -> NaN
        intens_vec = np.nanmedian(shifted_values, axis=0)
    intens_vec[nonan_counts < min_nonan] = np.nan
    return intens_vec


import pandas as pd
from numba import njit


class ProtvalCutter:
    def __init__(self, protvals_df, maximum_df_length=100):
        self._protvals_df = protvals_df
        self._maximum_df_length = maximum_df_length
        self._dataframe_too_long = None
        self._sorted_idx = None
        self._check_if_df_too_long_and_sort_index_if_so()

    def _check_if_df_too_long_and_sort_index_if_so(self):
        self._dataframe_too_long = (
            len(self._protvals_df.index) > self._maximum_df_length
        )
        if self._dataframe_too_long:
            self._determine_nansorted_df_index()

    def _determine_nansorted_df_index(self):
        """Sorts the dataframe index primarily by number of NaN values (ascending) and secondarily by summed intensity (descending). Sorting by intensties in case multiple ions have identical missing value counts. We expect initial sorting by ion name (which is done in the run_lfq module) to be deterministic.

        The sorting prioritizes:
        1. Rows with fewer NaN values come first
        2. For rows with equal number of NaNs, higher intensity sums come first
        """
        idxs = self._protvals_df.index
        self._sorted_idx = sorted(
            idxs,
            key=lambda idx: (
                sum(
                    np.isnan(self._protvals_df.loc[idx].to_numpy())
                ),  # First by number of NaNs (ascending)
                -np.nansum(
                    self._protvals_df.loc[idx].to_numpy()
                ),  # Then by sum of intensities (descending)
            ),
        )

    @staticmethod
    @njit
    def _get_num_nas_in_row(row):
        sum = 0
        isnans = np.isnan(row)
        for is_nan in isnans:
            sum += is_nan
        return sum

    def get_dataframe(self):
        if self._dataframe_too_long:
            return self._get_shortened_dataframe()
        else:
            return self._protvals_df

    def _get_shortened_dataframe(self):
        shortened_index = self._sorted_idx[: self._maximum_df_length]
        return self._protvals_df.loc[shortened_index]


def get_ion_intensity_dataframe_from_list_of_shifted_peptides(
    list_of_tuple_w_protein_profiles_and_shifted_peptides: list[
        tuple[Optional[np.ndarray], str, np.ndarray, np.ndarray]
    ],
    column_names: list[str],
) -> pd.DataFrame:
    """Assemble the per-ion intensity table from the workers' result tuples.

    Reads the ``(protein_profile, protein_name, ion_names, shifted_values)`` shape,
    converts the shifted log2 values back to linear space (NaNs to 0), and returns
    a DataFrame indexed by ``(protein, ion)`` with one column per sample.
    """
    ion_names = []
    ion_vals = []
    protein_names = []
    for (
        _,
        protein_name,
        inames,
        shifted_values,
    ) in list_of_tuple_w_protein_profiles_and_shifted_peptides:
        ion_names.extend(inames.tolist())
        ion_vals.append(shifted_values)
        protein_names.extend([protein_name] * len(inames))
    merged_ions = np.nan_to_num(2 ** np.concatenate(ion_vals))
    ion_df = pd.DataFrame(merged_ions)
    ion_df.columns = column_names
    ion_df["ion"] = ion_names
    ion_df["protein"] = protein_names
    return ion_df.set_index(["protein", "ion"])


def add_protein_names_to_ion_ints(ion_ints, allprots):
    ion_ints = [
        add_protein_name_to_ion_df(ion_ints[idx], allprots[idx])
        for idx in range(len(ion_ints))
    ]
    return ion_ints


def add_protein_name_to_ion_df(ion_df, protein):
    ion_df[config.PROTEIN_ID] = protein
    ion_df = ion_df.reset_index().set_index([config.PROTEIN_ID, config.QUANT_ID])
    return ion_df


def get_protein_dataframe_from_list_of_protein_profiles(
    list_of_tuple_w_protein_profiles_and_shifted_peptides: list[
        tuple[Optional[np.ndarray], str, np.ndarray, np.ndarray]
    ],
    normed_df: pd.DataFrame,
) -> pd.DataFrame:
    """Assemble the per-protein intensity table from the workers' result tuples.

    Reads the ``(protein_profile, protein_name, ion_names, shifted_values)`` shape,
    drops proteins whose profile is ``None``, converts the log2 profiles back to
    linear space (NaNs to 0), and returns a DataFrame with a ``protein`` column and
    one column per sample.
    """
    index_list = []
    profile_list = []

    list_of_protein_profiles = [
        x[0] for x in list_of_tuple_w_protein_profiles_and_shifted_peptides
    ]
    allprots = [x[1] for x in list_of_tuple_w_protein_profiles_and_shifted_peptides]

    for idx in range(len(allprots)):
        if list_of_protein_profiles[idx] is None:
            continue
        index_list.append(allprots[idx])
        profile_list.append(list_of_protein_profiles[idx])

    index_for_protein_df = pd.Index(data=index_list, name=config.PROTEIN_ID)
    protein_df = 2 ** pd.DataFrame(
        profile_list, index=index_for_protein_df, columns=normed_df.columns
    )
    protein_df = protein_df.replace(np.nan, 0)
    protein_df = protein_df.reset_index()
    return protein_df
