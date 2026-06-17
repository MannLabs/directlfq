__all__ = [
    "estimate_protein_intensities",
    "get_list_of_tuple_w_protein_profiles_and_shifted_peptides",
    "get_list_with_sequential_processing",
    "get_list_with_multiprocessing",
    "get_configured_multiprocessing_pool",
    "get_protein_workitems",
    "get_ion_intensity_dataframe_from_list_of_shifted_peptides",
    "add_protein_names_to_ion_ints",
    "add_protein_name_to_ion_df",
    "get_protein_dataframe_from_list_of_protein_profiles",
    "calculate_peptide_and_protein_intensities",
    "get_protein_profile_from_shifted_peptides",
    "get_list_with_protein_value_for_each_sample",
]

import pandas as pd
import numpy as np
import directlfq.normalization as lfqnorm
import multiprocess
import logging
import warnings
import directlfq.config as config

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
    protein_workitems = get_protein_workitems(
        normed_df, num_samples_quadratic, min_nonan
    )

    if num_cores is not None and num_cores <= 1:
        list_of_tuple_w_protein_profiles_and_shifted_peptides = (
            get_list_with_sequential_processing(protein_workitems)
        )
    else:
        list_of_tuple_w_protein_profiles_and_shifted_peptides = (
            get_list_with_multiprocessing(protein_workitems, num_cores)
        )
    return list_of_tuple_w_protein_profiles_and_shifted_peptides


def get_protein_workitems(normed_df, num_samples_quadratic, min_nonan):
    """Yield one work item per protein as contiguous numpy slices (no per-protein DataFrame).

    ``normed_df`` is sorted by (protein, ion), so each protein is a contiguous row block.
    Each item is the 6-tuple expected positionally by
    ``calculate_peptide_and_protein_intensities``; ``peptide_values`` has rows = ions,
    columns = samples.
    """
    protein_names = normed_df.index.get_level_values(0).to_numpy()
    ion_names = normed_df.index.get_level_values(1).to_numpy()
    normed_array = normed_df.to_numpy()
    switch_indices = find_nameswitch_indices(protein_names)
    for idx in range(len(switch_indices) - 1):
        start = switch_indices[idx]
        end = switch_indices[idx + 1]
        yield (
            idx,
            protein_names[start],
            ion_names[start:end],
            normed_array[start:end],
            num_samples_quadratic,
            min_nonan,
        )


def find_nameswitch_indices(arr):
    change_indices = np.where(arr[:-1] != arr[1:])[0] + 1

    # Add the index 0 for the start of the first element
    start_indices = np.insert(change_indices, 0, 0)

    # Append the index of the last element
    start_indices = np.append(start_indices, len(arr))

    return start_indices


def get_list_with_sequential_processing(protein_workitems):
    list_of_tuple_w_protein_profiles_and_shifted_peptides = list(
        map(
            lambda x: calculate_peptide_and_protein_intensities(*x),
            protein_workitems,
        )
    )
    return list_of_tuple_w_protein_profiles_and_shifted_peptides


def get_list_with_multiprocessing(protein_workitems, num_cores):
    pool = get_configured_multiprocessing_pool(num_cores)
    list_of_tuple_w_protein_profiles_and_shifted_peptides = pool.starmap(
        calculate_peptide_and_protein_intensities,
        protein_workitems,
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
    idx, protein_name, ion_names, peptide_values, num_samples_quadratic, min_nonan
):
    if peptide_values.shape[0] > 1:
        peptide_values, ion_names = _cut_peptide_values(
            peptide_values, ion_names, maximum=100
        )

    if config.LOG_PROCESSED_PROTEINS and (
        idx % config.LOG_PROCESSED_PROTEINS_INTERVAL == 0
    ):
        LOGGER.info(f"lfq-object {idx}")
    # Fortran order matches the column-major sum of the former homogeneous DataFrame (Gotcha 1)
    summed_pepint = np.nansum(np.asfortranarray(2**peptide_values))

    if peptide_values.shape[1] < 2:
        shifted_values = peptide_values  # single sample: skip normalization, keep values as-is
    else:
        shifted_values = lfqnorm.normalize_protein_ion_values(
            peptide_values, num_samples_quadratic
        )

    protein_profile = get_protein_profile_from_shifted_peptides(
        shifted_values, summed_pepint, min_nonan
    )

    return protein_profile, protein_name, ion_names, shifted_values


def get_protein_profile_from_shifted_peptides(
    normalized_peptide_profiles, summed_pepints, min_nonan
):
    intens_vec = get_list_with_protein_value_for_each_sample(
        normalized_peptide_profiles, min_nonan
    )
    intens_vec = np.array(intens_vec)
    summed_intensity = np.nansum(2**intens_vec)
    if summed_intensity == 0:  # this means all elements in intens vec are nans
        return None
    intens_conversion_factor = summed_pepints / summed_intensity
    scaled_vec = intens_vec + np.log2(intens_conversion_factor)
    return scaled_vec


def get_list_with_protein_value_for_each_sample(
    normalized_peptide_profiles: np.ndarray, min_nonan: int
) -> np.ndarray:
    """Collapse a protein's normalized peptide profiles into one value per sample.

    Each sample (column) is summarized by the nanmedian of its ion values. Samples
    with fewer than ``min_nonan`` finite ions are set to NaN. Input rows are ions,
    columns are samples.
    """
    arr = normalized_peptide_profiles
    nonan_counts = np.sum(~np.isnan(arr), axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter(
            "ignore", category=RuntimeWarning
        )  # all-NaN columns -> NaN
        intens_vec = np.nanmedian(arr, axis=0)
    intens_vec[nonan_counts < min_nonan] = np.nan
    return intens_vec


def _cut_peptide_values(peptide_values, ion_names, maximum=100):
    """Reduce a protein to its ``maximum`` most informative ions (rows = ions).

    Keeps ions sorted by NaN count ascending, then summed intensity descending --
    the numpy equivalent of ``ProtvalCutter``. Reproduces ``ProtvalCutter``'s stable
    ``sorted()`` tie-break: ions tied on both keys keep their original (ion-name) order.
    """
    if peptide_values.shape[0] <= maximum:
        return peptide_values, ion_names
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)  # all-NaN rows -> NaN
        neg_summed = -np.nansum(peptide_values, axis=1)
    nan_counts = np.isnan(peptide_values).sum(axis=1)
    # last key is primary and lexsort is stable -> full ties keep original order (Gotcha 2)
    order = np.lexsort((neg_summed, nan_counts))[:maximum]
    return peptide_values[order], ion_names[order]


def get_ion_intensity_dataframe_from_list_of_shifted_peptides(
    list_of_tuple_w_protein_profiles_and_shifted_peptides, column_names
):
    ion_names = []
    ion_vals = []
    protein_names = []
    for idx in range(len(list_of_tuple_w_protein_profiles_and_shifted_peptides)):
        _, protein_name, ion_names_arr, shifted_values = (
            list_of_tuple_w_protein_profiles_and_shifted_peptides[idx]
        )
        ion_names += ion_names_arr.tolist()
        ion_vals.append(shifted_values)
        protein_names.extend([protein_name] * len(ion_names_arr))
    merged_ions = 2 ** np.concatenate(ion_vals)
    merged_ions = np.nan_to_num(merged_ions)
    ion_df = pd.DataFrame(merged_ions)
    ion_df.columns = column_names
    ion_df["ion"] = ion_names
    ion_df["protein"] = protein_names
    ion_df = ion_df.set_index(["protein", "ion"])
    return ion_df


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
    list_of_tuple_w_protein_profiles_and_shifted_peptides, normed_df
):
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
