"""Tests for directlfq.normalization, salvaged from nbdev_nbs/02_normalization.ipynb."""

import os

import numpy as np
import pandas as pd
import pytest

import directlfq.normalization as lfq_norm
import directlfq.test_utils as lfq_test_utils
import directlfq.utils as lfq_utils


def _create_input_df_from_input_vals(list_of_vals):
    index_vals = [("A", f"ion{x}") for x in range(len(list_of_vals))]
    index = pd.Index(index_vals, name=("protein", "ion"))
    return pd.DataFrame(list_of_vals, index=index)


def _generate_randarrays(number_arrays, size_of_array):
    randarray = []
    for _ in range(number_arrays):
        shift = np.random.uniform(low=-10, high=+10)
        randarray.append(np.random.normal(loc=shift, size=size_of_array))
    return np.array(randarray)


def test_merged_distribs():
    anchor_distrib = np.array([1, 1, 1, 1, 1])
    shift_distrib = np.array([2, 2, 2, 2, 2])
    counts_anchor_distrib = 4
    counts_shifted_distib = 1
    assert (
        lfq_norm.merge_distribs(
            anchor_distrib, shift_distrib, counts_anchor_distrib, counts_shifted_distib
        )
        == np.array([1.2, 1.2, 1.2, 1.2, 1.2])
    ).any()


def test_order_of_shifts():
    vals1 = [1, np.nan, 1.5]
    vals2 = [1, 1, np.nan]
    vals3 = [3.2, 1, 2.8]
    vals4 = [4.2, 2, 3.8]
    list_of_vals = [vals1, vals2, vals3, vals4]
    protein_profile_df = _create_input_df_from_input_vals(list_of_vals)
    protein_profile_numpy = protein_profile_df.to_numpy()
    sample2shift = lfq_norm.get_normfacts(protein_profile_numpy)
    assert sample2shift == {1: 0.0, 2: -1.2999999999999998, 3: -2.3}


def test_sampleshift():
    np.random.seed(42)
    randarray = _generate_randarrays(5, 1000)
    sample2shift = lfq_norm.get_normfacts(randarray)
    normalized_randarray = lfq_norm.apply_sampleshifts(randarray, sample2shift)

    merged_sample = []
    for i in range(normalized_randarray.shape[0]):
        merged_sample.extend(normalized_randarray[i])
    stdev = np.std(merged_sample)
    assert stdev <= 1.2


def _assert_that_results_scatter_around_zero(input_df_normalized):
    median_intensities = input_df_normalized.median(axis=1)
    input_df_subtracted = input_df_normalized.subtract(median_intensities, axis=0)
    median_of_medians = input_df_subtracted.median(axis=0)
    assert (median_of_medians < 0.1).all()


@pytest.mark.parametrize("num_samples_quadratic", [100, 3, 1])
def test_normalizing_between_samples(test_data_dir, num_samples_quadratic):
    input_file = os.path.join(
        test_data_dir,
        "unit_tests",
        "protein_normalization",
        "peptides.txt.maxquant_peptides_benchmarking.aq_reformat.tsv",
    )
    input_df = pd.read_csv(input_file, sep="\t")
    input_df = lfq_utils.index_and_log_transform_input_df(input_df)
    input_df = input_df[[x for x in input_df.columns if "Shotgun" in x]]
    input_df_normalized = lfq_norm.NormalizationManagerSamples(
        input_df, num_samples_quadratic=num_samples_quadratic
    ).complete_dataframe
    _assert_that_results_scatter_around_zero(input_df_normalized)


def test_that_profiles_without_noise_are_shifted_exactly_on_top_of_each_other():
    peptide1 = lfq_test_utils.PeptideProfile(
        protein_name="protA",
        fraction_zeros_in_profile=0.1,
        systematic_peptide_shift=3000,
        add_noise=False,
    )
    peptide2 = lfq_test_utils.PeptideProfile(
        protein_name="protA",
        fraction_zeros_in_profile=0.9,
        systematic_peptide_shift=3,
        add_noise=False,
    )
    peptide3 = lfq_test_utils.PeptideProfile(
        protein_name="protA",
        fraction_zeros_in_profile=0.1,
        systematic_peptide_shift=0.1,
        add_noise=False,
    )
    peptide4 = lfq_test_utils.PeptideProfile(
        protein_name="protA",
        fraction_zeros_in_profile=0.9,
        systematic_peptide_shift=100,
        add_noise=False,
    )
    protein_df = lfq_test_utils.ProteinProfileGenerator(
        [peptide1, peptide2, peptide3, peptide4]
    ).protein_profile_dataframe
    normed_ion_profile = lfq_norm.normalize_ion_profiles(protein_df)
    column_from_shifted = normed_ion_profile.iloc[:, 11].dropna().to_numpy()
    assert np.allclose(column_from_shifted, column_from_shifted[0])


def test_that_profiles_with_noise_are_close():
    peptide1 = lfq_test_utils.PeptideProfile(
        protein_name="protA",
        fraction_zeros_in_profile=0,
        systematic_peptide_shift=3000,
        add_noise=True,
    )
    peptide2 = lfq_test_utils.PeptideProfile(
        protein_name="protA",
        fraction_zeros_in_profile=0,
        systematic_peptide_shift=3,
        add_noise=True,
    )
    peptide3 = lfq_test_utils.PeptideProfile(
        protein_name="protA",
        fraction_zeros_in_profile=0,
        systematic_peptide_shift=0.1,
        add_noise=True,
    )
    peptide4 = lfq_test_utils.PeptideProfile(
        protein_name="protA",
        fraction_zeros_in_profile=0,
        systematic_peptide_shift=100,
        add_noise=True,
    )
    protein_df = lfq_test_utils.ProteinProfileGenerator(
        [peptide1, peptide2, peptide3, peptide4]
    ).protein_profile_dataframe
    normed_ion_profile = lfq_norm.normalize_ion_profiles(protein_df)
    column_from_shifted = normed_ion_profile.iloc[:, 9].dropna().to_numpy()
    assert np.allclose(
        column_from_shifted, column_from_shifted[0], rtol=0.01, atol=0.01
    )


def _calc_distance(samples_1, samples_2):
    distrib = lfq_norm.get_fcdistrib(samples_1, samples_2)
    is_all_nan = np.all(np.isnan(distrib))
    if is_all_nan:
        return np.nan
    return np.nanmedian(distrib)


def test_calc_distance():
    # One array is entirely NaN
    samples_1 = np.array([np.nan, np.nan, np.nan])
    samples_2 = np.array([1, 2, 3])
    assert np.isnan(lfq_norm.SampleShifterLinear._calc_distance(samples_1, samples_2))

    # Both arrays are non-NaN and identical
    samples_1 = np.array([1, 2, 3])
    samples_2 = np.array([1, 2, 3])
    assert not np.isnan(_calc_distance(samples_1, samples_2))
    assert lfq_norm.SampleShifterLinear._calc_distance(samples_1, samples_2) == 0

    # Arrays with some NaN values
    samples_1 = np.array([1, np.nan, 3])
    samples_2 = np.array([13, 2, np.nan])
    assert not np.isnan(_calc_distance(samples_1, samples_2))
    assert lfq_norm.SampleShifterLinear._calc_distance(samples_1, samples_2) == -12

    # Arrays with different values but no NaNs
    samples_1 = np.array([1, 4, 7])
    samples_2 = np.array([2, 5, 8])
    assert lfq_norm.SampleShifterLinear._calc_distance(samples_1, samples_2) != 0

    # Empty arrays
    samples_1 = np.array([])
    samples_2 = np.array([])
    assert np.isnan(lfq_norm.SampleShifterLinear._calc_distance(samples_1, samples_2))
