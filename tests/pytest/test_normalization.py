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


# ============================================================================
# calc_nanmedian / calc_nanvar numba kernels (optimization step 7)
# ============================================================================
# These kernels gain @njit(cache=True) in step 7 (persisting the one-time JIT).
# Caching is numerically inert, so the kernels must keep matching numpy.


def test_calc_nanmedian_matches_numpy_skipping_nans():
    arr = np.array([1.0, np.nan, 3.0, 5.0])
    assert lfq_norm.calc_nanmedian(arr) == np.nanmedian(arr)


def test_calc_nanvar_matches_numpy_skipping_nans():
    arr = np.array([1.0, np.nan, 3.0, 5.0])
    assert lfq_norm.calc_nanvar(arr) == pytest.approx(np.nanvar(arr))


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


# ============================================================================
# set_samples_with_only_single_intensity_to_nan (optimization step 6)
# ============================================================================
# Contract locked in before replacing the per-row Python sum() loop with
# np.count_nonzero(axis=1): in place, any row with fewer than 2 finite values is
# set entirely to NaN; a row with exactly 2 finite values is kept.


def test_set_samples_with_only_single_intensity_to_nan_masks_rows_below_two():
    # given - rows with 3, 1, 0 and exactly 2 finite values
    samples = np.array(
        [
            [1.0, 2.0, 3.0],
            [1.0, np.nan, np.nan],
            [np.nan, np.nan, np.nan],
            [1.0, 2.0, np.nan],
        ]
    )

    # when - mutates in place
    lfq_norm.set_samples_with_only_single_intensity_to_nan(samples)

    # then - only the single-/zero-intensity rows are blanked
    expected = np.array(
        [
            [1.0, 2.0, 3.0],
            [np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan],
            [1.0, 2.0, np.nan],
        ]
    )
    assert np.array_equal(samples, expected, equal_nan=True)


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


# ============================================================================
# NormalizationManager._determine_sorted_rows (optimization step 2)
# ============================================================================
# Contract locked in before replacing the per-row .loc sort key with a numpy
# nan-count + stable argsort: index labels ordered by ascending NaN-count, with
# ties broken stably (original order preserved for equal counts).


def _norm_manager_with_rows(rows):
    """NormalizationManager over rows (ions) with a (protein, ion) MultiIndex."""
    index = pd.MultiIndex.from_tuples(
        [("P", f"ion{i}") for i in range(len(rows))], names=["protein", "ion"]
    )
    df = pd.DataFrame(rows, index=index)
    return lfq_norm.NormalizationManager(df, num_samples_quadratic=100)


def test_determine_sorted_rows_orders_by_ascending_nan_count_stably():
    # given - rows with nan-counts 0, 2, 1, 3, 1 (the two 1s are a tie)
    rows = [
        [1.0, 2.0, 3.0],
        [1.0, np.nan, np.nan],
        [1.0, 2.0, np.nan],
        [np.nan, np.nan, np.nan],
        [4.0, 5.0, np.nan],
    ]
    mgr = _norm_manager_with_rows(rows)

    # when
    mgr._determine_sorted_rows()

    # then - ascending nan-count; ion2 before ion4 (stable tie-break)
    assert mgr._rows_sorted_by_number_valid_values == [
        ("P", "ion0"),
        ("P", "ion2"),
        ("P", "ion4"),
        ("P", "ion1"),
        ("P", "ion3"),
    ]


def test_determine_sorted_rows_preserves_original_order_when_no_nans():
    # given - all rows fully finite (all keys equal -> stable order is identity)
    rows = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
    mgr = _norm_manager_with_rows(rows)

    # when
    mgr._determine_sorted_rows()

    # then
    assert mgr._rows_sorted_by_number_valid_values == [
        ("P", "ion0"),
        ("P", "ion1"),
        ("P", "ion2"),
    ]


# ============================================================================
# NormalizationManagerProtein._normalize_quadratic_and_linear (optimization step 4)
# ============================================================================
# Contract locked in before replacing the label-based .loc[tuples,:]=df
# round-trips with numpy. For proteins with more than num_samples_quadratic
# ions: the k rows with fewest NaNs (stable) are the quadratic subset normalized
# via get_normfacts/apply_sampleshifts; the remaining rows are each shifted onto
# the column-median of the normalized subset by nanmedian(reference - row).


def _normalize_quadratic_and_linear_oracle(complete_dataframe, num_samples_quadratic):
    """Executable spec via the public primitives. Validated against the baseline
    pandas implementation, then used to guard the numpy override."""
    arr = complete_dataframe.to_numpy(dtype=float, copy=True)
    k = num_samples_quadratic
    nan_counts = np.isnan(arr).sum(axis=1)
    order = np.argsort(nan_counts, kind="stable")
    q_pos = order[:k]
    linear_mask = np.ones(arr.shape[0], dtype=bool)
    linear_mask[q_pos] = False
    lin_pos = np.flatnonzero(linear_mask)
    q_vals = arr[q_pos].copy()
    sample2shift = lfq_norm.get_normfacts(q_vals)
    q_normed = lfq_norm.apply_sampleshifts(q_vals, sample2shift)
    arr[q_pos] = q_normed
    reference = np.nanmedian(q_normed, axis=0)
    shifts = np.nanmedian(reference - arr[lin_pos], axis=1)
    arr[lin_pos] = arr[lin_pos] + shifts[:, None]
    return arr


def test_normalize_quadratic_and_linear_matches_spec_on_mixed_nan_rows():
    # given - 5 ions x 4 samples; nan-counts 0,0,1,0,2 -> quadratic = rows 0,1,3
    df = _create_input_df_from_input_vals(
        [
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 3.0, 4.0, 5.0],
            [10.0, np.nan, 12.0, 13.0],
            [3.0, 4.0, 5.0, 6.0],
            [np.nan, np.nan, 20.0, 21.0],
        ]
    )
    expected = _normalize_quadratic_and_linear_oracle(df, num_samples_quadratic=3)

    # when
    result = lfq_norm.NormalizationManagerProtein(
        df.copy(), num_samples_quadratic=3
    ).complete_dataframe

    # then - values match the spec exactly and index/columns are preserved
    assert np.array_equal(result.to_numpy(), expected, equal_nan=True)
    assert list(result.index) == list(df.index)
    assert list(result.columns) == list(df.columns)


def test_normalize_quadratic_and_linear_overlaps_noisefree_profiles():
    # given - 6 noise-free peptides (exact scaled copies) forcing the
    # quadratic+linear path with num_samples_quadratic=3
    peptides = [
        lfq_test_utils.PeptideProfile(
            protein_name="protA",
            fraction_zeros_in_profile=0,
            systematic_peptide_shift=shift,
            add_noise=False,
        )
        for shift in [1, 4, 16, 64, 256, 1024]
    ]
    protein_df = lfq_test_utils.ProteinProfileGenerator(
        peptides
    ).protein_profile_dataframe

    # when
    normed = lfq_norm.NormalizationManagerProtein(
        protein_df, num_samples_quadratic=3
    ).complete_dataframe

    # then - all ions collapse onto a single profile (each column constant)
    values = normed.to_numpy()
    assert np.allclose(values, values[0])


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


# ============================================================================
# SampleShifterLinear._shift_columns_to_reference_sample (optimization step 3)
# ============================================================================
# Contract locked in before replacing the per-row `iloc[r,:] +=` shift with a
# single numpy block write: every row is shifted by nanmedian(reference - row),
# and an all-NaN row stays all-NaN.


def test_sample_shifter_shifts_each_row_by_nanmedian_distance():
    # given - rows that are fully finite, partially NaN, and fully NaN
    ion_dataframe = pd.DataFrame(
        [[1.0, 2.0, 3.0], [5.0, np.nan, 5.0], [np.nan, np.nan, np.nan]]
    )
    reference = pd.Series([10.0, 20.0, 30.0])

    # when
    shifted = lfq_norm.SampleShifterLinear(ion_dataframe, reference).ion_dataframe

    # then - row shifts are nanmedian([9,18,27])=18, nanmedian([5,25])=15, NaN
    expected = np.array(
        [[19.0, 20.0, 21.0], [20.0, np.nan, 20.0], [np.nan, np.nan, np.nan]]
    )
    assert np.array_equal(shifted.to_numpy(), expected, equal_nan=True)


def test_sample_shifter_leaves_row_unchanged_when_already_on_reference():
    # given - a row identical to the reference (distance 0)
    ion_dataframe = pd.DataFrame([[10.0, 20.0, 30.0]])
    reference = pd.Series([10.0, 20.0, 30.0])

    # when
    shifted = lfq_norm.SampleShifterLinear(ion_dataframe, reference).ion_dataframe

    # then
    assert np.array_equal(shifted.to_numpy(), np.array([[10.0, 20.0, 30.0]]))
