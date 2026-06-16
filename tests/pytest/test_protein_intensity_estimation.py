"""Tests for directlfq.protein_intensity_estimation, salvaged from nbdev_nbs/03_protein_intensity_estimation.ipynb."""

import numpy as np
import pandas as pd
import pytest

import directlfq.normalization as lfq_norm
import directlfq.protein_intensity_estimation as lfq_protint
import directlfq.test_utils as lfq_testutils


def test_that_protein_intensities_are_retained():
    peptide1 = lfq_testutils.PeptideProfile(
        protein_name="protA",
        fraction_zeros_in_profile=0.1,
        systematic_peptide_shift=3000,
        add_noise=True,
    )
    peptide2 = lfq_testutils.PeptideProfile(
        protein_name="protA",
        fraction_zeros_in_profile=0,
        systematic_peptide_shift=3,
        add_noise=True,
    )
    peptide3 = lfq_testutils.PeptideProfile(
        protein_name="protA",
        fraction_zeros_in_profile=0,
        systematic_peptide_shift=0.1,
        add_noise=True,
    )
    peptide4 = lfq_testutils.PeptideProfile(
        protein_name="protA",
        fraction_zeros_in_profile=0,
        systematic_peptide_shift=100,
        add_noise=True,
    )
    peptide_profiles = [peptide1, peptide2, peptide3, peptide4]
    summed_intensity_protein = sum(
        [np.nansum(x.peptide_profile_vector) for x in peptide_profiles]
    )

    protein_df = lfq_testutils.ProteinProfileGenerator(
        peptide_profiles
    ).protein_profile_dataframe
    protein_df_normed, _ = lfq_protint.estimate_protein_intensities(
        protein_df, min_nonan=1, num_samples_quadratic=100, num_cores=1
    )
    summed_lfq_intensities = np.sum(protein_df_normed.iloc[0, 1:].to_numpy())
    assert np.allclose(summed_lfq_intensities, summed_intensity_protein)


def test_that_multiple_proteins_are_handled():
    peptide1 = lfq_testutils.PeptideProfile(
        protein_name="protA",
        fraction_zeros_in_profile=0.1,
        systematic_peptide_shift=3000,
        add_noise=True,
    )
    peptide2 = lfq_testutils.PeptideProfile(
        protein_name="protA",
        fraction_zeros_in_profile=0,
        systematic_peptide_shift=3,
        add_noise=True,
    )
    peptide3 = lfq_testutils.PeptideProfile(
        protein_name="protA",
        fraction_zeros_in_profile=0,
        systematic_peptide_shift=0.1,
        add_noise=True,
    )
    peptide4 = lfq_testutils.PeptideProfile(
        protein_name="protB",
        fraction_zeros_in_profile=0,
        systematic_peptide_shift=100,
        add_noise=True,
    )
    peptide5 = lfq_testutils.PeptideProfile(
        protein_name="protC",
        fraction_zeros_in_profile=0,
        systematic_peptide_shift=100,
        add_noise=True,
    )
    peptide6 = lfq_testutils.PeptideProfile(
        protein_name="protD",
        fraction_zeros_in_profile=0,
        systematic_peptide_shift=100,
        add_noise=True,
    )
    peptide7 = lfq_testutils.PeptideProfile(
        protein_name="protD",
        fraction_zeros_in_profile=0,
        systematic_peptide_shift=100,
        add_noise=True,
    )
    peptide8 = lfq_testutils.PeptideProfile(
        protein_name="protD",
        fraction_zeros_in_profile=0,
        systematic_peptide_shift=100,
        add_noise=True,
    )

    peptide_profiles = [
        peptide1,
        peptide2,
        peptide3,
        peptide4,
        peptide5,
        peptide6,
        peptide7,
        peptide8,
    ]
    protein_df = lfq_testutils.ProteinProfileGenerator(
        peptide_profiles
    ).protein_profile_dataframe
    protein_df_normed, _ = lfq_protint.estimate_protein_intensities(
        protein_df, min_nonan=1, num_samples_quadratic=100, num_cores=1
    )
    assert len(protein_df_normed.index) == 4


def test_nameswitch_indices_are_recovered_correctly():
    array = np.array(["a", "a", "a", "b", "c", "c", "d"])
    start_indices = lfq_protint.find_nameswitch_indices(array)

    assert (array[start_indices[0] : start_indices[1]] == "a").all()
    assert len(array[start_indices[0] : start_indices[1]]) == 3

    assert (array[start_indices[1] : start_indices[2]] == "b").all()
    assert len(array[start_indices[1] : start_indices[2]]) == 1

    assert (array[start_indices[2] : start_indices[3]] == "c").all()
    assert len(array[start_indices[2] : start_indices[3]]) == 2

    assert (array[start_indices[3] : start_indices[4]] == "d").all()
    assert len(array[start_indices[3] : start_indices[4]]) == 1


# ============================================================================
# get_list_with_protein_value_for_each_sample (optimization steps 1 & 5)
# ============================================================================
# Contract: for each sample (column) return the nanmedian of its ion values
# when at least ``min_nonan`` are finite, otherwise NaN. Step 1 vectorized the
# per-sample Series loop; step 5 changed the input from a DataFrame to a numpy
# array (rows = ions, columns = samples).


def _profile_array(rows):
    """Build an ion-profile array (rows = ions, columns = samples)."""
    return np.array(rows, dtype=float)


def test_protein_value_per_sample_returns_column_nanmedians_when_all_finite():
    # given - three samples, all values finite
    arr = _profile_array([[1.0, 4.0, 10.0], [3.0, 6.0, 20.0], [5.0, 8.0, 30.0]])

    # when
    result = lfq_protint.get_list_with_protein_value_for_each_sample(arr, min_nonan=1)

    # then - per-column median, in column order
    assert np.array_equal(np.asarray(result), np.array([3.0, 6.0, 20.0]))


def test_protein_value_per_sample_skips_nans_in_median():
    # given - a column with an even number of finite values (median is averaged)
    arr = _profile_array([[4.0], [np.nan], [6.0]])

    # when
    result = lfq_protint.get_list_with_protein_value_for_each_sample(arr, min_nonan=1)

    # then - nanmedian over {4, 6}
    assert np.array_equal(np.asarray(result), np.array([5.0]))


def test_protein_value_per_sample_all_nan_column_is_nan():
    # given - one all-NaN sample alongside a finite one
    arr = _profile_array([[1.0, np.nan], [3.0, np.nan]])

    # when
    result = lfq_protint.get_list_with_protein_value_for_each_sample(arr, min_nonan=1)

    # then - finite column keeps its median, all-NaN column -> NaN
    assert np.array_equal(np.asarray(result), np.array([2.0, np.nan]), equal_nan=True)


@pytest.mark.parametrize(
    "min_nonan,expected",
    [
        (1, [3.0, 5.0, np.nan]),
        (2, [3.0, 5.0, np.nan]),
        (3, [3.0, np.nan, np.nan]),
    ],
)
def test_protein_value_per_sample_applies_min_nonan_threshold(min_nonan, expected):
    # given - sample finite-counts of 3, 2 and 0 respectively
    arr = _profile_array(
        [[1.0, 4.0, np.nan], [3.0, np.nan, np.nan], [5.0, 6.0, np.nan]]
    )

    # when
    result = lfq_protint.get_list_with_protein_value_for_each_sample(arr, min_nonan)

    # then - a column is NaN-ed out only when its finite count < min_nonan
    assert np.array_equal(np.asarray(result), np.array(expected), equal_nan=True)


# ============================================================================
# Per-protein numpy helpers (optimization step 5)
# ============================================================================
# Step 5 ships numpy slices to workers instead of per-protein DataFrames. The
# normalization helper is validated against the retained NormalizationManagerProtein
# that it replaces on the hot path; the cutting helper is checked against hardcoded
# expectations (the old ProtvalCutter has been removed).
# Gotcha 1 (Fortran-order summed_pepint on >100-ion proteins) is covered by the
# bit-exact reference check (optbench), which the unit data is too small to show.


def _multiindex_df(values, n_ions):
    index = pd.MultiIndex.from_tuples(
        [("P", f"i{j}") for j in range(n_ions)], names=["protein", "ion"]
    )
    return pd.DataFrame(values, index=index)


def test_cut_peptide_values_including_ties():
    # given - 4 ions, two of which tie on both nan-count and summed intensity
    ion_names = np.array(["i0", "i1", "i2", "i3"])
    values = np.array(
        [[2.0, 2.0], [2.0, 2.0], [10.0, 10.0], [np.nan, 5.0]]
    )  # nan-counts 0,0,0,1; sums 4,4,20,5

    # when
    cut_values, cut_names = lfq_protint._cut_peptide_values(
        values, ion_names, maximum=3
    )

    # then - lowest nan-count + highest-sum first; full tie keeps i0 before i1
    assert list(cut_names) == ["i2", "i0", "i1"]
    expected_values = values[[2, 0, 1]]
    assert np.array_equal(cut_values, expected_values, equal_nan=True)


def test_cut_peptide_values_is_noop_within_limit():
    # given - fewer ions than the maximum
    ion_names = np.array(["i0", "i1"])
    values = np.array([[1.0, 2.0], [3.0, 4.0]])

    # when
    cut_values, cut_names = lfq_protint._cut_peptide_values(
        values, ion_names, maximum=100
    )

    # then - returned unchanged (same order, same values)
    assert list(cut_names) == ["i0", "i1"]
    assert np.array_equal(cut_values, values)


def test_normalize_protein_values_matches_manager_quadratic_linear():
    # given - 5 ions > k=3 (quadratic+linear path), with mixed NaNs
    values = np.array(
        [
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 3.0, 4.0, 5.0],
            [10.0, np.nan, 12.0, 13.0],
            [3.0, 4.0, 5.0, 6.0],
            [np.nan, np.nan, 20.0, 21.0],
        ]
    )
    df = _multiindex_df(values, n_ions=5)
    expected = lfq_norm.NormalizationManagerProtein(
        df.copy(), num_samples_quadratic=3
    ).complete_dataframe.to_numpy()

    # when
    result = lfq_protint._normalize_protein_values(
        values.copy(), num_samples_quadratic=3
    )

    # then
    assert np.array_equal(result, expected, equal_nan=True)


def test_normalize_protein_values_matches_manager_quadratic_only():
    # given - 3 ions <= k=5 (quadratic-only path)
    values = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
    df = _multiindex_df(values, n_ions=3)
    expected = lfq_norm.NormalizationManagerProtein(
        df.copy(), num_samples_quadratic=5
    ).complete_dataframe.to_numpy()

    # when
    result = lfq_protint._normalize_protein_values(
        values.copy(), num_samples_quadratic=5
    )

    # then
    assert np.array_equal(result, expected, equal_nan=True)


def test_calculate_peptide_and_protein_intensities_returns_named_numpy_tuple():
    # given - a 2-ion protein over 3 samples
    ion_names = np.array(["i0", "i1"])
    peptide_values = np.array([[1.0, 2.0, 3.0], [2.0, 3.0, 4.0]])

    # when
    profile, name, names_out, shifted = (
        lfq_protint.calculate_peptide_and_protein_intensities(
            0, "protA", ion_names, peptide_values, num_samples_quadratic=10, min_nonan=1
        )
    )

    # then - the 4-tuple carries the protein name, ion names, and shifted values
    assert name == "protA"
    assert list(names_out) == ["i0", "i1"]
    assert np.array_equal(
        shifted,
        lfq_protint._normalize_protein_values(peptide_values, num_samples_quadratic=10),
        equal_nan=True,
    )
    assert profile.shape == (3,)


def test_ion_intensity_dataframe_built_from_tuples():
    # given - per-protein 4-tuples (profile, name, ion_names, shifted_values)
    tuples = [
        (None, "protA", np.array(["a1", "a2"]), np.array([[1.0, 2.0], [np.nan, 3.0]])),
        (None, "protB", np.array(["b1"]), np.array([[0.0, 1.0]])),
    ]

    # when
    ion_df = lfq_protint.get_ion_intensity_dataframe_from_list_of_shifted_peptides(
        tuples, column_names=["s1", "s2"]
    )

    # then - values are 2**shifted with NaN -> 0, indexed by (protein, ion)
    expected = pd.DataFrame(
        np.nan_to_num(2 ** np.array([[1.0, 2.0], [np.nan, 3.0], [0.0, 1.0]])),
        columns=["s1", "s2"],
    )
    expected["ion"] = ["a1", "a2", "b1"]
    expected["protein"] = ["protA", "protA", "protB"]
    expected = expected.set_index(["protein", "ion"])
    pd.testing.assert_frame_equal(ion_df, expected)
