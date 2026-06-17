"""Tests for directlfq.protein_intensity_estimation, salvaged from nbdev_nbs/03_protein_intensity_estimation.ipynb."""

import numpy as np
import pandas as pd
import pytest

import directlfq.protein_intensity_estimation as lfq_protint
import directlfq.test_utils as lfq_testutils


def test_sorting_by_num_nans():
    vals1 = np.array([9, np.nan, np.nan, np.nan])
    vals2 = np.array([5, 6, np.nan, np.nan])
    vals3 = np.array([1, 2, 3, np.nan])

    df = pd.DataFrame([vals1, vals2, vals3], index=[["P", "P", "P"], ["A", "B", "C"]])
    pcutter = lfq_protint.ProtvalCutter(df, maximum_df_length=2)
    sorted_idx = pcutter._sorted_idx
    df_sorted = df.loc[sorted_idx]

    assert np.allclose(df_sorted.iloc[2].to_numpy(), vals1, equal_nan=True)
    assert np.allclose(df_sorted.iloc[0].to_numpy(), vals3, equal_nan=True)


def test_cutting_of_df():
    vals1 = np.array([9, np.nan, np.nan, np.nan])
    vals2 = np.array([5, 6, np.nan, np.nan])
    vals3 = np.array([1, 2, 3, np.nan])

    df = pd.DataFrame([vals1, vals2, vals3], index=[["A", "B", "C"]])
    pcutter = lfq_protint.ProtvalCutter(df, maximum_df_length=2)
    cut_df = pcutter.get_dataframe()
    ion_idx = [x[0] for x in cut_df.index]
    assert ion_idx == ["C", "B"]


# ============================================================================
# _cut_peptide_values (optimization step 5) -- proven against live ProtvalCutter
# ============================================================================
# Strangler check: the numpy cutter must produce the SAME ion order (not just the
# same set) as ProtvalCutter, including the stable tie-break on full ties.


def _make_ion_df(rows, ion_names):
    index = pd.MultiIndex.from_arrays([["P"] * len(ion_names), ion_names])
    return pd.DataFrame(rows, index=index)


def test_cut_peptide_values_matches_protvalcutter_order_with_full_tie():
    # given - ion0 and ion1 are a full tie on both keys (nan=0, sum=6); the cut
    # keeps the top 3, so the tie-break (original order) is observable
    ion_names = ["ion0", "ion1", "ion2", "ion3", "ion4"]
    rows = [
        [1.0, 2.0, 3.0],  # ion0: nan=0 sum=6
        [1.0, 2.0, 3.0],  # ion1: nan=0 sum=6 (full tie with ion0)
        [10.0, 20.0, 30.0],  # ion2: nan=0 sum=60
        [5.0, np.nan, np.nan],  # ion3: nan=2 sum=5
        [np.nan, np.nan, np.nan],  # ion4: nan=3 sum=0
    ]
    df = _make_ion_df(rows, ion_names)

    cut_df = lfq_protint.ProtvalCutter(
        df.copy(), maximum_df_length=3
    ).get_dataframe()
    pc_ions = [t[1] for t in cut_df.index]

    # when
    new_vals, new_ions = lfq_protint._cut_peptide_values(
        df.to_numpy(), np.array(ion_names), maximum=3
    )

    # then - identical ion order and values
    assert list(new_ions) == pc_ions
    assert np.array_equal(new_vals, cut_df.to_numpy(), equal_nan=True)


def test_cut_peptide_values_is_noop_within_limit():
    ion_names = ["ion0", "ion1", "ion2"]
    rows = [[1.0, 2.0], [3.0, np.nan], [5.0, 6.0]]
    df = _make_ion_df(rows, ion_names)

    cut_df = lfq_protint.ProtvalCutter(
        df.copy(), maximum_df_length=100
    ).get_dataframe()
    pc_ions = [t[1] for t in cut_df.index]

    new_vals, new_ions = lfq_protint._cut_peptide_values(
        df.to_numpy(), np.array(ion_names), maximum=100
    )

    assert list(new_ions) == pc_ions == ion_names
    assert np.array_equal(new_vals, cut_df.to_numpy(), equal_nan=True)


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
# get_list_with_protein_value_for_each_sample (optimization step 1)
# ============================================================================
# Contract locked in before vectorizing the per-sample Series loop into a single
# numpy column-nanmedian pass: for each sample (column) return the nanmedian of
# its ion values when at least ``min_nonan`` are finite, otherwise NaN.


def _profile_df(rows):
    """Build an ion-profile DataFrame (rows = ions, columns = samples)."""
    return pd.DataFrame(rows)


def test_protein_value_per_sample_returns_column_nanmedians_when_all_finite():
    # given - three samples, all values finite
    df = _profile_df([[1.0, 4.0, 10.0], [3.0, 6.0, 20.0], [5.0, 8.0, 30.0]])

    # when
    result = lfq_protint.get_list_with_protein_value_for_each_sample(df, min_nonan=1)

    # then - per-column median, in column order
    assert np.array_equal(np.asarray(result), np.array([3.0, 6.0, 20.0]))


def test_protein_value_per_sample_skips_nans_in_median():
    # given - a column with an even number of finite values (median is averaged)
    df = _profile_df([[4.0], [np.nan], [6.0]])

    # when
    result = lfq_protint.get_list_with_protein_value_for_each_sample(df, min_nonan=1)

    # then - nanmedian over {4, 6}
    assert np.array_equal(np.asarray(result), np.array([5.0]))


def test_protein_value_per_sample_all_nan_column_is_nan():
    # given - one all-NaN sample alongside a finite one
    df = _profile_df([[1.0, np.nan], [3.0, np.nan]])

    # when
    result = lfq_protint.get_list_with_protein_value_for_each_sample(df, min_nonan=1)

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
    df = _profile_df([[1.0, 4.0, np.nan], [3.0, np.nan, np.nan], [5.0, 6.0, np.nan]])

    # when
    result = lfq_protint.get_list_with_protein_value_for_each_sample(df, min_nonan)

    # then - a column is NaN-ed out only when its finite count < min_nonan
    assert np.array_equal(np.asarray(result), np.array(expected), equal_nan=True)
