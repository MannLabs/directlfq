"""Tests for directlfq.protein_intensity_estimation, salvaged from nbdev_nbs/03_protein_intensity_estimation.ipynb."""

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

import directlfq.protein_intensity_estimation as lfq_protint
import directlfq.test_utils as lfq_testutils


# ============================================================================
# _cut_peptide_values (optimization step 5)
# ============================================================================
# The numpy cutter keeps the <= maximum ions sorted by NaN count asc then summed
# intensity desc, reproducing the former ProtvalCutter's stable sorted() tie-break
# (ions tied on both keys keep their original order). Expectations are hardcoded.


def test_cut_peptide_values_orders_by_nan_then_intensity_with_full_tie():
    # given - ion0 and ion1 are a full tie on both keys (nan=0, sum=6); the cut
    # keeps the top 3, so the tie-break (original order) is observable
    ion_names = np.array(["ion0", "ion1", "ion2", "ion3", "ion4"])
    peptide_values = np.array(
        [
            [1.0, 2.0, 3.0],  # ion0: nan=0 sum=6
            [1.0, 2.0, 3.0],  # ion1: nan=0 sum=6 (full tie with ion0)
            [10.0, 20.0, 30.0],  # ion2: nan=0 sum=60
            [5.0, np.nan, np.nan],  # ion3: nan=2 sum=5
            [np.nan, np.nan, np.nan],  # ion4: nan=3 sum=0
        ]
    )

    # when
    new_vals, new_ions = lfq_protint._cut_peptide_values(
        peptide_values, ion_names, maximum=3
    )

    # then - highest-sum first, then the full tie in original (ion-name) order
    assert list(new_ions) == ["ion2", "ion0", "ion1"]
    assert np.array_equal(
        new_vals, np.array([[10.0, 20.0, 30.0], [1.0, 2.0, 3.0], [1.0, 2.0, 3.0]])
    )


def test_cut_peptide_values_is_noop_within_limit():
    ion_names = np.array(["ion0", "ion1", "ion2"])
    peptide_values = np.array([[1.0, 2.0], [3.0, np.nan], [5.0, 6.0]])

    new_vals, new_ions = lfq_protint._cut_peptide_values(
        peptide_values, ion_names, maximum=100
    )

    assert list(new_ions) == ["ion0", "ion1", "ion2"]
    assert np.array_equal(new_vals, peptide_values, equal_nan=True)


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


def _profile_array(rows):
    """Build an ion-profile array (rows = ions, columns = samples)."""
    return np.array(rows, dtype=float)


def test_protein_value_per_sample_returns_column_nanmedians_when_all_finite():
    # given - three samples, all values finite
    df = _profile_array([[1.0, 4.0, 10.0], [3.0, 6.0, 20.0], [5.0, 8.0, 30.0]])

    # when
    result = lfq_protint.get_list_with_protein_value_for_each_sample(df, min_nonan=1)

    # then - per-column median, in column order
    assert np.array_equal(np.asarray(result), np.array([3.0, 6.0, 20.0]))


def test_protein_value_per_sample_skips_nans_in_median():
    # given - a column with an even number of finite values (median is averaged)
    df = _profile_array([[4.0], [np.nan], [6.0]])

    # when
    result = lfq_protint.get_list_with_protein_value_for_each_sample(df, min_nonan=1)

    # then - nanmedian over {4, 6}
    assert np.array_equal(np.asarray(result), np.array([5.0]))


def test_protein_value_per_sample_all_nan_column_is_nan():
    # given - one all-NaN sample alongside a finite one
    df = _profile_array([[1.0, np.nan], [3.0, np.nan]])

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
    df = _profile_array([[1.0, 4.0, np.nan], [3.0, np.nan, np.nan], [5.0, 6.0, np.nan]])

    # when
    result = lfq_protint.get_list_with_protein_value_for_each_sample(df, min_nonan)

    # then - a column is NaN-ed out only when its finite count < min_nonan
    assert np.array_equal(np.asarray(result), np.array(expected), equal_nan=True)


# ============================================================================
# calculate_peptide_and_protein_intensities single-sample guard (step 5)
# ============================================================================
# With a single sample every ion has one intensity, so get_normfacts would NaN
# every row and drop the protein. The guard skips normalization to keep values.


def test_single_sample_values_are_kept_and_protein_retained():
    # given - 3 ions, 1 sample
    peptide_values = np.array([[5.0], [6.0], [7.0]])
    ion_names = np.array(["ion0", "ion1", "ion2"])

    # when
    profile, name, ions, shifted = (
        lfq_protint.calculate_peptide_and_protein_intensities(
            0, "protA", ion_names, peptide_values, 10, 1
        )
    )

    # then - values kept as-is (not normalized to NaN) and the protein is retained
    assert np.array_equal(shifted, peptide_values, equal_nan=True)
    assert profile is not None
    assert name == "protA"
    assert list(ions) == ["ion0", "ion1", "ion2"]


# ============================================================================
# get_ion_intensity_dataframe_from_list_of_shifted_peptides (memory optim. A2)
# ============================================================================
# Contract locked in before switching the ion-table nan_to_num to in-place: each
# per-protein (profile, name, ion_names, shifted_values) tuple contributes its
# ions to one (protein, ion)-indexed frame; shifted log2 values are raised to
# 2**x back into linear space and NaNs are replaced with 0.


def test_ion_intensity_dataframe_delogs_values_and_zeroes_nans():
    # given - two proteins; shifted values are in log2 space, one NaN per protein
    tuples = [
        (
            None,
            "protA",
            np.array(["ion0", "ion1"]),
            np.array([[1.0, 2.0], [3.0, np.nan]]),
        ),
        (None, "protB", np.array(["ion2"]), np.array([[0.0, np.nan]])),
    ]

    # when
    result = lfq_protint.get_ion_intensity_dataframe_from_list_of_shifted_peptides(
        tuples, column_names=["S1", "S2"]
    )

    # then - 2**x back to linear space, NaN -> 0, indexed by (protein, ion)
    expected = pd.DataFrame(
        {"S1": [2.0, 8.0, 1.0], "S2": [4.0, 0.0, 0.0]},
        index=pd.MultiIndex.from_arrays(
            [["protA", "protA", "protB"], ["ion0", "ion1", "ion2"]],
            names=["protein", "ion"],
        ),
    )
    assert_frame_equal(result, expected)
