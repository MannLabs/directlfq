"""Tests for directlfq.protein_intensity_estimation, salvaged from nbdev_nbs/03_protein_intensity_estimation.ipynb."""

import numpy as np
import pandas as pd

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
