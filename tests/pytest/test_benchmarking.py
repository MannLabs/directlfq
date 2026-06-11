"""Tests for directlfq.benchmarking, salvaged from nbdev_nbs/06_benchmarking.ipynb."""

import os

import numpy as np
import pandas as pd
import pytest

import directlfq.benchmarking as lfq_benchmarking


@pytest.fixture
def quant_df():
    return pd.DataFrame({
        "sample_list": [1, 1, 1, 2, 2, 2, 3],
        "asd": ["a", "b", "c", "d", "e", "f", "g"],
        "cfs": [11, 23, 4, 5, 7, 4, 9],
    })


@pytest.fixture
def samplelist_df():
    return pd.DataFrame({"sample_list": ["s1", "s2", "s3"]})


@pytest.fixture
def template_df(test_data_dir):
    template_df_location = os.path.join(
        test_data_dir, "unit_tests", "protein_normalization", "example_proteins.tsv"
    )
    return pd.read_csv(template_df_location, index_col=["protein", "ion"], sep="\t")


@pytest.mark.parametrize("desired_num_samples", [1, 7, 13])
def test_that_scaled_numbers_of_samples_are_as_expected(quant_df, samplelist_df, desired_num_samples):
    scaled_df_creator = lfq_benchmarking.ScaledDFCreatorIQFormat(quant_df, samplelist_df, desired_num_samples)
    assert len(set(scaled_df_creator.scaled_quant_df["sample_list"])) == desired_num_samples
    assert len(scaled_df_creator.scaled_sample_list_df.index) == desired_num_samples


def test_that_repetition_worked_out(quant_df, samplelist_df):
    scaled_df_creator = lfq_benchmarking.ScaledDFCreatorIQFormat(quant_df, samplelist_df, 7)
    assert np.all(scaled_df_creator.scaled_quant_df["sample_list"][6] == [3])


@pytest.mark.parametrize("num_samples", [1, 100, 10000])
def test_that_shape_is_as_expected(template_df, num_samples):
    size_adjusted_df = lfq_benchmarking.ScaledDFCreatorDirectLFQFormat(
        template_df=template_df, desired_number_of_samples=num_samples
    ).scaled_df
    assert len(size_adjusted_df.columns) == num_samples
    assert len(size_adjusted_df.index) == len(template_df.index)


@pytest.mark.parametrize("num_samples", [5, 100])
def test_that_values_are_as_expected(template_df, num_samples):
    size_adjusted_df = lfq_benchmarking.ScaledDFCreatorDirectLFQFormat(
        template_df=template_df, desired_number_of_samples=num_samples
    ).scaled_df
    assert np.allclose(template_df.loc[:, "BoxCar_02-01_2"], size_adjusted_df.loc[:, "BoxCar_02-01_2_AND_remainder"])
