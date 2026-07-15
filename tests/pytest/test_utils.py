"""Tests for directlfq.utils."""

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

import directlfq.config as config
import directlfq.utils as lfqutils


@pytest.fixture
def _restore_intensity_dtype():
    original = config.INTENSITY_DTYPE
    yield
    config.INTENSITY_DTYPE = original


# ============================================================================
# index_and_log_transform_input_df (memory optimization A1b)
# ============================================================================
# Contract locked in before switching the pandas replace().log2() chain to an
# in-place numpy transform: set a (protein, ion) MultiIndex on the sample
# columns, map zeros to NaN, and return the log2 of every intensity (log2 of a
# NaN stays NaN). Output is float64 regardless of the input column dtype.


def _wide_input_df(sample_columns):
    """Wide input frame: protein/ion id columns plus the given sample columns."""
    data = {
        config.PROTEIN_ID: ["protA", "protA", "protB"],
        config.QUANT_ID: ["ion0", "ion1", "ion2"],
    }
    data.update(sample_columns)
    return pd.DataFrame(data)


def test_index_and_log_transform_sets_multiindex_and_log2s_values():
    # given - zeros in S1, a NaN in S2, otherwise plain intensities
    df = _wide_input_df({"S1": [1.0, 0.0, 2.0], "S2": [4.0, 8.0, np.nan]})

    # when
    result = lfqutils.index_and_log_transform_input_df(df)

    # then - (protein, ion) index, sample columns log2-ed, 0 -> NaN, NaN kept
    expected = pd.DataFrame(
        {"S1": [0.0, np.nan, 1.0], "S2": [2.0, 3.0, np.nan]},
        index=pd.MultiIndex.from_arrays(
            [["protA", "protA", "protB"], ["ion0", "ion1", "ion2"]],
            names=[config.PROTEIN_ID, config.QUANT_ID],
        ),
    )
    assert_frame_equal(result, expected)


def test_index_and_log_transform_returns_float64_for_integer_input():
    # given - integer intensity columns
    df = _wide_input_df({"S1": [1, 2, 4], "S2": [8, 16, 32]})

    # when
    result = lfqutils.index_and_log_transform_input_df(df)

    # then - integer input is transformed to float64 log2 values
    expected = pd.DataFrame(
        {"S1": [0.0, 1.0, 2.0], "S2": [3.0, 4.0, 5.0]},
        index=pd.MultiIndex.from_arrays(
            [["protA", "protA", "protB"], ["ion0", "ion1", "ion2"]],
            names=[config.PROTEIN_ID, config.QUANT_ID],
        ),
    )
    assert_frame_equal(result, expected)
    assert all(result.dtypes == np.float64)


def test_index_and_log_transform_does_not_mutate_input_df():
    # given
    df = _wide_input_df({"S1": [1.0, 0.0, 2.0]})
    df_before = df.copy()

    # when
    lfqutils.index_and_log_transform_input_df(df)

    # then - the caller's frame is untouched
    assert_frame_equal(df, df_before)


def test_index_and_log_transform_follows_config_intensity_dtype(
    _restore_intensity_dtype,
):  # noqa: ARG001
    # given - the opt-in float32 dtype is selected
    df = _wide_input_df({"S1": [1.0, 4.0, 2.0], "S2": [8.0, 16.0, 32.0]})
    config.INTENSITY_DTYPE = np.float32

    # when
    result = lfqutils.index_and_log_transform_input_df(df)

    # then - output is float32 and stays close to the float64 log2 values
    assert all(result.dtypes == np.float32)
    config.INTENSITY_DTYPE = np.float64
    expected_float64 = lfqutils.index_and_log_transform_input_df(df)
    assert np.allclose(
        result.to_numpy(), expected_float64.to_numpy(), rtol=1e-6, equal_nan=True
    )
