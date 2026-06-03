"""Tests for directlfq.protein_intensity_estimation.

These tests exercise ``get_normed_dfs`` only. ``find_nameswitch_indices`` and
``get_subdf`` are covered implicitly, since ``get_normed_dfs`` is their sole caller.
"""

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

import directlfq.config as config
from directlfq.protein_intensity_estimation import get_normed_dfs


# ============================================================================
# HELPERS
# ============================================================================


def _make_normed_df(index_tuples, values, columns=("s1", "s2")):
    """Build a normed_df with the [PROTEIN_ID, QUANT_ID] multi-index used as input."""
    index = pd.MultiIndex.from_tuples(
        index_tuples, names=[config.PROTEIN_ID, config.QUANT_ID]
    )
    return pd.DataFrame(values, index=index, columns=list(columns))


def _make_expected_subdf(index_tuples, values):
    """Build an expected sub-dataframe as produced by get_subdf.

    get_subdf creates the dataframe without passing column labels, so the
    columns are the default integer RangeIndex regardless of the input columns.
    """
    index = pd.MultiIndex.from_tuples(
        index_tuples, names=[config.PROTEIN_ID, config.QUANT_ID]
    )
    return pd.DataFrame(values, index=index)


# ============================================================================
# HAPPY PATH TESTS
# ============================================================================


def test_get_normed_dfs_splits_by_protein_with_mixed_block_sizes() -> None:
    """Test that get_normed_dfs splits a sorted multi-protein frame into one
    sub-dataframe per protein, covering single-ion and multi-ion blocks at the
    first, middle and last positions."""
    # given - 3 proteins with 2 / 1 / 3 ions, 2 samples, one NaN to check passthrough
    normed_df = _make_normed_df(
        index_tuples=[
            ("protA", "ionA1"),
            ("protA", "ionA2"),
            ("protB", "ionB1"),
            ("protC", "ionC1"),
            ("protC", "ionC2"),
            ("protC", "ionC3"),
        ],
        values=[
            [1.0, 2.0],
            [3.0, 4.0],
            [5.0, 6.0],
            [7.0, 8.0],
            [9.0, 10.0],
            [11.0, np.nan],
        ],
    )

    # when
    result = get_normed_dfs(normed_df)

    # then
    expected = [
        _make_expected_subdf(
            [("protA", "ionA1"), ("protA", "ionA2")],
            [[1.0, 2.0], [3.0, 4.0]],
        ),
        _make_expected_subdf(
            [("protB", "ionB1")],
            [[5.0, 6.0]],
        ),
        _make_expected_subdf(
            [("protC", "ionC1"), ("protC", "ionC2"), ("protC", "ionC3")],
            [[7.0, 8.0], [9.0, 10.0], [11.0, np.nan]],
        ),
    ]
    assert len(result) == len(expected)
    for result_df, expected_df in zip(result, expected):
        assert_frame_equal(result_df, expected_df)