"""Tests for directlfq.lfq_manager."""

from unittest.mock import patch

import pytest

import directlfq.lfq_manager as lfq_manager


# ============================================================================
# run_lfq use_float32 wiring (memory optimization A3)
# ============================================================================
# run_lfq must forward its use_float32 flag to config.set_intensity_dtype so the
# opt-in path takes effect. import_data is stubbed to stop the pipeline right
# after the config block, keeping the test independent of any input file.


class _StopPipeline(Exception):
    """Sentinel raised to abort run_lfq once the config block has executed."""


@patch("directlfq.lfq_manager.lfqutils.import_data", side_effect=_StopPipeline)
@patch("directlfq.lfq_manager.config.set_intensity_dtype")
def test_run_lfq_forwards_use_float32_true(mock_set_dtype, mock_import):  # noqa: ARG001
    with pytest.raises(_StopPipeline):
        lfq_manager.run_lfq("dummy_input.tsv", use_float32=True)

    mock_set_dtype.assert_called_once_with(use_float32=True)


@patch("directlfq.lfq_manager.lfqutils.import_data", side_effect=_StopPipeline)
@patch("directlfq.lfq_manager.config.set_intensity_dtype")
def test_run_lfq_defaults_use_float32_to_false(mock_set_dtype, mock_import):  # noqa: ARG001
    with pytest.raises(_StopPipeline):
        lfq_manager.run_lfq("dummy_input.tsv")

    mock_set_dtype.assert_called_once_with(use_float32=False)
