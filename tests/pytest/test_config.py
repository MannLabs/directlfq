"""Tests for directlfq.config."""

import numpy as np
import pytest

import directlfq.config as config


# ============================================================================
# set_intensity_dtype (memory optimization A3)
# ============================================================================
# The opt-in float32 path flips a single module-global that every downstream
# intensity matrix inherits. Default must stay float64 to keep the bit-exact path.


@pytest.fixture(autouse=True)
def _restore_intensity_dtype():
    original = config.INTENSITY_DTYPE
    yield
    config.INTENSITY_DTYPE = original


def test_set_intensity_dtype_defaults_to_float64():
    config.set_intensity_dtype()

    assert config.INTENSITY_DTYPE == np.float64


def test_set_intensity_dtype_true_selects_float32():
    config.set_intensity_dtype(use_float32=True)

    assert config.INTENSITY_DTYPE == np.float32


def test_set_intensity_dtype_false_selects_float64():
    config.set_intensity_dtype(use_float32=True)

    config.set_intensity_dtype(use_float32=False)

    assert config.INTENSITY_DTYPE == np.float64
