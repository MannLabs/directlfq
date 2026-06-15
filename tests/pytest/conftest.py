"""Shared fixtures for the directlfq pytest suite."""

import os

import matplotlib
import pytest

# Use a non-interactive backend so visualization tests run headless in CI.
matplotlib.use("Agg")


@pytest.fixture
def test_data_dir():
    """Absolute path to the repository's ``test_data`` directory."""
    here = os.path.dirname(os.path.realpath(__file__))
    return os.path.abspath(os.path.join(here, "..", "..", "test_data"))
