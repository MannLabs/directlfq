"""Tests for directlfq.tracefilter, salvaged from nbdev_nbs/08_tracefilter.ipynb."""

import numpy as np

import directlfq.tracefilter as lfq_trace_filter


def test_empty_matrix():
    lower_matrix = np.array([[]])
    expected = np.array([[]])
    actual = lfq_trace_filter.convert_lower_to_full_matrix(lower_matrix)
    assert np.all(actual == expected), "Failed on empty matrix"


def test_large_matrix_fully_connected():
    lower_matrix = np.tril(np.full((100, 100), np.inf))
    lower_matrix[lower_matrix == 0] = 1

    actual = lfq_trace_filter.get_unconnected_sample_idxs(lower_matrix)
    expected = []
    assert np.all(expected == actual), "Failed on large matrix"


def test_large_matrix():
    lower_matrix = np.full((100, 100), np.inf)

    lower_matrix[1, 0] = 1
    lower_matrix[2, 1] = 1

    actual = lfq_trace_filter.get_unconnected_sample_idxs(lower_matrix)
    expected = list(range(lower_matrix.shape[0]))
    expected = [item for item in expected if item not in [0, 1, 2]]
    assert np.all(expected == actual), "Failed on large matrix"


def test_patterned_matrix():
    lower_matrix = np.array([
        [np.inf, np.inf, np.inf],
        [3, np.inf, np.inf],
        [np.inf, 2, np.inf]
    ])
    expected = np.array([
        [np.inf, 3, np.inf],
        [3, np.inf, 2],
        [np.inf, 2, np.inf]
    ])
    actual = lfq_trace_filter.convert_lower_to_full_matrix(lower_matrix)
    assert np.all(actual == expected), "Failed on patterned matrix"


def test_convert_lower_to_full_matrix():
    lower_matrix = np.array([
        [np.inf, np.inf, np.inf, np.inf],
        [1, np.inf, np.inf, np.inf],
        [np.inf, 1, np.inf, np.inf],
        [np.inf, np.inf, np.inf, np.inf]
    ])
    expected_full_matrix = np.array([
        [np.inf, 1, np.inf, np.inf],
        [1, np.inf, 1, np.inf],
        [np.inf, 1, np.inf, np.inf],
        [np.inf, np.inf, np.inf, np.inf]
    ])
    actual_full_matrix = lfq_trace_filter.convert_lower_to_full_matrix(lower_matrix)
    assert np.all(actual_full_matrix == expected_full_matrix)


def test_get_unconnected_sample_idxs():
    lower_matrix = np.array([
        [np.inf, np.inf, np.inf, np.inf],
        [1, np.inf, np.inf, np.inf],
        [np.inf, 1, np.inf, np.inf],
        [np.inf, np.inf, np.inf, np.inf]
    ])
    expected = [3]
    actual = lfq_trace_filter.get_unconnected_sample_idxs(lower_matrix)
    assert actual == expected


def test_get_unconnected_sample_idxs2():
    lower_matrix = np.array([
        [np.inf, np.inf, np.inf, np.inf],
        [np.inf, np.inf, np.inf, np.inf],
        [np.inf, np.inf, np.inf, np.inf],
        [np.inf, np.inf, np.inf, np.inf]
    ])
    expected = [0, 1, 2, 3]
    actual = lfq_trace_filter.get_unconnected_sample_idxs(lower_matrix)
    assert np.all(actual == expected)


def test_execution_time():
    """Smoke test that exclude_unconnected_samples runs on a large random matrix."""
    lower_matrix = np.tril(np.full((1000, 1000), np.inf))
    for i in range(1, 1000):
        lower_matrix[i, np.random.randint(0, i)] = np.random.random()

    lfq_trace_filter.exclude_unconnected_samples(lower_matrix)


def test_exclusion_of_unconnected_samples():
    lower_matrix = np.array([
        [np.inf, np.inf, np.inf, np.inf, np.inf],
        [1, np.inf, np.inf, np.inf, np.inf],
        [np.inf, 1, np.inf, np.inf, np.inf],
        [np.inf, np.inf, np.inf, np.inf, np.inf],
        [np.inf, np.inf, np.inf, 1, np.inf]
    ])
    lfq_trace_filter.exclude_unconnected_samples(lower_matrix)
    assert (lower_matrix[4, :] == np.array([np.inf, np.inf, np.inf, np.inf, np.inf])).all()
