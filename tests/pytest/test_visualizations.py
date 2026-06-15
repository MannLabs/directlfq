"""Tests for directlfq.visualizations, salvaged from nbdev_nbs/05_visualizations.ipynb.

These assert on matplotlib line counts and are inherently coupled to the plotting
implementation; treat failures as a prompt to check whether the plotting changed.
"""

import os

import matplotlib.pyplot as plt
import pandas as pd

import directlfq.utils as lfq_utils
import directlfq.visualizations as lfq_viz


def test_that_iontracevisualized_produces_desired_plot(test_data_dir):
    example_prots = os.path.join(
        test_data_dir, "unit_tests", "protein_normalization", "example_proteins.tsv"
    )
    protein_df = pd.read_csv(example_prots, sep="\t")
    protein_df = lfq_utils.index_and_log_transform_input_df(protein_df)
    protein_df = protein_df.loc["A0A024R4E5"]

    ax = plt.subplot()
    pviz = lfq_viz.IonTraceVisualizer(protein_df, ax)
    pviz.add_median_trace([17 for _ in range(len(protein_df.columns))])

    assert len(ax.lines) == len(protein_df.index) + 1


def test_that_iontracecomparisonplotter_produces_desired_plots(test_data_dir):
    example_prots = os.path.join(
        test_data_dir, "unit_tests", "protein_normalization", "example_proteins.tsv"
    )
    _, axes = plt.subplots(1, 2)
    protein_df = pd.read_csv(example_prots, sep="\t")
    protein_df = lfq_utils.index_and_log_transform_input_df(protein_df)
    complotter = lfq_viz.IonTraceCompararisonPlotter(
        protein_df,
        selected_protein="A0A024R4E5",
        axis_normed=axes[0],
        axis_unnormed=axes[1],
    )
    assert len(complotter.axis_unnormed.lines) == len(
        protein_df.loc["A0A024R4E5"].index
    )
    assert (
        len(complotter.axis_normed.lines) == len(protein_df.loc["A0A024R4E5"].index) + 1
    )
