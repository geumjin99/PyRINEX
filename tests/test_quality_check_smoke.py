"""Smoke tests for the quality-check pipeline.

These exercise the data-flow end-to-end on the synthetic fixtures. They
do not assert numerical correctness against published reference values
(the synthetic fixtures aren't physically meaningful enough for that)
but they catch shape errors, indexing bugs, and wrong return types.
"""
from __future__ import annotations

import os

import matplotlib
import pytest

matplotlib.use("Agg")  # headless plotting

from PyRINEX import quality_check


def test_extract_l_c_shape(v2_obs_path):
    array = quality_check.extract_l_c(v2_obs_path)
    assert array.ndim == 3
    assert array.shape[2] == 4  # L1, L2, P1/C1, P2/C2
    # G01..G03 are present in every epoch with non-zero L1.
    assert (array[:, :, 0] != 0).any()


def test_frequency_table_matches_prn_count(v2_obs_path):
    from PyRINEX.reader import read_obs_header

    header = read_obs_header(v2_obs_path)
    table = quality_check._frequency_table(header.PRNS)
    assert table.shape == (len(header.PRNS), 4)


def test_satellite_signal_plot_writes_image(v2_obs_path, tmp_path):
    quality_check.satellite_signal_plot(v2_obs_path)
    image = str(v2_obs_path)[:-4] + "SignalPlot.jpg"
    assert os.path.exists(image)


def test_multipath_iod_writes_outputs(v2_obs_path):
    result = quality_check.multipath_iod(v2_obs_path)
    assert result.shape[2] == 4
    base = str(v2_obs_path)[:-4]
    assert os.path.exists(base + "mp1mp2.csv")
    assert os.path.exists(base + "IonIod.csv")
    assert os.path.exists(base + "MP1_plot.png")


def test_quality_check_runs_end_to_end(v2_pair):
    obs_path, _nav_path = v2_pair
    quality_check.quality_check(obs_path)
    base = str(obs_path)[:-4]
    assert os.path.exists(base + "SignalPlot.jpg")
    assert os.path.exists(base + "mp1mp2.csv")


def test_legacy_qulitycheck_module_emits_warning(v2_pair):
    import warnings

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        import importlib

        if "PyRINEX.QulityCheck" in __import__("sys").modules:
            del __import__("sys").modules["PyRINEX.QulityCheck"]
        importlib.import_module("PyRINEX.QulityCheck")
        assert any(issubclass(rec.category, DeprecationWarning) for rec in w)
