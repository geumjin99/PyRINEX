"""PyRINEX — Python toolkit for RINEX 2/3 batch processing and quality check.

Public API:

* :func:`reader.read_obs_header`, :func:`reader.read_obs`,
  :func:`reader.read_nav` — RINEX parsers returning plain Python dicts.
* :func:`data_management.find_rinex`, :func:`data_management.clean_rinex`
  — bulk discovery and rewrite helpers.
* :func:`quality_check.quality_check`,
  :func:`quality_check.batch_quality_check` — full QC pipeline driving
  multipath, ionosphere, cycle-slip, sky-plot and signal-plot products.

Legacy ``oheader`` / ``observations`` / ``navigations`` / ``QualityCheck``
names continue to work and emit ``DeprecationWarning``. The misspelled
:mod:`PyRINEX.QulityCheck` module re-exports the same symbols for callers
that were imported by following the v3.x manual.
"""
from __future__ import annotations

__version__ = "4.0.0"

from . import coordinates, data_management, frequencies, orbit, plotting, quality_check, reader, time_utils
from .data_management import (
    DataCleaning,
    DataFinding,
    clean_rinex,
    find_rinex,
)
from .quality_check import (
    QualityCheck,
    azimuth_elevation,
    batch_quality_check,
    batchQC,
    cycle_slip,
    multipath_iod,
    quality_check as run_quality_check,
    satellite_signal_plot,
)
from .reader import (
    ObsHeader,
    navigations,
    observations,
    oheader,
    read_nav,
    read_nav_iono,
    read_obs,
    read_obs_header,
)

__all__ = [
    "__version__",
    "ObsHeader",
    "DataCleaning",
    "DataFinding",
    "QualityCheck",
    "azimuth_elevation",
    "batch_quality_check",
    "batchQC",
    "clean_rinex",
    "coordinates",
    "cycle_slip",
    "data_management",
    "find_rinex",
    "frequencies",
    "multipath_iod",
    "navigations",
    "observations",
    "oheader",
    "orbit",
    "plotting",
    "quality_check",
    "read_nav",
    "read_nav_iono",
    "read_obs",
    "read_obs_header",
    "reader",
    "run_quality_check",
    "satellite_signal_plot",
    "time_utils",
]
