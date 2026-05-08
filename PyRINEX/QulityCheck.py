"""Backwards-compatibility shim for the misspelled v3.x module name.

The original module name had a typo (``Qulity`` instead of ``Quality``).
The published manual and any user code that followed it imports from
``PyRINEX.QulityCheck``, so this module re-exports every public symbol
from :mod:`PyRINEX.quality_check` and emits a :class:`DeprecationWarning`
on import.
"""
from __future__ import annotations

import warnings

warnings.warn(
    "PyRINEX.QulityCheck is a deprecated alias for PyRINEX.quality_check; "
    "import from PyRINEX.quality_check (correct spelling) instead.",
    DeprecationWarning,
    stacklevel=2,
)

from .quality_check import (  # noqa: E402, F401
    ION_MP,
    LandC,
    QualityCheck,
    SatelliteSignalPlot,
    azi_ele,
    azimuth_elevation,
    batch_quality_check,
    batchQC,
    cycle_slip,
    cycleslip,
    extract_l_c,
    frelist,
    gps_pos,
    gps_positions,
    multipath_iod,
    plot,
    quality_check,
    satellite_signal_plot,
    time_serises,
)
