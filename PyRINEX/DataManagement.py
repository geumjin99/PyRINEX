"""Backwards-compatibility shim for the v3.x ``DataManagement`` module
name. New code should import from :mod:`PyRINEX.data_management`."""
from __future__ import annotations

import warnings

warnings.warn(
    "PyRINEX.DataManagement is a deprecated alias for PyRINEX.data_management; "
    "import from PyRINEX.data_management (PEP-8 name) instead.",
    DeprecationWarning,
    stacklevel=2,
)

from .data_management import (  # noqa: E402, F401
    DataCleaning,
    DataFinding,
    clean_rinex,
    date2doy,
    date_to_doy,
    detectCode,
    file_name_extesion_judgement,
    file_name_extesion_judgement_file,
    filter_by_extension,
    find_rinex,
    generate_unique_filename,
    readtxt,
    show_files,
    unique_filename,
    xyz2BLH,
)
