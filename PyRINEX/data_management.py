"""Bulk-rename / coordinate-translate / non-ASCII-strip helpers for RINEX
observation and navigation files.

Compared to the v3.x ``DataManagement.py`` this module:

* uses :func:`os.path.join` instead of hard-coded ``"\\\\"`` (so it actually
  runs on Linux/macOS),
* fixes the indentation bug that made the extension filters return after
  the first directory level only,
* fixes the no-op ``new_marker.replace(" ", "a")`` (which never did
  anything because strings are immutable),
* fixes the ``TIME OF FIRST OBS / TIME OF LAST OBS`` key collision in the
  reader,
* drops the unused ``codecs`` import,
* exposes WGS-84 conversions through the dedicated :mod:`coordinates`
  module instead of duplicating the formulae.
"""
from __future__ import annotations

import csv
import logging
import os
import re
from pathlib import Path
from typing import Iterable, List, Optional, Sequence

from .coordinates import xyz_to_blh
from .reader import (
    PathLike,
    _detect_encoding,
    list_files,
    read_lines,
    read_obs_header,
)
from .time_utils import date_to_doy

__all__ = [
    "find_rinex",
    "filter_by_extension",
    "clean_rinex",
    "unique_filename",
    "DataFinding",
    "DataCleaning",
    "show_files",
    "readtxt",
    "detectCode",
    "date2doy",
    "xyz2BLH",
]

log = logging.getLogger(__name__)

CSV_FIELDNAMES = (
    "origin path",
    "version",
    "non English",
    "origin marker",
    "origin rec",
    "origin ant",
    "new path",
    "marker",
    "longitude",
    "latitude",
    "rec type",
    "ant type",
)


# ---------------------------------------------------------------------------
# File system helpers
# ---------------------------------------------------------------------------

def filter_by_extension(file_dir: PathLike, extension: str) -> List[str]:
    """Walk ``file_dir`` recursively and return absolute paths to every file
    whose extension matches ``extension`` (case-insensitive)."""
    ext_l = extension.lower()
    matches: List[str] = []
    for root, _dirs, files in os.walk(file_dir):
        for f in files:
            if os.path.splitext(f)[1].lower() == ext_l:
                matches.append(os.path.join(root, f))
    return matches


def find_rinex(
    root_path: PathLike,
    keywords: Iterable[str],
    extension: str,
) -> List[str]:
    """Return every file under ``root_path`` whose folder path contains all
    ``keywords`` and whose name ends with ``extension`` (case-insensitive).

    This is the corrected version of the ``DataFinding`` function from
    v3.x — the same intent, but ``os.path`` instead of string mangling.
    """
    keyword_list = list(keywords)
    ext_l = extension.lower()
    matches: List[str] = []
    for folder, _dirs, files in os.walk(root_path):
        if not all(kw in folder for kw in keyword_list):
            continue
        for f in files:
            if f.lower().endswith(ext_l):
                matches.append(os.path.join(folder, f))
    return matches


def unique_filename(path: PathLike) -> str:
    """If ``path`` exists, append ``_1`` / ``_2`` / … until a free name is
    found; otherwise return ``path`` unchanged."""
    path_str = str(path)
    if not os.path.exists(path_str):
        return path_str
    stem, ext = os.path.splitext(path_str)
    counter = 1
    while os.path.exists(path_str):
        path_str = f"{stem}_{counter}{ext}"
        counter += 1
    return path_str


# ---------------------------------------------------------------------------
# Library lookup tables
# ---------------------------------------------------------------------------

def _load_library(path: PathLike) -> dict[str, str]:
    """Parse a ReceiverLibrary / AntennaLibrary CSV.

    File format is one entry per line: ``<wrong>?<right>``. Commas inside
    either column are stripped to keep the resulting string CSV-safe.
    """
    out: dict[str, str] = {}
    for line in read_lines(path):
        if "?" not in line:
            continue
        wrong, right = line.split("?", 1)
        out[wrong.replace(",", "")] = right.rstrip("\n").replace(",", "")
    return out


# ---------------------------------------------------------------------------
# Cleaning workflow
# ---------------------------------------------------------------------------

def _strip_non_english(line: str) -> tuple[str, bool]:
    """Replace non-Latin / non-Hebrew characters with ``-`` (matches the
    legacy behaviour of stripping CJK characters from corrupt headers)."""
    chars = list(line)
    cleaned = [c if not re.findall(r"[^\u0000-\u05C0\u2100-\u214F]+", c) else "-" for c in chars]
    return "".join(cleaned), cleaned != chars


def clean_rinex(
    rinex_files: Sequence[PathLike],
    receiver_library_path: PathLike,
    antenna_library_path: PathLike,
    output_root: PathLike,
    radome: bool = False,
) -> None:
    """Rewrite RINEX observation files into a marker/DOY-organised tree.

    For each ``.??o`` file in ``rinex_files`` this function:

    1. parses the header,
    2. shortens the marker name to four characters (right-padded with ``a``
       if needed),
    3. looks up the receiver and antenna types in the library CSVs and
       rewrites the corresponding header lines,
    4. optionally fills in ``NONE`` as the antenna radome,
    5. strips non-English characters from the first 20 header lines,
    6. writes the result to ``<output_root>/<yy><doy>/<MARKER><doy>0.<yy>o``,
    7. copies the matching navigation file alongside it,
    8. appends a row to ``<output_root>/report.csv`` summarising the
       transformation.
    """
    output_root = str(output_root)
    os.makedirs(output_root, exist_ok=True)

    receivers = _load_library(receiver_library_path)
    antennas = _load_library(antenna_library_path)

    rows: List[dict] = []
    for opath in rinex_files:
        opath = str(opath)
        nav_letter = "n" if opath[-1] == "o" else "N"
        npath = opath[:-1] + nav_letter
        record = {"origin path": opath}
        try:
            lines = read_lines(opath)
            header = read_obs_header(opath)
            record["version"] = str(header.version)

            first_time = header.TIME_OF_FIRST_OBS[0].replace("-", " ") if header.TIME_OF_FIRST_OBS else ""
            try:
                tokens = first_time.split()
                doy = date_to_doy(int(tokens[0]), int(tokens[1]), int(tokens[2]))
                doy_str = f"{doy:03d}"
            except (IndexError, ValueError):
                doy_str = "xxx"

            # Marker name
            marker_index = header.MARKER_NAME[0] if header.MARKER_NAME else None
            origin_marker = header.MARKER_NAME[1] if header.MARKER_NAME else "x"
            record["origin marker"] = origin_marker
            new_marker = (origin_marker[0:4] or "aaaa").replace(" ", "a")
            new_marker = (new_marker + "aaaa")[:4]
            record["marker"] = new_marker
            if marker_index is not None:
                lines[marker_index] = "{0:<60}".format(new_marker) + "MARKER NAME\n"

            # Coordinates → lon/lat
            x, y, z = header.APPROX_POSITION_XYZ
            lon, lat, _h = xyz_to_blh(x, y, z, radians=False)
            record["longitude"] = float(lon)
            record["latitude"] = float(lat)

            # Strip non-English characters in the first 20 header lines
            non_english = False
            for k in range(min(20, len(lines))):
                cleaned, hit = _strip_non_english(lines[k])
                if hit:
                    non_english = True
                    lines[k] = cleaned
            if non_english:
                record["non English"] = "Non-English characters present"

            # Receiver / antenna corrections
            rec_idx = header.receiver_type[0] if header.receiver_type else None
            ant_idx = header.antenna_type[0] if header.antenna_type else None
            if rec_idx is not None:
                record["origin rec"] = lines[rec_idx][20:40]
                rec_field = header.receiver_type[1]
                for wrong, right in receivers.items():
                    if "{0:<20}".format(wrong) == rec_field:
                        lines[rec_idx] = lines[rec_idx].replace(rec_field, "{0:<20}".format(right))
                        break
                record["rec type"] = lines[rec_idx][20:40]
            if ant_idx is not None:
                record["origin ant"] = lines[ant_idx][20:40]
                ant_field = header.antenna_type[1]
                for wrong, right in antennas.items():
                    if "{0:<20}".format(wrong) == ant_field:
                        lines[ant_idx] = lines[ant_idx].replace(ant_field, "{0:<20}".format(right))
                        break
                if radome:
                    ant_rad = header.antenna_type[1][16:20] if len(header.antenna_type[1]) >= 20 else ""
                    if not ant_rad.strip():
                        line = lines[ant_idx]
                        if len(line) >= 40:
                            lines[ant_idx] = line[:36] + "NONE" + line[40:]
                record["ant type"] = lines[ant_idx][20:40]

            # Output path
            new_folder = os.path.join(output_root, opath[-3:-1] + doy_str)
            os.makedirs(new_folder, exist_ok=True)
            out_o = os.path.join(
                new_folder, f"{new_marker}{doy_str}0.{opath[-3:]}"
            )
            out_o = unique_filename(out_o)
            record["new path"] = out_o
            with open(out_o, "w", encoding="utf-8") as fp:
                fp.writelines(lines)

            if os.path.exists(npath):
                nav_lines = read_lines(npath)
                nav_letter_out = "n" if out_o[-1] == "o" else "N"
                out_n = out_o[:-1] + nav_letter_out
                with open(out_n, "w", encoding="utf-8") as fp:
                    fp.writelines(nav_lines)

        except UnicodeDecodeError:
            record["new path"] = "UnicodeDecodeError"
            log.warning("UnicodeDecodeError on %s", opath)
        except (KeyError, AttributeError) as exc:
            record["new path"] = type(exc).__name__
            log.warning("%s on %s: %s", type(exc).__name__, opath, exc)
        rows.append(record)

    report_path = os.path.join(output_root, "report.csv")
    with open(report_path, "w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# v3.x backwards-compatible aliases
# ---------------------------------------------------------------------------

DataFinding = find_rinex
DataCleaning = clean_rinex
show_files = list_files
readtxt = read_lines
detectCode = _detect_encoding
date2doy = date_to_doy


def xyz2BLH(x, y, z, rad: bool = True):  # noqa: D401, N802 — legacy name
    """Deprecated alias of :func:`PyRINEX.coordinates.xyz_to_blh`. Returns
    ``[lon, lat, h]`` to match the v3.x ordering."""
    lon, lat, h = xyz_to_blh(x, y, z, radians=rad)
    return [lon, lat, h]


def generate_unique_filename(path):  # noqa: D401 — legacy name
    """Deprecated alias of :func:`unique_filename`."""
    return unique_filename(path)


def file_name_extesion_judgement(file_dir, extension_name):  # noqa: D401, N802
    """Deprecated alias of :func:`filter_by_extension` (typo preserved)."""
    return filter_by_extension(file_dir, extension_name)


def file_name_extesion_judgement_file(file_dir, extension_name):  # noqa: D401, N802
    """Like :func:`file_name_extesion_judgement` but returns bare filenames."""
    ext_l = extension_name.lower()
    out: List[str] = []
    for _root, _dirs, files in os.walk(file_dir):
        for f in files:
            if os.path.splitext(f)[1].lower() == ext_l:
                out.append(f)
    return out
