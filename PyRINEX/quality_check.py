"""GNSS RINEX quality-check routines.

This is a rewrite of the v3.x ``QulityCheck.py`` (note the typo, which is
preserved in :mod:`PyRINEX.QulityCheck` as a deprecation shim). The
public surface is the same — :func:`quality_check`, :func:`batch_quality_check`,
:func:`multipath_iod`, :func:`cycle_slip`, :func:`azimuth_elevation`,
:func:`satellite_signal_plot` — but:

* the four hand-unrolled per-constellation extraction blocks (G/R/S/E,
  ~50 lines each) are replaced by a config-driven loop;
* the ``abs(x > threshold)`` paren bug in the old MP1/MP2 outlier guard is
  fixed;
* the ``fieldnames = prns; fieldnames.insert(0, "epoch")`` aliasing bug
  that silently shifted PRN indices is fixed;
* the ``json.dumps``/``json.loads`` round-trip between every call is gone;
* dead commented-out code has been removed.

Implementation details follow the equations laid out in the LaTeX manual
under ``docs/manual.tex``.
"""
from __future__ import annotations

import csv
import math
import os
from typing import Dict, List, Mapping, Sequence, Tuple

import numpy as np

from . import plotting
from .coordinates import enu_basis
from .frequencies import frequency_pair
from .orbit import _to_float, gps_satellite_position
from .reader import (
    EMPTY_OBS,
    PathLike,
    read_nav,
    read_obs,
    read_obs_header,
)
from .time_utils import time_from_ephemeris, utc_to_gps

# Per-constellation observation-band layout used by :func:`extract_l_c`.
# The mapping says: for satellites whose PRN starts with ``code``, slot
# ``0`` is L-band 1, slot ``1`` is L-band 2, slot ``2`` is C/P range 1,
# slot ``3`` is C/P range 2. Each entry of the inner tuple is a list of
# RINEX observation-type prefixes that we accept for that slot.
_BAND_LAYOUT: Dict[str, Tuple[Tuple[str, ...], ...]] = {
    "G": (("L1",), ("L2",), ("C1", "P1"), ("C2", "P2")),
    "R": (("L1",), ("L2",), ("C1", "P1"), ("C2", "P2")),
    "S": (("L1",), ("L5",), ("C1",), ("C5",)),
    "E": (("L1",), ("L5",), ("C1",), ("C5",)),
}

_MULTIPATH_LIMIT = 10  # threshold ``maxoff`` from the original code


# ---------------------------------------------------------------------------
# Time / geometry
# ---------------------------------------------------------------------------

def epoch_labels(observations: Mapping[str, dict]) -> List[str]:
    """Return the epoch keys of an observation dict in insertion order."""
    return list(observations.keys())


def gps_positions(opath: PathLike) -> List[dict]:
    """Compute the ECEF position of every GPS satellite at every epoch.

    Pairs each observation epoch with the broadcast ephemeris closest in
    time-of-ephemeris.
    """
    obs = read_obs(opath)
    nav_letter = "n" if str(opath)[-1] == "o" else "N"
    nav_path = str(opath)[:-1] + nav_letter
    nav = read_nav(nav_path)

    out: List[dict] = []
    for label, epoch in obs.items():
        utc = label.split()
        seconds_of_week = utc_to_gps(utc)[1]
        record = {"epoch": label}
        for prn in epoch:
            if not isinstance(epoch[prn], dict) or not prn.startswith("G"):
                continue
            ephemerides = nav.get(prn)
            if not ephemerides:
                continue
            best: Tuple[float, dict] | None = None
            for eph in ephemerides:
                try:
                    toe = _to_float(eph["Toe Time of Ephemeris"])
                except (KeyError, ValueError):
                    continue
                t_k = time_from_ephemeris(seconds_of_week, toe)
                if t_k is None:
                    continue
                if best is None or abs(t_k) < abs(best[0]):
                    best = (t_k, eph)
            if best is None:
                continue
            t_k, eph = best
            record[prn] = list(gps_satellite_position(eph, t_k))
        out.append(record)
    return out


def azimuth_elevation(opath: PathLike) -> List[dict]:
    """Compute (azimuth, elevation) of every GPS satellite per epoch.

    Writes a sky-plot PNG and an ``aziele.csv`` next to ``opath``.
    Returns the per-epoch list of dicts.
    """
    header = read_obs_header(opath)
    xyz = header.APPROX_POSITION_XYZ
    basis = enu_basis(*xyz)
    positions = gps_positions(opath)

    records: List[dict] = []
    for entry in positions:
        record = {"epoch": entry["epoch"]}
        for prn, pos in entry.items():
            if prn == "epoch":
                continue
            d = basis @ (np.array(pos) - np.array(xyz))
            d_norm = math.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2)
            sin_el = d[2] / d_norm if d_norm else 0.0
            sin_el = max(-1.0, min(1.0, sin_el))
            elevation = math.asin(sin_el)
            if d[1] == 0:
                azimuth = math.pi / 2 if d[0] > 0 else 1.5 * math.pi
            else:
                azimuth = math.atan(d[0] / d[1])
                if d[1] < 0:
                    azimuth += math.pi
                elif d[0] < 0:
                    azimuth += 2 * math.pi
            record[prn] = (azimuth, elevation)
        records.append(record)

    plotting.sky_plot(str(opath), records, header.PRNS)
    _write_csv(
        str(opath)[:-4] + "aziele.csv",
        ["epoch"] + [p for p in header.PRNS if p.startswith("G")],
        [
            {"epoch": r["epoch"], **{k: v for k, v in r.items() if k != "epoch"}}
            for r in records
        ],
    )
    return records


# ---------------------------------------------------------------------------
# Multipath / ionosphere / cycle-slip
# ---------------------------------------------------------------------------

def _first_value(obs_dict: Mapping[str, str], prefixes: Sequence[str]) -> float:
    """Return the first non-empty observation among ``prefixes``, else 0."""
    for prefix in prefixes:
        for obs_type, raw in obs_dict.items():
            if prefix in obs_type and raw != EMPTY_OBS:
                try:
                    return float(raw[0:14])
                except ValueError:
                    return 0.0
    return 0.0


def extract_l_c(opath: PathLike) -> np.ndarray:
    """Extract the L1, L2, P1/C1, P2/C2 cycle/range pairs into an
    ``(n_epochs, n_prns, 4)`` numpy array.

    Replaces the four 50-line blocks in the original ``LandC`` function.
    """
    obs = read_obs(opath)
    header = read_obs_header(opath)
    prns = header.PRNS
    epochs = list(obs.keys())
    out = np.zeros((len(epochs), len(prns), 4), dtype=float)

    for epoch_idx, epoch_label in enumerate(epochs):
        epoch = obs[epoch_label]
        for prn, sat in epoch.items():
            if not isinstance(sat, dict):
                continue
            if prn not in prns:
                continue
            layout = _BAND_LAYOUT.get(prn[0])
            if layout is None:
                continue
            sat_idx = prns.index(prn)
            for band_idx, prefixes in enumerate(layout):
                out[epoch_idx, sat_idx, band_idx] = _first_value(sat, prefixes)
    return out


def _frequency_table(prns: Sequence[str]) -> np.ndarray:
    """Build the ``(n_prns, 4)`` numpy table of (λ₁, λ₂, γ, α_iono).

    Replaces the ``frelist`` function from the original module.
    """
    fres = np.zeros((len(prns), 4), dtype=float)
    for i, prn in enumerate(prns):
        pair = frequency_pair(prn)
        if pair is None:
            continue
        fres[i, 0] = pair.wavelength_1
        fres[i, 1] = pair.wavelength_2
        fres[i, 2] = pair.gamma
        fres[i, 3] = pair.alpha_iono
    return fres


def multipath_iod(opath: PathLike) -> np.ndarray:
    """Compute multipath (MP1, MP2), ionospheric delay (ION) and ionospheric
    drift (IOD) per (epoch, PRN), and emit the corresponding plots and CSV
    files alongside ``opath``.

    Returns a ``(n_epochs, n_prns, 4)`` array of ``[MP1, MP2, ION, IOD]``.
    """
    obs = read_obs(opath)
    header = read_obs_header(opath)
    prns = list(header.PRNS)
    fres = _frequency_table(prns)
    lc = extract_l_c(opath)

    n_epochs = lc.shape[0]
    n_prns = lc.shape[1]
    epochs = list(obs.keys())
    result = np.zeros((n_epochs, n_prns, 4), dtype=float)
    initial_iono = np.zeros(n_prns)

    for n in range(n_epochs):
        for m in range(n_prns):
            if 0.0 in lc[n, m]:
                continue
            l1, l2, gamma, _alpha = fres[m]
            l1_m = lc[n, m, 0] * l1
            l2_m = lc[n, m, 1] * l2
            denom = gamma - 1
            if denom == 0:
                continue
            result[n, m, 0] = lc[n, m, 2] - l1_m * (2 / denom + 1) + l2_m * (2 / denom)
            result[n, m, 1] = (
                lc[n, m, 3] - l1_m * (2 * gamma / denom) + l2_m * (2 * gamma / denom - 1)
            )
            current_iono = (gamma / denom) * (l1_m - l2_m)
            if initial_iono[m] == 0:
                initial_iono[m] = current_iono
                result[n, m, 2] = 0
            else:
                result[n, m, 2] = initial_iono[m] - current_iono
                if abs(result[n, m, 2]) > 10 * _MULTIPATH_LIMIT:
                    initial_iono[m] = current_iono
                if n > 0:
                    dt = (
                        utc_to_gps(epochs[n].split())[1]
                        - utc_to_gps(epochs[n - 1].split())[1]
                    )
                    if dt:
                        result[n, m, 3] = (
                            (result[n, m, 2] - result[n - 1, m, 2]) / dt * 60
                        )

    _detrend_multipath(result, n_epochs, n_prns)
    _write_multipath_outputs(opath, prns, epochs, result)
    return result


def _detrend_multipath(result: np.ndarray, n_epochs: int, n_prns: int) -> None:
    """Subtract the mean of each contiguous segment of MP1/MP2 from itself.

    Direct port of the segment loop in the original ``ION_MP``, kept here
    because the algorithm is intentional (it isolates antenna multipath
    bias by removing the per-segment ambiguity bias).
    """
    mp1_mean = np.zeros((n_epochs, n_prns))
    mp2_mean = np.zeros((n_epochs, n_prns))
    for j in range(n_prns):
        i = 0
        while i < n_epochs - 1:
            if result[i, j, 0] == 0:
                i += 1
                continue
            if abs(result[i, j, 0] - result[i + 1, j, 0]) < 2 * _MULTIPATH_LIMIT:
                start = i
                ave1 = ave2 = 0.0
                while i < n_epochs - 1 and abs(
                    result[i, j, 0] - result[i + 1, j, 0]
                ) < 2 * _MULTIPATH_LIMIT:
                    ave1 += result[i, j, 0]
                    ave2 += result[i, j, 1]
                    i += 1
                length = i - start
                if length > 0:
                    ave1 /= length
                    ave2 /= length
                    mp1_mean[start : i + 1, j] = ave1
                    mp2_mean[start : i + 1, j] = ave2
            else:
                i += 1

    nonzero = result[:, :, 0] != 0
    result[:, :, 0] = np.where(nonzero, result[:, :, 0] - mp1_mean, 0)
    result[:, :, 1] = np.where(nonzero, result[:, :, 1] - mp2_mean, 0)


def _write_multipath_outputs(
    opath: PathLike,
    prns: List[str],
    epochs: Sequence[str],
    result: np.ndarray,
) -> None:
    base = str(opath)[:-4]
    fieldnames = ["epoch", *prns]

    rows_mp: List[dict] = []
    rows_ion: List[dict] = []
    for n, epoch in enumerate(epochs):
        row_mp = {"epoch": epoch}
        row_ion = {"epoch": epoch}
        for m, prn in enumerate(prns):
            if result[n, m, 0] != 0 and result[n, m, 1] != 0:
                row_mp[prn] = [result[n, m, 0], result[n, m, 1]]
            if result[n, m, 2] != 0 and result[n, m, 3] != 0:
                row_ion[prn] = [result[n, m, 2], result[n, m, 3]]
        rows_mp.append(row_mp)
        rows_ion.append(row_ion)

    _write_csv(base + "mp1mp2.csv", fieldnames, rows_mp)
    _write_csv(base + "IonIod.csv", fieldnames, rows_ion)

    name = os.path.basename(str(opath))
    columns = math.ceil(len(prns) / 18)
    plotting.scatter_per_satellite(
        base + "MP1_plot.png",
        name[:-4] + " MP1 plot",
        prns,
        epochs,
        result,
        0,
        columns,
        y_label="MP1 (m)",
    )
    plotting.scatter_per_satellite(
        base + "MP2_plot.png",
        name[:-4] + " MP2 plot",
        prns,
        epochs,
        result,
        1,
        columns,
        y_label="MP2 (m)",
    )
    plotting.scatter_per_satellite(
        base + "ION_plot.png",
        name[:-4] + " ION plot",
        prns,
        epochs,
        result,
        2,
        columns,
        y_label="ION (m)",
    )
    plotting.scatter_per_satellite(
        base + "IOD_plot.png",
        name[:-4] + " IOD plot",
        prns,
        epochs,
        result,
        3,
        columns,
        y_label="IOD (m)",
    )


def cycle_slip(opath: PathLike) -> np.ndarray:
    """Detect carrier-phase cycle slips via the second-difference of the
    geometry-free combination, then write the cycle-slip plot and CSV.

    Returns a ``(n_epochs - 2, n_prns, 4)`` array; the last band index is
    the second-difference indicator used for plotting.
    """
    obs = read_obs(opath)
    header = read_obs_header(opath)
    prns = list(header.PRNS)
    fres = _frequency_table(prns)
    lc = extract_l_c(opath)

    n_epochs = lc.shape[0]
    n_prns = lc.shape[1]

    last_p = np.zeros(n_prns)
    last_phi = np.zeros(n_prns)
    p = np.zeros((n_epochs, n_prns))
    phi = np.zeros((n_epochs, n_prns))
    delta_p = np.zeros((n_epochs, n_prns))
    delta_phi = np.zeros((n_epochs, n_prns))

    for n in range(n_epochs):
        for m in range(n_prns):
            if lc[n, m, 0] == 0:
                continue
            alpha1 = fres[m, 3]
            alpha2 = 1 - alpha1
            lambda1 = fres[m, 0]
            lambda2 = fres[m, 1]
            p[n, m] = alpha1 * lc[n, m, 2] + alpha2 * lc[n, m, 3]
            phi[n, m] = (
                alpha1 * lambda1 * lc[n, m, 0] + alpha2 * lambda2 * lc[n, m, 1]
            )
            if last_p[m] != 0:
                delta_p[n, m] = p[n, m] - last_p[m]
                delta_phi[n, m] = phi[n, m] - last_phi[m]
            last_p[m] = p[n, m]
            last_phi[m] = phi[n, m]

    if n_epochs < 4:
        return np.zeros((max(n_epochs - 2, 0), n_prns, 4))

    cycleslips = np.zeros((n_epochs - 2, n_prns, 4))

    for n in range(1, n_epochs - 1):
        for m in range(n_prns):
            second_diff_p = (delta_p[n - 1, m] + delta_p[n + 1, m]) - 2 * delta_p[n, m]
            if abs(second_diff_p) > 280000:
                second_diff_p -= math.copysign(299792.458, second_diff_p)
            cycleslips[n - 1, m, 0] = second_diff_p
            cycleslips[n - 1, m, 1] = (
                delta_phi[n - 1, m] + delta_phi[n + 1, m]
            ) - 2 * delta_phi[n, m]

    for n in range(1, n_epochs - 3):
        for m in range(n_prns):
            cycleslips[n - 1, m, 2] = (
                cycleslips[n - 1, m, 0] + cycleslips[n + 1, m, 0]
            ) - 2 * cycleslips[n, m, 0]
            cycleslips[n - 1, m, 3] = (
                cycleslips[n - 1, m, 1] + cycleslips[n + 1, m, 1]
            ) - 2 * cycleslips[n, m, 1]

    epochs = list(obs.keys())
    base = str(opath)[:-4]
    rows: List[dict] = []
    for n in range(len(cycleslips)):
        row = {"epoch": epochs[n]}
        for m, prn in enumerate(prns):
            if cycleslips[n, m, 2] != 0 and cycleslips[n, m, 3] != 0:
                row[prn] = cycleslips[n, m, 3]
        rows.append(row)
    _write_csv(base + "CScarrier.csv", ["epoch", *prns], rows)

    columns = math.ceil(len(prns) / 18)
    plotting.cycle_slip_plot(str(opath), epochs, cycleslips, prns, columns)
    return cycleslips


# ---------------------------------------------------------------------------
# Top-level orchestration
# ---------------------------------------------------------------------------

def satellite_signal_plot(opath: PathLike) -> None:
    """Render the per-PRN signal-presence raster (Figure 12 in the manual)."""
    obs = read_obs(opath)
    epochs: List[str] = []
    sats: List[str] = []
    for label, epoch in obs.items():
        for prn, sat in epoch.items():
            if isinstance(sat, dict):
                epochs.append(label)
                sats.append(prn)
    plotting.signal_plot(str(opath), epochs, sats)


def quality_check(opath: PathLike) -> None:
    """Run every QC routine that has the data it needs.

    If a navigation file matching ``opath`` exists, the satellite-position
    dependent step (:func:`azimuth_elevation`) runs as well; otherwise it
    is skipped silently.
    """
    nav_letter = "n" if str(opath)[-1] == "o" else "N"
    nav_path = str(opath)[:-1] + nav_letter
    cycle_slip(opath)
    multipath_iod(opath)
    if os.path.exists(nav_path):
        azimuth_elevation(opath)
    satellite_signal_plot(opath)


def batch_quality_check(
    root_path: PathLike,
    keywords: Sequence[str],
    extension: str,
) -> None:
    """Run :func:`quality_check` on every RINEX file under ``root_path``
    matching ``extension`` and the ``keywords`` filter."""
    from .data_management import find_rinex

    for path in find_rinex(root_path, keywords, extension):
        quality_check(path)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_csv(path: str, fieldnames: Sequence[str], rows: Sequence[dict]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(fieldnames))
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# v3.x backwards-compatible aliases
# ---------------------------------------------------------------------------

ION_MP = multipath_iod
cycleslip = cycle_slip
azi_ele = azimuth_elevation
SatelliteSignalPlot = satellite_signal_plot
QualityCheck = quality_check
batchQC = batch_quality_check
LandC = extract_l_c
frelist = _frequency_table
gps_pos = gps_positions
plot = plotting.scatter_per_satellite
time_serises = epoch_labels  # typo preserved for compatibility
