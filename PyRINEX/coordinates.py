"""Geodetic coordinate conversions on the WGS-84 ellipsoid."""
from __future__ import annotations

import math
from typing import Tuple

import numpy as np

WGS84_A = 6378137.0
WGS84_B = 6356752.3142
WGS84_E2 = 1.0 - (WGS84_B / WGS84_A) ** 2


def xyz_to_blh(
    x: float | np.ndarray,
    y: float | np.ndarray,
    z: float | np.ndarray,
    radians: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert ECEF (X, Y, Z) to geodetic (longitude, latitude, height).

    Returns ``(lon, lat, h)``. Angles in radians unless ``radians=False``.
    Height is in metres above the WGS-84 ellipsoid.

    Iterates the standard closed-form Bowring formulation; convergence
    threshold is 0.1 m on height.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)

    p = np.sqrt(x ** 2 + y ** 2)
    lon = np.arctan2(y, x)

    def _solve(zi: float, pi: float) -> Tuple[float, float]:
        b = math.atan2(zi, pi * (1.0 - WGS84_E2))
        h_prev = 1e9
        for _ in range(20):
            sin_b = math.sin(b)
            n = WGS84_A / math.sqrt(1.0 - WGS84_E2 * sin_b * sin_b)
            h = pi / math.cos(b) - n
            b = math.atan2(zi, pi * (1.0 - WGS84_E2 * n / (n + h)))
            if abs(h - h_prev) < 0.1:
                break
            h_prev = h
        return h, b

    h, lat = np.vectorize(_solve)(z, p)
    if radians:
        return lon, lat, h
    return np.degrees(lon), np.degrees(lat), h


def blh_to_xyz(
    lon: float | np.ndarray,
    lat: float | np.ndarray,
    h: float | np.ndarray,
    radians: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Inverse of :func:`xyz_to_blh`."""
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    h = np.asarray(h, dtype=float)
    if not radians:
        lon = np.radians(lon)
        lat = np.radians(lat)
    n = WGS84_A / np.sqrt(1.0 - WGS84_E2 * np.sin(lat) ** 2)
    x = (n + h) * np.cos(lat) * np.cos(lon)
    y = (n + h) * np.cos(lat) * np.sin(lon)
    z = (n * (1.0 - WGS84_E2) + h) * np.sin(lat)
    return x, y, z


def enu_basis(x: float, y: float, z: float) -> np.ndarray:
    """Return the local ENU rotation matrix at receiver position ``(x, y, z)``.

    Rows are East, North, Up unit vectors expressed in ECEF.
    """
    p = math.sqrt(x * x + y * y)
    r = math.sqrt(x * x + y * y + z * z)
    return np.array(
        [
            [-y / p, x / p, 0.0],
            [-x * z / (p * r), -y * z / (p * r), p / r],
            [x / r, y / r, z / r],
        ]
    )
