"""GPS broadcast-ephemeris orbit propagation (Keplerian → ECEF)."""
from __future__ import annotations

import math
from typing import Mapping, Tuple

EARTH_GM = 3.986005e14  # m^3/s^2 — WGS-84
EARTH_OMEGA_DOT = 7.2921151467e-5  # rad/s — WGS-84

_KEPLER_MAX_ITER = 20
_KEPLER_TOL = 1e-12


def _to_float(value: str) -> float:
    """Parse a Fortran-style ``D``-exponent number from RINEX nav records."""
    return float(value.replace("D", "E").replace("d", "e"))


def gps_satellite_position(
    eph: Mapping[str, str],
    t_k: float,
) -> Tuple[float, float, float]:
    """Return GPS satellite ECEF position at time-from-ephemeris ``t_k`` (s).

    ``eph`` is one record as returned by :func:`PyRINEX.reader.read_nav` —
    a mapping from human-readable keys ("sqrt(A)", "Crs", "Cuc", …) to the
    raw 19-character Fortran-style string fields parsed from the navigation
    file.

    Implements IS-GPS-200, Table 20-IV (User Algorithm for Ephemeris
    Determination) without relativistic clock correction (clients should
    apply that separately if needed).
    """
    sqrt_a = _to_float(eph["sqrt(A)"])
    e = _to_float(eph["e Eccentricity"])
    i_0 = _to_float(eph["i0"])
    omega_0 = _to_float(eph["OMEGA0"])
    omega = _to_float(eph["omega"])
    m_0 = _to_float(eph["M0"])
    delta_n = _to_float(eph["Delta n"])
    i_dot = _to_float(eph["IDOT"])
    omega_dot = _to_float(eph["OMEGA DOT"])
    cuc = _to_float(eph["Cuc"])
    cus = _to_float(eph["Cus"])
    crc = _to_float(eph["Crc"])
    crs = _to_float(eph["Crs"])
    cic = _to_float(eph["Cic"])
    cis = _to_float(eph["CIS"])
    toe = _to_float(eph["Toe Time of Ephemeris"])

    a = sqrt_a ** 2
    n_0 = math.sqrt(EARTH_GM / a ** 3)
    n = n_0 + delta_n
    m_k = (m_0 + n * t_k) % (2 * math.pi)

    # Solve Kepler's equation iteratively (Newton-Raphson on E - e*sin E = M).
    e_k = m_k
    for _ in range(_KEPLER_MAX_ITER):
        f = e_k - e * math.sin(e_k) - m_k
        f_prime = 1.0 - e * math.cos(e_k)
        delta = f / f_prime
        e_k -= delta
        if abs(delta) < _KEPLER_TOL:
            break

    sin_e = math.sin(e_k)
    cos_e = math.cos(e_k)
    nu_k = math.atan2(math.sqrt(1.0 - e * e) * sin_e, cos_e - e)
    phi_k = nu_k + omega

    sin_2phi = math.sin(2.0 * phi_k)
    cos_2phi = math.cos(2.0 * phi_k)
    delta_u = cus * sin_2phi + cuc * cos_2phi
    delta_r = crs * sin_2phi + crc * cos_2phi
    delta_i = cis * sin_2phi + cic * cos_2phi

    u_k = phi_k + delta_u
    r_k = a * (1.0 - e * cos_e) + delta_r
    i_k = i_0 + i_dot * t_k + delta_i

    x_p = r_k * math.cos(u_k)
    y_p = r_k * math.sin(u_k)
    omega_k = omega_0 + (omega_dot - EARTH_OMEGA_DOT) * t_k - EARTH_OMEGA_DOT * toe

    sin_omega = math.sin(omega_k)
    cos_omega = math.cos(omega_k)
    cos_i = math.cos(i_k)
    sin_i = math.sin(i_k)

    x = x_p * cos_omega - y_p * cos_i * sin_omega
    y = x_p * sin_omega + y_p * cos_i * cos_omega
    z = y_p * sin_i
    return x, y, z
