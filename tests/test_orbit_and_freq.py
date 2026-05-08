"""Orbit-propagation and frequency-table sanity checks."""
from __future__ import annotations

import math

import pytest

from PyRINEX.frequencies import GLONASS_SLOT, frequency_pair
from PyRINEX.orbit import _to_float, gps_satellite_position


def test_frequency_pair_gps():
    pair = frequency_pair("G01")
    assert pair is not None
    assert math.isclose(pair.wavelength_1, 0.190293673, abs_tol=1e-6)
    assert pair.gamma > 1
    assert pair.alpha_iono > 0


def test_frequency_pair_glonass_uses_slot_table():
    pair = frequency_pair("R03")
    assert pair is not None
    # GLONASS slot 5 → f1 ≈ (1602 + 5 * 9/16) MHz
    f1 = (1602e6 + GLONASS_SLOT["R03"] * 9.0 / 16.0 * 1e6)
    assert math.isclose(0.299792458e9 / pair.wavelength_1, f1, rel_tol=1e-9)


def test_frequency_pair_unknown_returns_none():
    assert frequency_pair("X99") is None


def test_to_float_handles_fortran_d_notation():
    assert math.isclose(_to_float("0.515365475006D+04"), 5153.65475006, rel_tol=1e-9)
    assert math.isclose(_to_float("-.182539225817D-05"), -1.82539225817e-6, rel_tol=1e-9)


def test_gps_satellite_position_yields_orbital_radius():
    # Plausible GPS broadcast ephemeris (matches the synthetic nav fixture).
    eph = {
        "sqrt(A)": ".515365475006D+04",
        "e Eccentricity": ".120000000000D-01",
        "i0": ".967452236144D+00",
        "OMEGA0": "-.143119876958D+01",
        "omega": ".714159265358D+00",
        "M0": ".268174000000D+01",
        "Delta n": ".457961000000D-08",
        "IDOT": "-.815000000000D-08",
        "OMEGA DOT": "-.815000000000D-08",
        "Cuc": "-.182539225817D-05",
        "Cus": ".102236867130D-04",
        "Crc": ".180812500000D+03",
        "Crs": "-.343125000000D+02",
        "Cic": ".819563865662D-07",
        "CIS": "-.819563865662D-07",
        "Toe Time of Ephemeris": ".172800000000D+06",
    }
    x, y, z = gps_satellite_position(eph, t_k=0.0)
    radius = math.sqrt(x * x + y * y + z * z)
    # GPS satellites orbit at about 26 600 km from Earth's centre.
    assert 25e6 < radius < 28e6


def test_kepler_solver_converges_for_high_eccentricity():
    eph = {
        "sqrt(A)": ".515365475006D+04",
        "e Eccentricity": ".200000000000D+00",  # High ecc — non-physical but
        "i0": ".967452236144D+00",
        "OMEGA0": "-.143119876958D+01",
        "omega": ".714159265358D+00",
        "M0": ".268174000000D+01",
        "Delta n": ".457961000000D-08",
        "IDOT": "-.815000000000D-08",
        "OMEGA DOT": "-.815000000000D-08",
        "Cuc": "0.0D+00",
        "Cus": "0.0D+00",
        "Crc": "0.0D+00",
        "Crs": "0.0D+00",
        "Cic": "0.0D+00",
        "CIS": "0.0D+00",
        "Toe Time of Ephemeris": ".172800000000D+06",
    }
    # Should not loop forever or crash for any t_k in the validity window.
    for t_k in (-7200, -1000, 0, 1000, 7200):
        x, y, z = gps_satellite_position(eph, t_k)
        assert math.isfinite(x) and math.isfinite(y) and math.isfinite(z)
