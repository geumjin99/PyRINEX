"""Round-trip and consistency checks for the WGS-84 coordinate helpers."""
from __future__ import annotations

import math

import numpy as np
import pytest

from PyRINEX.coordinates import blh_to_xyz, enu_basis, xyz_to_blh


@pytest.mark.parametrize(
    "x,y,z",
    [
        (-3173543.0196, 4134173.2686, 3666275.4904),  # Suwon, Korea
        (4642417.31, 121519.42, 4364805.20),  # Paris
        (-2710830.0, -4307229.0, 3851666.0),  # Boulder, US
    ],
)
def test_xyz_blh_round_trip(x, y, z):
    lon, lat, h = xyz_to_blh(x, y, z, radians=True)
    x_back, y_back, z_back = blh_to_xyz(lon, lat, h)
    assert math.isclose(float(x_back), x, abs_tol=1e-3)
    assert math.isclose(float(y_back), y, abs_tol=1e-3)
    assert math.isclose(float(z_back), z, abs_tol=1e-3)


def test_xyz_to_blh_returns_degrees_when_requested():
    lon, lat, _h = xyz_to_blh(-3173543.0196, 4134173.2686, 3666275.4904, radians=False)
    assert -180 <= float(lon) <= 180
    assert -90 <= float(lat) <= 90
    # The fixture XYZ resolves to roughly (127E, 35N), eastern China.
    assert math.isclose(float(lon), 127.5, abs_tol=2.0)
    assert math.isclose(float(lat), 35.5, abs_tol=2.0)


def test_enu_basis_is_orthonormal():
    basis = enu_basis(-3173543.0196, 4134173.2686, 3666275.4904)
    assert basis.shape == (3, 3)
    product = basis @ basis.T
    np.testing.assert_allclose(product, np.eye(3), atol=1e-9)
