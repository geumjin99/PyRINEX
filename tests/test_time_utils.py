"""Tests for time conversions."""
from __future__ import annotations

import math

import pytest

from PyRINEX.time_utils import (
    date_to_doy,
    time_from_ephemeris,
    utc_to_gps,
    utc_to_mjd,
)


def test_utc_to_mjd_known_value():
    # 2024-01-01 00:00:00 UTC = MJD 60310
    mjd = utc_to_mjd(["24", "1", "1", "0", "0", "0"])
    assert math.isclose(mjd, 60310.0, abs_tol=1e-6)


def test_utc_to_gps_known_value():
    # 2024-01-01 00:00:00 UTC ≈ GPS week 2295, day 1 → 86400 s offset.
    week, seconds = utc_to_gps(["24", "1", "1", "0", "0", "0"])
    assert week == 2295
    assert seconds == 86400


def test_utc_to_gps_two_digit_and_four_digit_year_match():
    a = utc_to_gps(["24", "6", "15", "12", "0", "0"])
    b = utc_to_gps(["2024", "6", "15", "12", "0", "0"])
    assert a == b


@pytest.mark.parametrize(
    "year,month,day,expected",
    [
        (2024, 1, 1, 1),
        (2024, 12, 31, 366),  # leap
        (2023, 12, 31, 365),
        (2024, 3, 1, 61),
    ],
)
def test_date_to_doy(year, month, day, expected):
    assert date_to_doy(year, month, day) == expected


def test_date_to_doy_rejects_bad_month():
    with pytest.raises(ValueError):
        date_to_doy(2024, 13, 1)


def test_time_from_ephemeris_within_window_returns_value():
    assert time_from_ephemeris(100.0, 0.0) == 100.0


def test_time_from_ephemeris_outside_window_returns_none():
    assert time_from_ephemeris(20000.0, 0.0) is None


def test_time_from_ephemeris_handles_week_wrap():
    week_seconds = 7 * 86400
    # ``t`` 5 s ago in week-N expressed as "near the end of week N-1":
    near_end = week_seconds - 5
    assert math.isclose(time_from_ephemeris(0.0, near_end), 5.0, abs_tol=1e-6)
