"""Time conversions used by the GNSS quality-check routines.

The original code used a poorly-named ``utf2*`` family of helpers operating
on a list ``[yy, mm, dd, hh, mi, ss]`` parsed out of the RINEX epoch line.
We keep the same input shape for backwards compatibility but expose
clearer names.
"""
from __future__ import annotations

import math
from typing import Sequence, Tuple

GPS_EPOCH_MJD = 44244.0  # 1980-01-06 00:00:00 UTC
SECONDS_PER_DAY = 86400
SECONDS_PER_WEEK = 7 * SECONDS_PER_DAY


def _normalize_year(year_token: str) -> int:
    """Expand a 2-digit RINEX year to 4 digits, leaving 4-digit values intact."""
    if len(year_token) == 2:
        year_token = "20" + year_token
    return int(year_token)


def utc_to_mjd(utc: Sequence) -> float:
    """Convert ``[yy, mm, dd, hh, mi, ss]`` to Modified Julian Date.

    Implements the standard Julian-day formula with a Gregorian/Julian
    cutoff at 1582-10-15.
    """
    year = _normalize_year(str(utc[0]))
    month = int(utc[1])
    day = int(utc[2]) + (
        int(utc[3]) * 3600 + int(utc[4]) * 60 + float(utc[5])
    ) / SECONDS_PER_DAY

    y = year + 4800
    m = month
    if m <= 2:
        m += 12
        y -= 1
    e = math.floor(30.6 * (m + 1))
    a = math.floor(y / 100)
    if (year < 1582) or (year == 1582 and month < 10) or (
        year == 1582 and month == 10 and day < 15
    ):
        b = -38
    else:
        b = math.floor((a / 4) - a)
    c = math.floor(365.25 * y)
    jd = b + c + e + day - 32167.5
    return jd - 2400000.5


def utc_to_gps(utc: Sequence) -> Tuple[int, int]:
    """Convert ``[yy, mm, dd, hh, mi, ss]`` to ``(gps_week, seconds_of_week)``."""
    mjd = utc_to_mjd(utc)
    elapsed_days = mjd - GPS_EPOCH_MJD
    week = int(math.floor(elapsed_days / 7.0))
    seconds = round((elapsed_days - week * 7) * SECONDS_PER_DAY)
    return week, seconds


def time_from_ephemeris(seconds_of_week: float, toe: float) -> float | None:
    """Return ``t_k = t - toe`` rolled into ``±302400`` s, or ``None`` if outside
    the broadcast ephemeris validity window of ±2 hours.
    """
    t_k = seconds_of_week - toe
    if t_k > SECONDS_PER_WEEK / 2:
        t_k -= SECONDS_PER_WEEK
    elif t_k < -SECONDS_PER_WEEK / 2:
        t_k += SECONDS_PER_WEEK
    if -7200.0 <= t_k <= 7200.0:
        return t_k
    return None


def date_to_doy(year: int, month: int, day: int) -> int:
    """Return the day-of-year (1–366) for a Gregorian calendar date."""
    if not 1 <= month <= 12:
        raise ValueError(f"month out of range: {month}")
    leap = year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)
    days = [31, 29 if leap else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    return sum(days[: month - 1]) + day
