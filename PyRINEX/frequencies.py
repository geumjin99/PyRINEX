"""Per-constellation carrier-frequency definitions and GLONASS slot table."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

SPEED_OF_LIGHT = 0.299792458e9  # m/s

# Frequencies (Hz)
GPS_L1 = 1.57542e9
GPS_L2 = 1.22760e9
GAL_E1 = 1.57542e9
GAL_E5a = 1.17645e9
SBAS_L1 = 1.57542e9
SBAS_L5 = 1.17645e9

# GLONASS L1/L2 use FDMA: f_k = f0 + k * delta_f.
# This table maps PRN to frequency channel number ``k``. Keep in sync with
# the ICD-defined slot assignments at the time of the original release;
# editable by users for newer almanacs.
GLONASS_SLOT: Dict[str, int] = {
    "R01": 1, "R02": -4, "R03": 5, "R04": 6,
    "R05": 1, "R06": -4, "R07": 5, "R08": 6,
    "R09": -2, "R10": -7, "R11": 0, "R12": -1,
    "R13": -2, "R14": -7, "R15": 0, "R16": -1,
    "R17": 4, "R18": -3, "R19": 3, "R20": 2,
    "R21": 4, "R22": -3, "R23": 3, "R24": 2,
}


@dataclass(frozen=True)
class FrequencyPair:
    """Two-frequency parameters used by ionosphere / multipath routines.

    Attributes:
        wavelength_1: λ₁, metres.
        wavelength_2: λ₂, metres.
        gamma:        (f₁/f₂)².
        alpha_iono:   f₁² / (f₁² − f₂²); the geometry-free combination
                      coefficient used to extract ionospheric delay.
    """

    wavelength_1: float
    wavelength_2: float
    gamma: float
    alpha_iono: float

    @classmethod
    def from_frequencies(cls, f1: float, f2: float) -> "FrequencyPair":
        return cls(
            wavelength_1=SPEED_OF_LIGHT / f1,
            wavelength_2=SPEED_OF_LIGHT / f2,
            gamma=(f1 / f2) ** 2,
            alpha_iono=f1 ** 2 / (f1 ** 2 - f2 ** 2),
        )


def frequency_pair(prn: str) -> FrequencyPair | None:
    """Return the two-frequency parameters for ``prn``, or ``None`` if the
    constellation isn't covered by this table.
    """
    system = prn[0]
    if system == "G":
        return FrequencyPair.from_frequencies(GPS_L1, GPS_L2)
    if system == "E":
        return FrequencyPair.from_frequencies(GAL_E1, GAL_E5a)
    if system == "S":
        return FrequencyPair.from_frequencies(SBAS_L1, SBAS_L5)
    if system == "R":
        k = GLONASS_SLOT.get(prn)
        if k is None:
            return None
        f1 = (1602e6 + k * 9.0 / 16.0 * 1e6)
        f2 = (1246e6 + k * 7.0 / 16.0 * 1e6)
        return FrequencyPair.from_frequencies(f1, f2)
    return None
