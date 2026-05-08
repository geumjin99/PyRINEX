"""Synthetic RINEX fixture generators for the test suite.

The published v3.x test data folder was removed from the repository, so
the test suite generates minimal-but-valid RINEX 2/3 obs and nav files
on the fly. The data is *just enough* to exercise the parsers and the
downstream math; numerical values are not physically meaningful.
"""
from __future__ import annotations

from pathlib import Path

# Tiny RINEX 2.11 GPS observation file. Two epochs, three satellites,
# observation types L1 C1 L2 P2 D1.
RINEX_V2_OBS = """\
     2.11           OBSERVATION DATA    G                   RINEX VERSION / TYPE
PYRINEX TEST                            2024-01-01 00:00:00 PGM / RUN BY / DATE
TESTSTATION                                                 MARKER NAME
TEST0001                                                    MARKER NUMBER
TEST_AGENCY         TEST_OBSERVER                           OBSERVER / AGENCY
0001                TRIMBLE 5700        1.0                 REC # / TYPE / VERS
0001                TRM39105.00         NONE                ANT # / TYPE
 -3173543.0196   4134173.2686   3666275.4904                APPROX POSITION XYZ
        0.0000        0.0000        0.0000                  ANTENNA: DELTA H/E/N
     1     1                                                WAVELENGTH FACT L1/2
     5    L1    C1    L2    P2    D1                        # / TYPES OF OBSERV
    30.0000                                                 INTERVAL
  2024     1     1     0     0    0.0000000     GPS         TIME OF FIRST OBS
                                                            END OF HEADER
 24  1  1  0  0  0.0000000  0  3G01G02G03
  -860455.086 6 25655170.640    -670354.860 7 25655178.040 6   1500.000
  -772302.852 6 25071945.022    -601876.150 7 25071950.000 6   1480.000
  -683993.539 6 25088749.890    -533090.450 7 25088755.000 6   1460.000
 24  1  1  0  0 30.0000000  0  3G01G02G03
  -595526.402 6 25105585.069    -464189.420 7 25105590.000 6   1500.000
  -506900.594 6 25122449.429    -395127.620 7 25122455.000 6   1480.000
  -418120.063 6 25139344.640    -325873.910 7 25139350.000 6   1460.000
"""

# Minimal RINEX 2 GPS navigation file with one ephemeris record per
# satellite. Values are physically reasonable broadcast ephemerides for
# GPS sat #1 mid-2024 (orbital period ≈ 12 h, eccentricity ≈ 0.012).
RINEX_V2_NAV = """\
     2.11           N: GPS NAV DATA                         RINEX VERSION / TYPE
PYRINEX TEST                            2024-01-01 00:00:00 PGM / RUN BY / DATE
   .1118D-07   .2235D-07  -.5960D-07  -.1192D-06          ION ALPHA
   .1106D+06   .9830D+05  -.1310D+06  -.1966D+06          ION BETA
                                                            END OF HEADER
 1 24  1  1  0  0  0.0 .619990751147D-04 .250111042648D-11 .000000000000D+00
    .900000000000D+02-.343125000000D+02 .457961000000D-08 .268174000000D+01
   -.182539225817D-05 .120000000000D-01 .102236867130D-04 .515365475006D+04
    .172800000000D+06 .819563865662D-07-.143119876958D+01-.819563865662D-07
    .967452236144D+00 .180812500000D+03 .714159265358D+00-.815000000000D-08
   -.342853000000D-09 .100000000000D+01 .229200000000D+04 .000000000000D+00
    .200000000000D+01 .000000000000D+00-.111758708954D-07 .900000000000D+02
    .172776000000D+06
 2 24  1  1  0  0  0.0 .244840000000D-03 .250111042648D-11 .000000000000D+00
    .900000000000D+02-.343125000000D+02 .457961000000D-08 .268174000000D+01
   -.182539225817D-05 .120000000000D-01 .102236867130D-04 .515365475006D+04
    .172800000000D+06 .819563865662D-07-.143119876958D+01-.819563865662D-07
    .967452236144D+00 .180812500000D+03 .714159265358D+00-.815000000000D-08
   -.342853000000D-09 .100000000000D+01 .229200000000D+04 .000000000000D+00
    .200000000000D+01 .000000000000D+00-.111758708954D-07 .900000000000D+02
    .172776000000D+06
 3 24  1  1  0  0  0.0 .977194570000D-04 .250111042648D-11 .000000000000D+00
    .900000000000D+02-.343125000000D+02 .457961000000D-08 .268174000000D+01
   -.182539225817D-05 .120000000000D-01 .102236867130D-04 .515365475006D+04
    .172800000000D+06 .819563865662D-07-.143119876958D+01-.819563865662D-07
    .967452236144D+00 .180812500000D+03 .714159265358D+00-.815000000000D-08
   -.342853000000D-09 .100000000000D+01 .229200000000D+04 .000000000000D+00
    .200000000000D+01 .000000000000D+00-.111758708954D-07 .900000000000D+02
    .172776000000D+06
"""

# Tiny RINEX 3.04 mixed-system observation file. Two epochs, two GPS
# satellites with three observation types (C1C, L1C, S1C).
RINEX_V3_OBS = """\
     3.04           OBSERVATION DATA    M                   RINEX VERSION / TYPE
PYRINEX TEST        TEST                2024-01-01 00:00:00 PGM / RUN BY / DATE
TESTSTATION                                                 MARKER NAME
GEODETIC                                                    MARKER TYPE
TEST_OBSERVER       TEST_AGENCY                             OBSERVER / AGENCY
0001                TRIMBLE NETR9       1.0                 REC # / TYPE / VERS
0001                TRM59800.00         NONE                ANT # / TYPE
 -3173543.0196   4134173.2686   3666275.4904                APPROX POSITION XYZ
        0.0000        0.0000        0.0000                  ANTENNA: DELTA H/E/N
G    3 C1C L1C S1C                                          SYS / # / OBS TYPES
    30.000                                                  INTERVAL
  2024     1     1     0     0    0.0000000     GPS         TIME OF FIRST OBS
                                                            END OF HEADER
> 2024  1  1  0  0  0.0000000  0  2
G01  25655170.640   -860455.086 6        45.000
G02  25071945.022   -772302.852 6        44.000
> 2024  1  1  0  0 30.0000000  0  2
G01  25105585.069   -595526.402 6        45.000
G02  25122449.429   -506900.594 6        44.000
"""


def write_v2_obs(folder: Path) -> Path:
    path = folder / "test0010.24o"
    path.write_text(RINEX_V2_OBS)
    return path


def write_v2_nav(folder: Path) -> Path:
    path = folder / "test0010.24n"
    path.write_text(RINEX_V2_NAV)
    return path


def write_v3_obs(folder: Path) -> Path:
    path = folder / "test_v3.rnx"
    path.write_text(RINEX_V3_OBS)
    return path
