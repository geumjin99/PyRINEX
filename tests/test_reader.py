"""Reader tests against the synthetic RINEX fixtures."""
from __future__ import annotations

import json
import warnings

import pytest

from PyRINEX import reader


def test_read_v2_obs_header(v2_obs_path):
    header = reader.read_obs_header(v2_obs_path)
    assert header.version == 2
    assert header.type == "G"
    assert header.MARKER_NAME[1] == "TESTSTATION"
    assert header.MARKER_NUMBER[1] == "TEST0001"
    assert header.APPROX_POSITION_XYZ == [-3173543.0196, 4134173.2686, 3666275.4904]
    assert header.TIME_OF_FIRST_OBS == ["2024-1-1-0-0-0.0000000"]
    assert isinstance(header.ObsTypes, list)
    assert header.ObsTypes == ["L1", "C1", "L2", "P2", "D1"]
    assert "G01" in header.PRNS
    assert "G02" in header.PRNS
    assert "G03" in header.PRNS
    assert header.END_OF_HEADER >= 0


def test_read_v2_obs_body(v2_obs_path):
    epochs = reader.read_obs(v2_obs_path)
    assert len(epochs) == 2
    first_epoch = list(epochs.values())[0]
    assert first_epoch["sat_num"] == 3
    assert "G01" in first_epoch
    g01 = first_epoch["G01"]
    assert isinstance(g01, dict)
    assert set(g01) == {"L1", "C1", "L2", "P2", "D1"}
    # Field strings remain 16 chars wide
    assert len(g01["L1"]) == 16
    assert "-860455.086" in g01["L1"]


def test_read_v2_nav(v2_nav_path):
    nav = reader.read_nav(v2_nav_path)
    assert "G01" in nav
    assert "G02" in nav
    assert "G03" in nav
    record = nav["G01"][0]
    # 30 keys (NAV_KEYS_V2)
    assert "sqrt(A)" in record
    assert "Toe Time of Ephemeris" in record


def test_read_nav_iono(v2_nav_path):
    iono = reader.read_nav_iono(v2_nav_path)
    assert iono != ["NO ION data."]
    assert len(iono) == 8


def test_read_v3_obs_header(v3_obs_path):
    header = reader.read_obs_header(v3_obs_path)
    assert header.version == 3
    assert header.MARKER_NAME[1] == "TESTSTATION"
    assert isinstance(header.ObsTypes, dict)
    assert header.ObsTypes["G"] == ["C1C", "L1C", "S1C"]


def test_read_v3_obs_body(v3_obs_path):
    epochs = reader.read_obs(v3_obs_path)
    assert len(epochs) == 2
    first = list(epochs.values())[0]
    assert first["sat_num"] == 2
    assert "G01" in first
    g01 = first["G01"]
    assert set(g01) == {"C1C", "L1C", "S1C"}


def test_legacy_oheader_still_returns_json(v2_obs_path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        payload = reader.oheader(v2_obs_path)
    parsed = json.loads(payload)
    assert parsed["version"] == "2"
    assert parsed["MARKER_NAME"][1] == "TESTSTATION"


def test_legacy_oheader_emits_deprecation_warning(v2_obs_path):
    with pytest.warns(DeprecationWarning):
        reader.oheader(v2_obs_path)


def test_time_of_last_obs_does_not_overwrite_first():
    # Inline header where both lines are present — in v3.x the LAST line
    # would overwrite TIME_OF_FIRST_OBS. Here they must be distinct.
    sample = (
        "     2.11           OBSERVATION DATA    G                   RINEX VERSION / TYPE\n"
        "  2024     1     1     0     0    0.0000000     GPS         TIME OF FIRST OBS\n"
        "  2024     1     2     0     0    0.0000000     GPS         TIME OF LAST OBS\n"
        "     1    L1                                                # / TYPES OF OBSERV\n"
        "                                                            END OF HEADER\n"
    )
    import tempfile, os
    with tempfile.NamedTemporaryFile("w", suffix=".24o", delete=False) as fp:
        fp.write(sample)
        path = fp.name
    try:
        header = reader.read_obs_header(path)
        assert header.TIME_OF_FIRST_OBS == ["2024-1-1-0-0-0.0000000"]
        assert header.TIME_OF_LAST_OBS == ["2024-1-2-0-0-0.0000000"]
    finally:
        os.unlink(path)
