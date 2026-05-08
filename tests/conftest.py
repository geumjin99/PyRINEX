"""Pytest fixtures shared across the test suite."""
from __future__ import annotations

from pathlib import Path

import pytest

from . import fixtures


@pytest.fixture()
def v2_obs_path(tmp_path: Path) -> Path:
    return fixtures.write_v2_obs(tmp_path)


@pytest.fixture()
def v2_nav_path(tmp_path: Path) -> Path:
    return fixtures.write_v2_nav(tmp_path)


@pytest.fixture()
def v3_obs_path(tmp_path: Path) -> Path:
    return fixtures.write_v3_obs(tmp_path)


@pytest.fixture()
def v2_pair(tmp_path: Path) -> tuple[Path, Path]:
    """A matching .24o / .24n pair under the same directory."""
    return fixtures.write_v2_obs(tmp_path), fixtures.write_v2_nav(tmp_path)
