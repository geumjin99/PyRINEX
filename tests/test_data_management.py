"""Tests for the file discovery / cleaning utilities."""
from __future__ import annotations

import os
from pathlib import Path

import pytest

from PyRINEX.data_management import (
    filter_by_extension,
    find_rinex,
    unique_filename,
)


def test_filter_by_extension_handles_subdirs(tmp_path: Path):
    (tmp_path / "a.txt").write_text("a")
    nested = tmp_path / "level1" / "level2"
    nested.mkdir(parents=True)
    (nested / "deep.txt").write_text("deep")
    (nested / "skip.csv").write_text("skip")

    matches = filter_by_extension(tmp_path, ".txt")
    matches_set = {os.path.basename(p) for p in matches}
    assert matches_set == {"a.txt", "deep.txt"}


def test_filter_by_extension_is_case_insensitive(tmp_path: Path):
    (tmp_path / "x.OBS").write_text("a")
    (tmp_path / "y.obs").write_text("b")
    matches = filter_by_extension(tmp_path, ".obs")
    assert len(matches) == 2


def test_find_rinex_keyword_filter(tmp_path: Path):
    a = tmp_path / "city_seoul" / "data"
    b = tmp_path / "city_busan" / "data"
    a.mkdir(parents=True)
    b.mkdir(parents=True)
    (a / "rover.24o").write_text("a")
    (b / "rover.24o").write_text("b")

    only_seoul = find_rinex(tmp_path, ["seoul"], ".24o")
    assert len(only_seoul) == 1
    assert "seoul" in only_seoul[0]


def test_unique_filename_appends_counter(tmp_path: Path):
    target = tmp_path / "report.csv"
    target.write_text("x")
    result = unique_filename(target)
    assert result.endswith("report_1.csv")
