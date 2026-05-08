"""RINEX 2.xx and 3.xx observation/navigation file readers.

Returns plain Python objects (no JSON round-trip). The original 3.x
behaviour — functions returning ``json.dumps(...)`` strings — is preserved
in :mod:`PyRINEX.QulityCheck` and :func:`oheader`/``observations``/
``navigations`` shims for callers that still call ``json.loads(...)`` on
the output.
"""
from __future__ import annotations

import json
import math
import os
import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import chardet

PathLike = os.PathLike[str] | str
_NAV_KEYS_BASE = (
    "EPOCH",
    "SV clock bias",
    "SV clock drift",
    "SV clock drift rate",
    "IODE Issue of Data, Ephemeris",
    "Crs",
    "Delta n",
    "M0",
    "Cuc",
    "e Eccentricity",
    "Cus",
    "sqrt(A)",
    "Toe Time of Ephemeris",
    "Cic",
    "OMEGA0",
    "CIS",
    "i0",
    "Crc",
    "omega",
    "OMEGA DOT",
    "IDOT",
    "Codes on L2 channel",
    "GPS Week #",
    "L2 P data flag",
    "SV accuracy",
    "SV health",
    "TGD",
    "IODC Issue of Data, Clock",
    "Transmission time of message",
    "Fit interval",
)
NAV_KEYS_V2 = _NAV_KEYS_BASE
NAV_KEYS_V3 = _NAV_KEYS_BASE + ("spare1", "spare2")

# A single RINEX-2 observation field is exactly 16 chars wide; the
# original code used this magic string as a "missing observation"
# sentinel. We keep the convention because downstream tools depend on it.
EMPTY_OBS = "         0.000  "
_BLANK_FIELD = "                "


# ---------------------------------------------------------------------------
# Low-level file reading
# ---------------------------------------------------------------------------

def _detect_encoding(path: PathLike) -> str:
    with open(path, "rb") as fp:
        return chardet.detect(fp.read(200000))["encoding"] or "utf-8"


def read_lines(path: PathLike) -> List[str]:
    """Read a text file as a list of lines, falling back to chardet on
    decode errors. Always preserves trailing newlines."""
    try:
        with open(path, "r") as fp:
            return fp.readlines()
    except UnicodeDecodeError:
        with open(path, "r", encoding=_detect_encoding(path)) as fp:
            return fp.readlines()


def list_files(path: PathLike) -> List[str]:
    """Recursively list every regular file under ``path``."""
    out: List[str] = []
    for root, _dirs, files in os.walk(path):
        out.extend(os.path.join(root, f) for f in files)
    return out


# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------

@dataclass
class ObsHeader:
    """Parsed RINEX observation-file header.

    The ``[line_number, value]`` shape used in v3.x of this library is
    preserved on a few fields (``MARKER_NAME``, ``MARKER_NUMBER``,
    ``receiver_type``, ``antenna_type``) because :mod:`data_management`
    rewrites those lines in place.
    """

    version: int  # 2 or 3
    type: str  # observation type letter, usually 'G', 'M', etc.
    APPROX_POSITION_XYZ: List[float] = field(default_factory=list)
    PRNS: List[str] = field(default_factory=list)
    ObsTypes: Any = field(default_factory=dict)  # list (v2) or dict (v3)
    END_OF_HEADER: int = -1
    MARKER_NAME: Optional[List[Any]] = None
    MARKER_NUMBER: Optional[List[Any]] = None
    receiver_type: Optional[List[Any]] = None
    antenna_type: Optional[List[Any]] = None
    TIME_OF_FIRST_OBS: Optional[List[str]] = None
    TIME_OF_LAST_OBS: Optional[List[str]] = None

    def to_dict(self) -> Dict[str, Any]:
        """Return the header as a plain dict mirroring the legacy keys."""
        # Legacy code used spaces in keys; preserve that for the JSON shim.
        return {
            "version": str(self.version),
            "type": self.type,
            "MARKER_NAME": self.MARKER_NAME,
            "MARKER_NUMBER": self.MARKER_NUMBER,
            "receiver_type": self.receiver_type,
            "antenna_type": self.antenna_type,
            "APPROX POSITION XYZ": self.APPROX_POSITION_XYZ,
            "TIME OF FIRST OBS": self.TIME_OF_FIRST_OBS,
            "TIME OF LAST OBS": self.TIME_OF_LAST_OBS,
            "END OF HEADER": self.END_OF_HEADER,
            "PRNS": self.PRNS,
            "ObsTypes": self.ObsTypes,
        }


def _normalize_prn(token: str) -> str:
    """Pad a 3-char PRN like ' G1' / 'G 1' to canonical 'G01' form."""
    if len(token) < 3:
        token = token.rjust(3)
    if token[0] == " ":
        token = "G" + token[1:]
    if token[1] == " ":
        token = token[0] + "0" + token[2]
    return token


def _scan_v2_prns(lines: List[str], start: int, year_suffix: str) -> List[str]:
    """Walk every epoch in a RINEX-2 obs body and collect PRN strings."""
    seen: List[str] = []
    n = start
    while n < len(lines):
        if year_suffix == lines[n][1:3] and lines[n][3:4] == " ":
            try:
                sat_num = int(lines[n][29:32])
            except ValueError:
                n += 1
                continue
            epoch_lines = math.ceil(sat_num / 12)
            sat_block = ""
            for m in range(n, min(n + epoch_lines, len(lines))):
                sat_block += lines[m][32:68].rstrip("\n")
            for token in re.findall(r".{3}", sat_block):
                prn = _normalize_prn(token)
                if prn.strip() and prn not in seen:
                    seen.append(prn)
            n += epoch_lines
        else:
            n += 1
    seen.sort()
    return seen


def _scan_v3_prns(lines: List[str], start: int) -> List[str]:
    """Walk every epoch in a RINEX-3 obs body and collect PRN strings."""
    seen: List[str] = []
    n = start
    while n < len(lines):
        if lines[n].startswith(">"):
            try:
                sat_num = int(lines[n][32:35])
            except ValueError:
                n += 1
                continue
            for m in range(n + 1, min(n + 1 + sat_num, len(lines))):
                prn = lines[m][:3]
                if prn.strip() and prn not in seen:
                    seen.append(prn)
            n += 1 + sat_num
        else:
            n += 1
    seen.sort()
    return seen


def read_obs_header(path: PathLike) -> ObsHeader:
    """Parse the header section of a RINEX observation file."""
    lines = read_lines(path)
    if not lines:
        raise ValueError(f"empty file: {path}")

    version_char = lines[0][5]
    type_char = lines[0][40] if len(lines[0]) > 40 else " "
    header = ObsHeader(version=int(version_char), type=type_char)

    obstypes_line_numbers: List[int] = []
    prns: List[str] = []

    for n, line in enumerate(lines):
        if "MARKER NAME" in line:
            header.MARKER_NAME = [n, line.split()[0] if line.split() else ""]
            if n + 1 < len(lines):
                header.MARKER_NUMBER = [
                    n + 1,
                    lines[n + 1].split()[0] if lines[n + 1].split() else "",
                ]
        elif "REC # / TYPE / VERS" in line:
            header.receiver_type = [n, line[20:40]]
        elif "ANT # / TYPE" in line:
            header.antenna_type = [n, line[20:40]]
        elif "APPROX POSITION XYZ" in line:
            tokens = line.split()
            header.APPROX_POSITION_XYZ = [
                float(tokens[0]),
                float(tokens[1]),
                float(tokens[2]),
            ]
        elif "TIME OF FIRST OBS" in line:
            header.TIME_OF_FIRST_OBS = ["-".join(line.split()[0:6])]
        elif "TIME OF LAST OBS" in line:
            header.TIME_OF_LAST_OBS = ["-".join(line.split()[0:6])]
        elif "# / TYPES OF OBSERV" in line or "SYS / # / OBS TYPES" in line:
            obstypes_line_numbers.append(n)
        elif "PRN / # OF OBS" in line:
            prn = line[3:6]
            if prn.strip():
                prns.append(prn)
        elif "END OF HEADER" in line:
            header.END_OF_HEADER = n

    if header.END_OF_HEADER < 0:
        raise ValueError(f"END OF HEADER not found in {path}")

    # Fall back to scanning the body if the header didn't list PRNs.
    if not prns:
        year_suffix = "".join(str(path)[-3:-1])
        if header.version == 2:
            prns = _scan_v2_prns(lines, header.END_OF_HEADER + 1, year_suffix)
        elif header.version == 3:
            prns = _scan_v3_prns(lines, header.END_OF_HEADER + 1)

    header.PRNS = [_normalize_prn(p) for p in prns]

    if obstypes_line_numbers:
        first_marker = lines[obstypes_line_numbers[0]]
        if "# / TYPES OF OBSERV" in first_marker:
            obs_types: List[str] = []
            for ln in obstypes_line_numbers:
                obs_types.extend(lines[ln][6:60].split())
            header.ObsTypes = obs_types
        else:  # SYS / # / OBS TYPES (v3)
            obs_types_per_sys: Dict[str, List[str]] = {}
            current_sys: Optional[str] = None
            for ln in obstypes_line_numbers:
                line = lines[ln]
                lead = line[0]
                if lead != " ":
                    current_sys = lead
                    obs_types_per_sys[current_sys] = []
                if current_sys is not None:
                    obs_types_per_sys[current_sys].extend(line[6:60].split())
            header.ObsTypes = obs_types_per_sys

    return header


# ---------------------------------------------------------------------------
# Observations
# ---------------------------------------------------------------------------

def _pad_field(field_str: str, width: int = 16) -> str:
    """Right-pad a possibly-truncated 16-char observation field."""
    if len(field_str) < width:
        return field_str + " " * (width - len(field_str))
    return field_str


def _parse_v2_obs_block(lines: List[str], start: int, n_obs: int) -> List[str]:
    """Concatenate ``n_obs`` 16-char fields from the consecutive RINEX-2
    observation lines beginning at index ``start`` (5 fields per line)."""
    fields_per_line = 5
    raw = ""
    for i in range(math.ceil(n_obs / fields_per_line)):
        line = lines[start + i][:80].rstrip("\n")
        raw += _pad_field(line, 80)
    fields = re.findall(r".{16}", raw)
    while len(fields) < n_obs:
        fields.append(_BLANK_FIELD)
    return [EMPTY_OBS if f == _BLANK_FIELD else f for f in fields[:n_obs]]


def _parse_v3_obs_line(line: str, n_obs: int) -> List[str]:
    """Slice the data portion (column 3 onward) of one RINEX-3 sat line
    into ``n_obs`` 16-character fields."""
    data = line[3:].rstrip("\n")
    needed = n_obs * 16
    if len(data) < needed:
        data = data + " " * (needed - len(data))
    fields = re.findall(r".{16}", data)
    return [EMPTY_OBS if f == _BLANK_FIELD else f for f in fields[:n_obs]]


def read_obs(path: PathLike) -> Dict[str, Dict[str, Any]]:
    """Parse a RINEX observation file body, keyed by epoch label.

    The epoch label format is preserved from earlier releases: a single
    string of ``"yy mm dd hh mi ss.sssssss"`` for v2 and the same for v3
    (with the leading ``>`` stripped). Each value is a dict containing
    ``sat_num`` and one entry per PRN whose value is a dict from
    observation-type to 16-char string.
    """
    header = read_obs_header(path)
    lines = read_lines(path)
    epochs: Dict[str, Dict[str, Any]] = {}
    end_of_header = header.END_OF_HEADER

    if header.version == 2:
        obs_types = header.ObsTypes  # type: ignore[assignment]
        if not isinstance(obs_types, list):
            raise ValueError("v2 ObsTypes must be a list")
        sat_lines = math.ceil(len(obs_types) / 5)
        year_suffix = "".join(str(path)[-3:-1])

        n = end_of_header + 1
        while n < len(lines):
            line = lines[n]
            if year_suffix != line[1:3] or line[3:4] != " ":
                n += 1
                continue
            try:
                sat_num = int(line[29:32])
            except ValueError:
                n += 1
                continue

            sat_lines_count = math.ceil(sat_num / 12)
            sat_block = ""
            for m in range(n, n + sat_lines_count):
                sat_block += lines[m][32:68].rstrip("\n")
            sats = [
                _normalize_prn(t)
                for t in re.findall(r".{3}", sat_block)
                if t.strip()
            ]

            epoch: Dict[str, Any] = {"sat_num": sat_num}
            for idx, prn in enumerate(sats):
                start = n + sat_lines_count + idx * sat_lines
                fields = _parse_v2_obs_block(lines, start, len(obs_types))
                epoch[prn] = dict(zip(obs_types, fields))

            label = " ".join(line.split()[0:6])
            epochs[label] = epoch
            n += sat_lines_count + sat_num * sat_lines

    elif header.version == 3:
        obs_types_per_sys = header.ObsTypes  # type: ignore[assignment]
        if not isinstance(obs_types_per_sys, dict):
            raise ValueError("v3 ObsTypes must be a dict")

        n = end_of_header + 1
        while n < len(lines):
            line = lines[n]
            if not line.startswith(">"):
                n += 1
                continue
            try:
                sat_num = int(line[32:35])
            except ValueError:
                n += 1
                continue

            epoch = {"sat_num": sat_num}
            for m in range(n + 1, n + 1 + sat_num):
                if m >= len(lines):
                    break
                sat_line = lines[m]
                system = sat_line[0]
                obs_types = obs_types_per_sys.get(system, [])
                fields = _parse_v3_obs_line(sat_line, len(obs_types))
                epoch[sat_line[0:3]] = dict(zip(obs_types, fields))

            label = " ".join(line.split()[1:7])
            epochs[label] = epoch
            n += 1 + sat_num
    else:
        raise ValueError(f"unsupported RINEX version: {header.version}")

    return epochs


# ---------------------------------------------------------------------------
# Navigation
# ---------------------------------------------------------------------------

def read_nav_iono(path: PathLike) -> List[str]:
    """Return the eight ionosphere model parameters (alpha + beta) from a
    GPS navigation file, or ``["NO ION data."]`` if absent.
    """
    contents = read_lines(path)
    if not contents:
        return ["NO ION data."]
    version = contents[0][5]
    iono: List[str] = []
    block = ""
    if version == "2":
        for n, line in enumerate(contents):
            if "ION ALPHA" in line:
                for m in range(n, n + 2):
                    block += contents[m].rstrip("\n")[2:50]
                iono = re.findall(r".{12}", block)
                break
    elif version == "3":
        for n, line in enumerate(contents):
            if "GPSA" in line:
                for m in range(n, n + 2):
                    block += contents[m].rstrip("\n")[5:53]
                iono = re.findall(r".{12}", block)
                break
    if not iono:
        return ["NO ION data."]
    return iono


def read_nav(path: PathLike) -> Dict[str, List[Dict[str, str]]]:
    """Parse a GPS broadcast-navigation file.

    Returns a mapping ``PRN → list-of-records``. Each record is a dict
    keyed by the human-readable parameter names defined in
    :data:`NAV_KEYS_V2` / :data:`NAV_KEYS_V3`. The 19-character
    Fortran-style field strings are returned verbatim — convert with
    :func:`PyRINEX.orbit._to_float` (or your own parser) when you need
    numbers.
    """
    contents = read_lines(path)
    if not contents:
        return {}

    version = contents[0][5]
    keys = NAV_KEYS_V2 if version == "2" else NAV_KEYS_V3

    end_of_header = -1
    for n, line in enumerate(contents):
        if "END OF HEADER" in line:
            end_of_header = n
            break
    if end_of_header < 0:
        raise ValueError(f"END OF HEADER not found in {path}")

    out: Dict[str, List[Dict[str, str]]] = {}
    m = end_of_header + 1
    while m < len(contents):
        line = contents[m]
        if version == "2":
            # v2 nav PRN is the first two characters: ' 1' or '23'.
            prn_token = line[0:2]
            if not prn_token.strip().isdigit():
                m += 1
                continue
            prn_num = int(prn_token.strip())
            prn = f"G{prn_num:02d}"
            data_offset = 3
        else:  # v3
            if line[0:3].strip() == "":
                m += 1
                continue
            prn = line[0:3]
            if prn[1] == " ":
                prn = prn[0] + "0" + prn[2]
            data_offset = 4

        if m + 7 >= len(contents):
            break

        block = "".join(contents[m + i][data_offset:].rstrip("\n") for i in range(8))
        fields = re.findall(r".{19}", block)
        record = dict(zip(keys, fields[: len(keys)]))
        out.setdefault(prn, []).append(record)
        m += 8

    return out


# ---------------------------------------------------------------------------
# Legacy JSON-string API
# ---------------------------------------------------------------------------

def oheader(path: PathLike) -> str:
    """Deprecated. Use :func:`read_obs_header` and operate on the dataclass.

    Kept for compatibility with the v3.x example code that calls
    ``json.loads(oheader(path))``.
    """
    warnings.warn(
        "oheader() returns a JSON string for compatibility; "
        "prefer read_obs_header() which returns an ObsHeader dataclass.",
        DeprecationWarning,
        stacklevel=2,
    )
    return json.dumps(read_obs_header(path).to_dict())


def observations(path: PathLike) -> str:
    """Deprecated JSON-string wrapper around :func:`read_obs`."""
    warnings.warn(
        "observations() returns a JSON string for compatibility; "
        "prefer read_obs() which returns a dict.",
        DeprecationWarning,
        stacklevel=2,
    )
    return json.dumps(read_obs(path))


def navigations(path: PathLike) -> str:
    """Deprecated JSON-string wrapper around :func:`read_nav`."""
    warnings.warn(
        "navigations() returns a JSON string for compatibility; "
        "prefer read_nav() which returns a dict.",
        DeprecationWarning,
        stacklevel=2,
    )
    return json.dumps(read_nav(path))


def navigaions(path: PathLike) -> str:  # noqa: D401 — matches manual.pdf typo
    """Alias preserved for the typo in the published manual."""
    return navigations(path)


def nheader(path: PathLike) -> List[str]:
    """Deprecated alias of :func:`read_nav_iono`."""
    warnings.warn(
        "nheader() is renamed to read_nav_iono() in 4.x.",
        DeprecationWarning,
        stacklevel=2,
    )
    return read_nav_iono(path)
