# PyRINEX

[![DOI](https://zenodo.org/badge/718226106.svg)](https://zenodo.org/doi/10.5281/zenodo.10140408)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-Apache%202.0-green.svg)](LICENSE)

A Python package for batch processing RINEX 2/3 GNSS files with quality-check
products (multipath, ionospheric drift, cycle slips, sky plot, signal plot).

## Installation

```bash
pip install PyRINEX
```

Python 3.9+ is required. Runtime dependencies: `numpy`, `matplotlib`, `seaborn`,
`chardet`, `python-dateutil`.

## Features

| Module                       | What it does                                                  |
| ---------------------------- | ------------------------------------------------------------- |
| `PyRINEX.reader`             | Parse RINEX observation / navigation files into Python dicts  |
| `PyRINEX.data_management`    | Bulk discovery, marker / receiver / antenna cleaning          |
| `PyRINEX.quality_check`      | Multipath, ionosphere, cycle slip, sky plot, signal plot      |
| `PyRINEX.coordinates`        | WGS-84 ECEF ↔ BLH                                              |
| `PyRINEX.orbit`              | GPS broadcast-ephemeris orbit propagation                     |

## Quick start

```python
from PyRINEX import reader, quality_check

header = reader.read_obs_header("test0010.24o")
print(header.MARKER_NAME, header.APPROX_POSITION_XYZ)

quality_check.quality_check("test0010.24o")
```

The full API is documented in [`docs/manual.tex`](docs/manual.tex). To build the PDF:

```bash
cd docs
pdflatex manual.tex && bibtex manual && pdflatex manual.tex && pdflatex manual.tex
```

## Migrating from 3.x

PyRINEX 4.0 corrects the misspelled module name (`QulityCheck` → `quality_check`)
and stops returning JSON-encoded strings from the parser functions. Old import
paths and old function names still work and emit `DeprecationWarning`:

| 3.x                                               | 4.x                                                |
| ------------------------------------------------- | -------------------------------------------------- |
| `from PyRINEX.QulityCheck import QualityCheck`    | `from PyRINEX.quality_check import quality_check`  |
| `from PyRINEX.DataManagement import DataCleaning` | `from PyRINEX.data_management import clean_rinex`  |
| `json.loads(reader.oheader(path))`                | `reader.read_obs_header(path)`                     |
| `json.loads(reader.observations(path))`           | `reader.read_obs(path)`                            |
| `json.loads(reader.navigations(path))`            | `reader.read_nav(path)`                            |

See `docs/manual.tex`, §Changelog for the full list of fixes.

## Citation

If you use PyRINEX in your research, please cite:

> Han J, Lee SJ, Yun HS, Kim KB, Bae SW. **PyRINEX: a new multi-purpose Python package for GNSS RINEX data.** *PeerJ Computer Science* 10:e1800. 2024. https://doi.org/10.7717/peerj-cs.1800

```bibtex
@article{han2024pyrinex,
  title   = {PyRINEX: a new multi-purpose Python package for GNSS RINEX data},
  author  = {Han, Jinzhen and Lee, Seung Jun and Yun, Hong Sik and Kim, Kwang Bae and Bae, Sang Won},
  journal = {PeerJ Computer Science},
  volume  = {10},
  pages   = {e1800},
  year    = {2024},
  publisher = {PeerJ Inc.}
}
```

## License

Apache License 2.0
