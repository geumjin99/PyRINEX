# PyRINEX

[![DOI](https://zenodo.org/badge/718226106.svg)](https://zenodo.org/doi/10.5281/zenodo.10140408)
[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-Apache%202.0-green.svg)](LICENSE)

A Python package for batch processing RINEX 2.0/3.0 files with data quality analysis.

## Installation

```bash
pip install PyRINEX
```

## Features

| Module | Functions |
|--------|-----------|
| **Reader** | Parse observation/navigation files |
| **Data Management** | Batch rename, filter, coordinate conversion |
| **Quality Check** | Multipath, ionosphere, cycle slip analysis |

## Quick Start

```python
from PyRINEX.reader import oheader, observations
from PyRINEX.QulityCheck import QualityCheck

# Read RINEX header
header = oheader('example.24o')

# Run quality check
QualityCheck('example.24o')
```

## Documentation

📖 See [PyRINEX manuel.pdf](PyRINEX%20manuel.pdf) for detailed usage.

## Citation

If you use PyRINEX in your research, please cite:

> Han J, Kim T, Song J, Kim S, Oh K. **PyRINEX: a new multi-purpose Python package for GNSS RINEX data.** *PeerJ Computer Science* 10:e1800. 2024. https://doi.org/10.7717/peerj-cs.1800

```bibtex
@article{han2024pyrinex,
  title={PyRINEX: a new multi-purpose Python package for GNSS RINEX data},
  author={Han, Jinzhen and Kim, Taehee and Song, Jungho and Kim, Sanghoon and Oh, Keunyeong},
  journal={PeerJ Computer Science},
  volume={10},
  pages={e1800},
  year={2024},
  doi={10.7717/peerj-cs.1800}
}
```

## License

Apache License 2.0
