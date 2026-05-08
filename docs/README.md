# PyRINEX manual

LaTeX source for the user manual.

## Build

```bash
cd docs
pdflatex manual.tex
bibtex manual
pdflatex manual.tex
pdflatex manual.tex
```

Or, if you have `latexmk`:

```bash
cd docs
latexmk -pdf manual.tex
```

The build artefacts (`manual.aux`, `manual.bbl`, `manual.pdf`, etc.) are
ignored by `.gitignore` — only `manual.tex` and `manual.bib` are tracked.

## Notes

The legacy v3.x manual `PyRINEX manuel.pdf` is kept at the repository
root for archival purposes (it is the version cited in the published
PeerJ paper). New documentation should be added to `manual.tex`.
