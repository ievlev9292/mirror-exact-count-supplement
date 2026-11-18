# Moving-mirror figures and local theorem

Companion materials for *Particle number and energy from the shape of relativistic mirror motion* by Michael R. R. Good, Evgenii Ievlev, and Edige Zholdasbay.

## Figures

- `rapidity.pdf`: signed rapidity versus laboratory time for the polynomial flights and returning motions.
- `penrose.pdf`: the same trajectories in compactified spacetime.
- `shape_radiation.pdf`: the count/energy comparison retained as a commented-out figure in this manuscript version.

Each PDF has its editable TeX source beside it. All curves use a common time scale, with units chosen so that tau=1. Symmetric examples have maximum speed V=0.97; the asymmetric cubic return has V_+=0.97 and V_-=0.65. The elliptic round trip uses L/tau=5.175133604750678, determined from its actual maximum speed V=0.97. Its turning point is t=z=0, and its asymptotic position is z=L. The curves have no additional position or time shifts.

Sampling is deterministic. The script includes small normalization, trajectory and calibration checks. Its tests support numerical reproduction; the analytic argument is in the proof.


## Proof of the theorem for ultrarelativistic asymptotics

`proof/proof.tex` is the single source, with `references.bib` and the compiled `proof.pdf` beside it. It introduces notation, states the separated-event theorem, and then proves the single positive event, its number coefficient and energy, followed by negative and multiple events. It states the physical meaning of the characteristic laboratory-time scale and the uniform hypotheses needed for the result.

The positive-number and energy representations and the exact polynomial-flight count are cited as inputs from the accompanying manuscript. Their derivations are in that manuscript. All counts and energies include both exterior regions for one real scalar field.


## Technical note

Use Python 3.12 with the versions in `requirements.txt`, plus pdfLaTeX and BibTeX. The figure sources require TikZ, PGFPlots 1.18 and `standalone`; the proof uses standard AMS packages, `geometry`, `hyperref` and Latin Modern fonts.

From this folder, run:

```text
python reproduce.py
```

This rebuilds the three PDFs in `figures/` and `proof/proof.pdf`. Generated data and TeX intermediate files go into a temporary working folder whose location is printed. No external data downloads or TeX shell escape are required.

Use `--figures-only` or `--proof-only` to select an output group. If the TeX tools are not on your path, supply their executable paths with `--tex-engine` and `--bibtex-engine`. The optional `--work-dir` argument retains all intermediate files in a chosen folder.