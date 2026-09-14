# TSP submission draft

This directory contains the manuscript and supplement prepared for
*IEEE Transactions on Signal Processing*.

The TSP version emphasizes the signal-processing problem: exact robust design
of a common local approximation to the rational graph filter `(I+L)^{-1}` when
topology and inputs beyond a modest visible port set are hidden. The hidden
network may be arbitrarily large; the exact algorithm is fixed-parameter in
the observed port count, not in the total physical-graph size. The mathematical
scope is frozen; edits in this directory concern positioning, exposition,
experiments, and submission layout.

## Build

From the repository root, the complete checked build is:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/build_tsp.ps1
```

This regenerates table/figure data, compiles the standalone vector figures,
runs the regression suite, and builds both PDFs. Use `-SkipTests` or
`-SkipData` only for an editing pass. The main manuscript is `main.tex`;
fixed-known-graph algorithms, external scale tests, numerical postsolve
details, and hardware-specific timing are in `supplement.tex`.

Current verified build: 9 manuscript pages, 3 supplementary pages, a
242-word abstract (246 including four inline-math units), 1 main theorem, and
3 corollaries. Both LaTeX logs are free
of unresolved references and overfull boxes; the main log has two benign
`Underfull \\vbox` notices on pages 1 and 8 that were checked visually.

## Current submission constraints

The working target follows the IEEE Signal Processing Society author guidance
checked on 2026-09-14:

- regular-paper initial submission: at most 13 double-column, 10-point,
  single-spaced pages, including figures, appendices, and references;
- abstract: 150--250 words;
- TSP uses single-anonymized review;
- novelty, experimental evidence, English quality, and adequate citation can
  be screened before external review.

Official sources:

- <https://signalprocessingsociety.org/publications-resources/information-authors>
- <https://signalprocessingsociety.org/publications-resources/ieee-transactions-signal-processing/about-ieee-transactions-signal-processing>
