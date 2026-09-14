# Figure QA for the TSP draft

Checked on 2026-09-14 in research-paper mode.

- The original seven and the three core-evidence standalone TikZ/PGFPlots
  sources in `paper/figures/tikz/` pass the static safety checker.
  `figure_style.tex` is a shared style fragment and was correctly excluded
  because it is not a standalone `tikzpicture`.
- The vector PDFs are included without rasterization. Their source PNG contact
  sheets and all nine rendered manuscript pages were inspected.
- Automated `title_band_collision` flags are false positives caused by panel
  titles, legends, and axis labels occupying their intended bands. No true
  text/plot overlap, clipping, or ambiguous arrow attachment was observed.
- The TSP layout enlarges `scenario_scaling.pdf` to improve single-column axis
  legibility. The theorem map, certificate comparison, budget curve,
  fixed-graph tradeoff, and external stress plots remain readable at their
  embedded sizes.
- `cpu_gpu_crossover.pdf` is retained in the separate supplement because its
  hardware-specific result is secondary to the theorem-matched certificate.
- The core-evidence figures were reviewed in `contact-sheet-core-v02.png`.
  Their remaining automated title-band warnings are false positives caused by
  the two separated group-plot panel titles; original-resolution inspection
  found no clipping or overlap. These figures replace the older sampled
  counterexample as the main experimental evidence because they use a fixed
  54-problem grid, exact reevaluation, and a 420-fit same-distribution sampling
  study.
