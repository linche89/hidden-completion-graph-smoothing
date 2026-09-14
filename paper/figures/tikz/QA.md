# TikZ figure QA record

- Target medium: mathematical paper.
- Presentation mode: `research`.
- Language: English labels, axes, legends and captions.
- Render: XeLaTeX; one-page vector PDF; Poppler PNG at 200 dpi.
- Palette: Okabe--Ito colours with marker and dash redundancy.
- Original seven-figure batch review: `contact-sheet-v03.png`.
- Core-evidence batch review: `contact-sheet-core-v02.png`.

## Gates

All ten paper `.tex` files passed `check_tikz_safety.py`.  Every PDF compiled, every
PNG was nonblank, and the final contact sheet was inspected at original
resolution.  No clipping, text-on-text overlap, unreadable label, legend/data
collision or inconsistent reading order was found.  Targeted text-on-plot QA
for the sampling-gap label and the `CPU = GPU` label reported no
`text_plot_overlap` finding.

The bundled rendered checker reports `title_band_collision` for all seven
standalone figures.  Manual inspection classifies these as heuristic false
positives: the checker selects a panel title, legend entry or rotated y-axis
label as a page title and then treats other legitimate text on the same
standalone page as an intrusion.  The raw JSON reports are retained rather
than rewritten.  They contain no genuine pairwise text-overlap or edge-clipping
finding.

The three core-evidence additions (`core_design_benchmark`,
`sampling_budget_curve`, and `path_tightness`) also pass the static checker,
compile as one-page vector PDFs, and were manually inspected together in
`contact-sheet-core-v02.png`. Two focused legend repairs were made before the
`v02` snapshots. The remaining automated `title_band_collision` item is the
same false positive: the checker interprets the two deliberately separated
group-plot panel titles as a collision. Manual review at original resolution
finds no overlap, clipping, legend/data ambiguity, or unreadable type.

## Logic and design review

| Figure | Math/diagram logic | Complexity decision | Checked invariant |
|---|---|---|---|
| `theorem_structure` | schematic geometry; exact dependencies | keep | one main theorem and three numbering-independent consequence blocks; upper and converse proof inputs point into the main theorem |
| `budget_pareto` | exact | keep | plotted risk gaps match (q/(q+h)) to numerical tolerance; enumerated and pruned counts share the same certified radii |
| `certificate_sampling_failure` | exact | keep | all bar heights and the 0.295 gap are read from the canonical robust-result record |
| `local_method_tradeoff` | exact plotted measurements | keep | four fixed-graph families are averaged at common round values; zero-communication points are omitted only from the logarithmic byte-hop panel |
| `scenario_scaling` | exact | keep | finite candidates, retained vertices and partial-partition counts use (h=2) and (q=2,\ldots,8) |
| `cpu_gpu_crossover` | exact plotted measurements | keep | medians and min/max band use the same four grid sizes in all three process runs |
| `external_stress` | exact plotted measurements | keep | PGLib and SuiteSparse transformations remain separate and no WLS/Kalman semantic claim is added |
| `core_design_benchmark` | exact theorem-matched optimization results | keep | all 54 fixed grid cells are reported; every sampled and PSD rule is reevaluated on the same exact finite scenario set before comparison |
| `sampling_budget_curve` | exact post-evaluation of randomized designs | keep | both cases use nested, same-distribution training samples for each seed; all 420 fitted rules are evaluated on the complete finite scenario set |
| `path_tightness` | exact canonical risks and constructed path completions | keep | each path is connected, has maximum degree two, contains exactly the stated hidden budget, and the final risk/block gaps are below `3e-6` |

## Commands used

```powershell
$tikzSkill = '<path-to-installed-tikz-diagrams-skill>'
python "$tikzSkill/scripts/check_tikz_safety.py" <seven .tex files>
python "$tikzSkill/scripts/compile_render.py" figure.tex --engine xelatex --dpi 200 --visual-check --visual-mode research
python "$tikzSkill/scripts/make_contact_sheet.py" . --output contact-sheet-v03.png --cols 2 --thumb-width 760 --thumb-height 500 --show-qa
```
