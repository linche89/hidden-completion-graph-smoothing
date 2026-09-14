# PDF QA for the TSP submission pair

Checked on 2026-09-14 after the complete `scripts/build_tsp.ps1` build.

- `tsp_manuscript.pdf`: 9 US-letter pages, a 242-word abstract (246 under
  TeXcount's text-plus-inline-math sum), one main theorem, and three
  corollaries.
- `tsp_supplement.pdf`: 3 US-letter pages.
- Both final LaTeX logs contain no unresolved reference, citation, package, or
  overfull-box warning. The main log has two `Underfull \\vbox` notices, on the
  top-matter page and a float-heavy experiment page; both were checked and
  contain no layout defect.
- All twelve pages were rendered with Poppler at 160 dpi and inspected at
  original resolution. No clipping, missing glyph, blank page, occluded text,
  or illegible embedded figure was found.
- The two experiment-heavy manuscript pages were also inspected in grayscale.
  Marker shapes and dash patterns preserve the curve distinctions without
  color.
- Extracted PDF text contains no unresolved `??` marker or unexpanded generated
  number macro.

Final SHA-256 values:

- manuscript: `C82810330EE46EF140F6C8E762A0669BC3DBCAED09D3A53D319148E0A6F41C2B`
- supplement: `D8488FEA42EDC976DDB8E8AA4CF3D742217C35C04275FC14684EFFB1500865D3`
