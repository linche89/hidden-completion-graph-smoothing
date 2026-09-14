# Public-release audit

Status: **not ready to publish the current Git history as a public repository**.

This audit is deliberately separate from the mathematical and experimental
acceptance checks.  It records redistribution, attribution, identity, and
repository-history blockers; it is not legal advice.

## Blocking items

1. The tracked history contains 69 third-party PDFs: 68 under `literature/`
   and one paper at the repository root.  Their presence in a private research
   workspace does not establish permission to redistribute them publicly.
2. The project has no author-selected source-code or manuscript license.  A
   license must not be inferred from dependency or dataset licenses.
3. Git currently records `ChopperLin <choppers@126.com>`.  The P5 plan calls
   for a publication email and an amended baseline, but the intended address
   has not been supplied.
4. No Git remote is configured.  Creating `github.com/linche89/...` also needs
   the final repository name and an explicit publication decision.

Do not solve the PDF issue by deleting the files only in a new commit: they
would remain downloadable from prior public history.

## Recommended publication shape

Create a clean public export with a new root commit, rather than pushing this
private research history.  Include only:

- `hidden_completion/`, the experiment modules and configs, tests, and build
  scripts;
- manuscript sources, original TikZ sources, generated vector figures, and
  checked result records;
- reproducibility documentation, environment pins, project license(s), and
  third-party notices;
- either documented data download instructions or the specifically licensed
  data subset with its upstream notices and metadata.

Exclude `literature/`, the root third-party PDF, scratch material, local virtual
environments, LaTeX intermediates, and credentials.  Keep the private repository
as the internal provenance record.

## External-data obligations already identified

- PGLib-OPF v23.07 data is marked CC BY 4.0; its bundled software is MIT.  A
  public artifact must retain the data notice, state the version, credit the
  original case sources in file headers, link the license, and say that this
  project extracts topology and positive branch weights.
- SuiteSparse matrices are marked CC BY 4.0.  A public artifact must retain
  source metadata and any matrix-specific citations, cite Davis and Hu (2011),
  link the license, and state the conversion to an undirected sparsity graph.
- The installed `tikz-diagrams` helper was used as a local QA/build workflow.
  No helper source code or template asset was copied into this repository.

## Checks completed

- Paper text distinguishes exact hidden-completion certificates from external
  fixed-graph stress tests.
- The paper names PGLib-OPF v23.07 and cites both benchmark collections.
- Exact CPU and GPU requirement files and a deterministic paper build script
  are present.
- Generated paper figures are original repository artifacts and retain their
  editable TikZ/CSV sources.

## Author decisions still required

- publication email;
- public repository name;
- source-code license;
- manuscript/preprint license and target venue's sharing policy;
- whether licensed raw benchmark data is vendored or downloaded on demand.

After those decisions, build the clean export, inspect its complete file list,
run a credential scan, compile the paper from the export, and only then create
and push the public remote.
