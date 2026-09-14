# TSP submission checklist

## Manuscript

- [x] IEEEtran double-column, 10-point journal layout.
- [x] TSP-facing title, abstract, keywords, and introduction.
- [x] Abstract is within the 150--250 word range.
- [x] Related work and precise distinction are stated before the model.
- [x] One main theorem and three application corollaries.
- [x] State explicitly that a modest port boundary does not restrict the
  hidden-network size; the unbounded model permits arbitrarily many hidden
  vertices.
- [x] The complete hidden-input operator is retained.
- [x] Exact SDP and theorem-matched numerical evidence are included.
- [x] Fixed 54-problem grid, same-distribution sampling curves, and physical
  connected-path witnesses are included without result-based filtering.
- [x] Every approximate baseline is re-evaluated on the exact finite scenario
  set; solver-scale entries are thresholded and recertified.
- [x] Report descriptive end-to-end timings for all 54 theorem-matched exact
  finite designs without turning them into a hardware-independent claim.
- [x] All 14 primary `optimal_inaccurate` cells are independently rerun with
  high-accuracy SCS; all return `optimal` and agree in radius to (2.4\times
  10^{-8}) relatively.
- [x] Fixed-graph MSE is not presented as general WLS/Kalman MSE.
- [x] Fixed-known-graph, external-scale, numerical-audit, and CPU/GPU results
  are separated into a supplement.
- [ ] Replace anonymous author placeholders with names, affiliations, ORCIDs,
  funding acknowledgments, and corresponding-author contact.
- [ ] Recheck the final page count after author metadata is inserted.
- [ ] Select the final TSP EDICS classification in the submission system.

## Files and compliance

- [ ] Insert the public repository/archival DOI after the release record is
  frozen; verify every command from a clean environment.
- [x] Prepare a Git-history-free supplementary software archive with a
  SHA-256 manifest for review upload.
- [ ] Confirm licenses and redistribution permissions for every external data
  file and for the final PDF/preprint route.
- [ ] Upload editable LaTeX sources, bibliography, tables, and vector figures.
- [x] Verify that figure text remains legible at its final printed size.
- [ ] Run IEEE PDF/eXpress or the submission-system PDF checks if requested.
- [ ] Disclose any prior journal rejection and include the reports/response if
  the SPS resubmission policy applies; a mere internal fallback draft does not
  require disclosure.
- [ ] Add acknowledgments, conflicts, and any required data/code availability
  statement before submission.

## Final editorial pass

- [x] Replace working numerical precision with the final reporting precision.
- [x] Cross-check every percentage against the archived JSON/CSV record.
- [ ] Verify every bibliography entry against its publisher landing page.
- [x] Read every rendered page at 100% and inspect experiment figures in
  grayscale.
- [x] Cite and generate the supplement as a separate PDF.
- [ ] Upload the supplement as a separate submission file.
