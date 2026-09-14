# Review-artifact QA

Checked on 2026-09-14 with `scripts/package_tsp_artifact.ps1`.

- Archive: `output/release/tsp_anonymous_artifact.zip`
- Packaged files: 159
- Archive size: 8,032,547 bytes
- SHA-256: `071DE3C5A891AC923A76876692BED0CE3CDBA60D1813068057DF684955B73B38`
- Every entry in the packaged `MANIFEST.sha256` was independently rehashed;
  zero mismatches were found.
- No `.git`, build directory, Python cache, LaTeX auxiliary, or bytecode file
  is present.
- A case-insensitive scan found no author name, GitHub handle, or user-profile
  path. The only email-shaped strings are two generic institutional contact
  templates retained verbatim in third-party PGLib data comments; neither
  identifies an author or local user.
- External PGLib and SuiteSparse inputs are packaged together with their
  retained license, source, and checksum notices.
