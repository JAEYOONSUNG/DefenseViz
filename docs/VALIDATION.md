# Release validation — 2.0.0

## Local execution

- **116 automated tests passed** across analysis (45), native IO (45), rendering (19), and end-to-end pipeline/extension integration (7).
- Built wheel installed into a fresh virtual environment, then demo, validate and run executed from outside the source tree. Relocated inputs reproduce the same analysis. The check emits `CLEAN_INSTALL_OK` only after these assertions.
- Real, read-only native imports: 121 DefenseFinder files, 71 PADLOC files and one DefensePredictor table loaded without import failure. These are file counts, not independent genomes or method-performance estimates.
- Synthetic examples exercise missing native output, empty calls, unknown families, duplicate/renamed identifiers, corrupt hashes, inconsistent source sequences, coherent versus stitched domains, independent-site-map disagreement, idempotent enrichment, custom-policy labeling and output preservation.

## Historical regression (local data, not redistributed)

The new public API evaluates all **968 current evidence records** exactly like revision027, including every decision field. Explicit sample/gene mappings preserve a **971-CDS input pool**: 968 scored candidates plus three native C-control candidates outside the R/M/S score scope. Retained CDS = **346**; existing structural-support records = **244**. Renaming every sample and gene identifier leaves decisions and numerical summaries unchanged.

A full primary-cohort run with **52 genomes** also reproduces the native caller union (**2,857 CDS**) and integrated union (**2,902 CDS**). Distinct caller counts are DefenseFinder **1,378**, PADLOC **2,001**, and DefensePredictor **1,831**. These overlap and must not be summed. Its generic native-family view contains 70 current mapped labels and paginates across six pages under the supplied layout options; it is not a claim of 70 genomes or 70 independent defense classes.

This verifies software migration and arithmetic, not biological activity or independently validated predictive accuracy. Some original policy criteria have historical calibration provenance. No thresholds were fitted during this release. No detector or structure-prediction engine was rerun.

## Figure and critique review

Synthetic small/large layouts and the actual-cohort first/last pages were visually inspected. Layout tests independently check actual text bounds, header alignment, pagination, data lengths, fixed bar heights, missing-data labels and editable SVG text. [Independent critique findings](REVIEW.md) were reproduced and addressed before release. Exact upstream module identities are in [source_lineage.json](source_lineage.json).

Remote CI covers Python 3.10/3.12/3.14 on Linux and Python 3.12 on macOS/Windows; current status is reported by the repository Actions checks rather than a hardcoded badge. The Docker recipe is supplied as an optional wrapper; the release's executed installation checks use clean Python virtual environments, not a Docker build.
