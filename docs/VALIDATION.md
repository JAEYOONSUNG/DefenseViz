# Release validation — 2.1.0

## Saved DNMBsuite input

The `from-suite` command imports native saved output directories and exports an atlas, R-M type tables, original/prioritized source overlaps and category/family inventories. It aggregates saved annotations and explicit evidence; it does not rerun defense detectors, recognition algorithms or structure prediction.

Validation separates two data sources: seven read-only native folders spanning multiple saved-output versions (a mixed primary/external set), and a separate 52-genome source-plus-review regression. The latter includes later source-bound functional/context reviews; it is not evidence that old native files contain all those reviews. Only synthetic examples are distributed.

## Executed checks

- **319 automated tests passed** locally on Python 3.14, including native input, legacy compatibility, scientific set boundaries, plotting and complete CLI execution.
- A built wheel installed into a fresh virtual environment passes legacy and native demos, inspection, standalone execution and analytical equivalence outside the checkout (`CLEAN_INSTALL_OK`).
- All **seven real native folders** import and render through `from-suite`; consumed input hashes remain unchanged. Existing metadata/export conflicts remain explicit, not silently repaired. Tests also cover strict rejection.
- The separate **52-genome source-plus-review regression passes 17,965 checks**. Its integrated union (2,929), effector CDSs (497) and annotation-based priorities (500) remain identical to the pinned reviewed inputs.
- Root visually inspected native, synthetic and private-cohort outputs. Geometry checks cover text/track collisions, right-aligned numbers, fixed bar dimensions, exact intersections, type states, nested glyph areas, pagination and editable text.


The private regression checks exact membership and all family/category/cell/source intersections, typed-pool partitions and case context/orphan subsets. It does not tune thresholds or assess biological accuracy. Final gene-union and priority quantities remain equal to the pinned reviewed inputs. Sources and producer versions remain explicit; unknown historical producer versions are not backfilled from a newer checkout.

## Review and limitations

Independent input, analysis and rendering reviewers identified and corrected native ID/schema variants, contradictory duplicate context, false zero states, non-RM contamination of typed totals, review/native provenance confusion and small-cohort layout collisions. Actual Claude Opus review was used for the code and rendered figures; its material findings were reproduced and fixed before publication. See [review resolutions](REVIEW.md).

Automatic effector transfer uses exact annotated components and complete required cores. Raw motif/domain/orthology observations remain separate. Older exports without saved curation or detailed functional review remain unassessed or lower bounds. Context completeness does not demonstrate activity. No broad accuracy metrics or superiority claims are made.

SVG retains editable text; PDF preserves vector geometry. Arial is used when installed, with actual font hashes and explicit fallbacks. Letter comparison pages paginate without shrinking numbers or dropping rows. Remote CI runs the documented matrix; the Actions result is authoritative for its status.

---

# Historical validation — 2.0.0

## Local execution

- **119 automated tests passed** across analysis (45), native IO (46), rendering (20), and end-to-end pipeline/extension integration (8).
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
