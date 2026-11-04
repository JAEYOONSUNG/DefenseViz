# Release validation — 2.1.3

## Saved R-M evidence and complete publication output

Legacy native exports do not always contain executed retained decisions. The
optional `--rm-evidence` bridge now evaluates complete, source-bound saved evidence
with the unchanged DNMB policy. It preserves original calls, all candidate audit
records and missing assessments. This is a fresh evaluation of stored evidence,
followed by aggregation and rendering; it is not a detector or recognition-engine
rerun. Newer exports can continue to use their saved decisions directly.

The final source passed **772 automated tests** and a clean installed-wheel check
(`CLEAN_INSTALL_OK`), with source hashes unchanged during both checks. Tests cover
multipart input, saved classification decisions, identity and scope conflicts,
retained context, partial assessment, explicit family mappings and figure layout.
Windows CI exposed text-mode hashing and default-codec mistakes in two test
fixtures. The fixtures now hash actual file bytes and read the UTF-8 layout
ledger explicitly; new CRLF cases also reject newline-only artifact tampering.
Production hashing and figure code were unchanged by these fixture corrections.
Remote CI is authoritative for the operating-system/Python matrix.

Three real native sample folders, including a separate external comparator, were
freshly imported with their own bound evidence. Their retained counts are 11, 9
and 9; adding separate reviewed annotations gives type-table totals 14, 11 and
12. Missing external functional and orphan assessments remain unassessed. No
primary-cohort review was borrowed for the external sample.

Historical **51-genome** and current **52-genome** inputs were imported and checked
separately. Independent source reconciliation passed **7,790** and **7,930** checks,
respectively. Integrated unions remain 2,886 / 2,929; effector counts 491 / 497;
annotation-based priority counts 494 / 500. Direct comparison against the original
figure source also verified all 676 R-M role/context cells and 135 type-case
numeric cells. These checks assess software and source reconciliation, not
independent biological accuracy.

Each cohort produces three pages: the detailed workflow/landscape a/b, R-M type
tables, and a Letter portrait comparison. All 57 category/family rows fit the
full-cohort comparison page at unchanged 7.8-point body size. Bars have equal
track dimensions, Total/Effector counts occupy separate right-aligned columns,
and observed PT partners are shown without inventing absent families. PDF and
editable SVG come from the same render. Only synthetic inputs are distributed.

Type-table assessment is conditional on the observed assigned-type pool; broader
candidate assessment remains separately exported. Untyped no-hit nominees remain
in the evaluation audit without inflating the R-M display pool. Partial counts
and unavailable assessments remain distinct from zero. Actual Claude critique
and its resolutions are recorded in [REVIEW.md](REVIEW.md).

---

# Release validation — 2.1.2

## Generic annotations and complete small-cohort figures

The default atlas has no target-strain or transformability legend or phenotype
stripe. Optional sample markers require an explicit field, label and color;
only literal `true` values are marked. Highlight keys, taxonomy keys and absence
keys describe items actually shown on the current page. Distinct family pools
with colliding display labels fall back to their declared IDs without changing
memberships, counts or analytical output.

Two and three real DNMBsuite sample folders were imported independently through
`from-suite`. Each produces three pages: the original five-card a / detailed b
layout, R-M type tables, and the original/prioritized overlaps with category/family
comparisons. Imported memberships match the prior saved-source import and all
consumed input hashes remain unchanged. These older exports contain only
unassessed retained-decision records; the corresponding R-M cells remain `NA`.
No later reviewed decisions are invented or copied from the 52-genome cohort.

The separate reviewed 52-genome replay produces four pages, including both
Letter comparison pages. All 57 comparison rows and 30 before/after intersection
regions match its saved analysis; integrated CDSs (2,929), effector CDSs (497)
and priorities (500) are unchanged. Only synthetic examples are distributed.

Independent review caught and verified the small-cohort label-collision fix.
**481 automated tests passed** on the final local source with unchanged source
hashes during execution. Installation and input receipts retain exact commands
and hashes; remote CI is authoritative for the current platform matrix.

---

# Release validation — 2.1.1

## Original manuscript a/b fidelity

The default saved-suite atlas now preserves the original five-card vector panel a and detailed panel b. The earlier compact publication view remains an explicit option. This is a saved-source display replay and aggregation, not a detector rerun or new biological validation.

Independent reference comparison verifies the five cards, 48 workflow labels, 16 icons, 18 arrows, 14 landscape blocks, 13 R-M columns and 2,496 family cells in the private reviewed cohort. Block and glyph coordinates/radii agree within 0.00005 mm (the reference ledger rounds coordinates). Counts, category bands, purple effector areas, six pink highlights, source membership symbols, taxonomy typography and legends match. Mini-bar tracks reserve measured digit space to prevent overlap; the type table remains a separate page rather than reproducing the original combined pagination.

The private source-bound display projection passes 16,094 comparisons against the original saved arrays without changing gene-level assignments or memberships. Later reviewed evidence is still required to reproduce later reviewed annotations; raw older native exports do not acquire invented reviews.

## General inputs and verification

- **461 automated tests passed** on the final local source, with source hashes unchanged throughout execution.
- A wheel installed in a fresh environment passes both CLI workflows outside the checkout (`CLEAN_INSTALL_OK`).
- All **seven real native folders** import and render with consumed input hashes unchanged.


Regression cases cover explicit excluded group members, unknown roles/types, HMM-only providers, missing versus partial endpoints, provider-specific general predictions, mixed native/reference identities, neutral bound display annotations, contradictory saved criteria, long names, small and large cohorts, partial zeros, outside-source memberships and taxonomy color consistency. Actual native folders, synthetic CLI execution and clean wheel installation were rerun. The original live text helper also passes on Matplotlib 3.7.5 and 3.11.2 in separate environments.

Exact original design review and actual read-only Claude code critique are separate checks. Review findings were reproduced with synthetic inputs before correction; private genomic inputs were not included in the Claude source snapshot. See [review resolutions](REVIEW.md). Remote CI results remain authoritative for the operating-system/Python matrix.

---

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
