# Use saved DNMBsuite results

DefenseViz reads saved DNMBsuite outputs and evaluates their evidence automatically.
It does not rerun DefenseFinder, PADLOC, DefensePredictor or structure prediction.
An explicitly configured current REBASE cache enables a fresh homology search
for the original reference nominees; its provenance remains separate from saved calls.
The public analysis pipeline name is **DNMB**; historical source paths and module
names remain in provenance.

## Install and run

```sh
python -m pip install '.[suite]'
defenseviz inspect-suite /path/to/DNMBsuite_outputs
defenseviz from-suite /path/to/DNMBsuite_outputs --out defenseviz-results
```

The input can be one genome's result directory or a parent containing nested
per-genome directories. Discovery walks subdirectories and stops at each native
DNMBsuite root, so its module and plasmid subfolders are not separate samples.
Hidden directories and verified DefenseViz result directories are excluded;
ordinary folders named `output` or `outputs` may still contain native inputs.
GenBank assembly metadata supplies stable sample identities, and source organism
and strain qualifiers supply display labels. Without assembly metadata the folder
basename remains the identity. Repeated distinct runs of the same identity require
an explicit manifest instead of silently selecting one. Original files are
read-only, and the output directory must be new. `inspect-suite` validates and
reports available endpoints without generating figures; it is not a required
preprocessing step for `from-suite`.

A [portable project file](suite-projects.md) can store these options for repeated single-command runs, including separately analyzed external type cases and combined Letter pages.

The full command automatically writes the classified workbook, retained R-M and
priority annotation lists, and the atlas, type tables and comparison figures.
PDF, editable SVG and PNG previews are generated together without a separate
plotting command. An optional layout file can disable PNG with `preview_png: false`.

R reads saved RDS completion/version metadata and, when requested, current
reference metadata through fixed packaged helpers. Set `--rscript /path/to/Rscript` when it is not on `PATH`.
Full mode sources the packaged fixed native core against your configured DNMB R installation. Stored evaluation uses fixed import helpers. Missing completion evidence must remain
explicit; a directory without results is not a successful zero-hit run.


## Full DNMB R-M execution

Configure local resources using [engine.example.json](engine.example.json). Paths may be absolute or relative to the configuration file. The R libraries must contain compatible DNMB and legacy DefenseViz runtime dependencies. The reference cache must contain `rebase_data.rds`, `rebase_db.fasta` and its BLAST indexes; profile metadata must describe the configured HMM library. Databases and installed runtime packages are read-only and are not distributed by this repository. BLAST+ must provide both `blastp` and `blastdbcmd` (the latter is found beside `blastp`, on PATH, or through an explicit configuration field). The index is checked against every FASTA entry before analysis.

```sh
defenseviz from-suite /path/to/DNMBsuite_outputs \
  --engine-config engine.json --evaluation full --out defenseviz-results
```

For subsequent runs, set `DEFENSEVIZ_ENGINE_CONFIG` or place the configuration at `~/.config/defenseviz/engine.json`. Default `auto` mode then selects full execution. Explicit `full` mode fails if resources are missing; it never silently falls back to stored evidence. `--rm-evidence` and `--reference-cache` belong to the separate stored-evidence path and cannot be combined with full mode.

Full mode searches all translated CDSs with the configured native R-M HMM library. It combines these nominations with native annotation nominees and qualifying saved DF/PAD profiles and exact source-bound previous R-M candidates, searches the complete configured REBASE database for nominated genes, evaluates motifs and Type I–IV context, and applies the unchanged selected policy. Saved comparator profiles nominate candidates; they do not substitute for native DNMB profile scoring. DF/PAD/DP outputs are preserved. This is a fresh R-M execution, not a rerun of those three tools or new structure prediction.

Reference type and component metadata are projected from explicit descriptions and numbered component prefixes. Repairs, search commands, exit codes, proteins, databases and outputs are hash-bound in the output evidence directory. Actual loaded R packages/DLLs and the frozen native source modules are verified. Previous candidate identities retain no inherited role, type or score credit. The original cache remains unchanged. Missing optional structure earns no support points; this is not a claim of biological inactivity. Retained annotations and component support are computational assessments.

Completed context inventories yield numeric MS/RMS/MR counts, including zero. An unavailable measurement is shown with an em dash. Qualified retained R component evidence is also propagated to putative effector counts and current-context support markers. Completed reviews report observed qualified annotation counts separately from unresolved biological functions; an observed zero does not establish biological absence. Included orphans and incomplete R contexts retain the compact `(n)` notation. When orphan or context status was not assessed, a positive count shows `(—)` instead, so it is not read as an assessed zero; Type IV R context is outside the Type I–III context rule and stays unmarked. Detailed coverage and source records remain in the tables and JSON; other uncertainty suffixes are not printed on the figure.

A failed full run preserves its logs in a separate `.failed-…` directory and produces no successful publication receipt. Each successful invocation requires a new output directory.

## Automatic evaluation

The default command collects source-bound candidate annotations, evaluates the
original DNMB motif and architecture rules on exact stored proteins, applies the
unchanged numerical policy, reassesses current retained contexts and reviews all
observed defense systems. It writes candidate, retained, not-retained, unresolved
and out-of-scope lists beside the original atlas and comparison figures.

Missing measurements stay unknown. `score_min` and `score_max` bound all
admissible completions of the missing evidence. A decision is confirmed only
when the entire interval gives the same retention outcome. This robustness
assessment is not experimental validation. A completed computational review
can still have unresolved functional evidence.

Some legacy exports lack reference lengths and the original DNMB profile search.
They cannot yield full-length homology support from identity alone. To evaluate
original reference nominees against a current DNMB reference cache in the same
invocation, provide its resource directory. This legacy partial-evaluation
example requires no configured full engine: neither `DEFENSEVIZ_ENGINE_CONFIG`
nor the default `~/.config/defenseviz/engine.json` (or its `XDG_CONFIG_HOME`
equivalent) may select one. With a configured engine, including the DNMBsuite
container, use the full command above and its configured references instead.

```sh
defenseviz from-suite /path/to/results --reference-cache /path/to/rebase/cache --out evaluated
```

This optional step requires BLAST+ `blastp` on `PATH` and reads the reference
FASTA plus metadata without changing them. The output records the current
database hashes, exact protein pairs, alignment command and search count.
All HSP intervals for each nominated pair are merged using the native coverage
and identity rules. A weak/no-hit result for one nominated subject cannot prove
that all homology evidence is absent. Original REBASE overlap memberships remain
unchanged. A current reference resource is never represented as the missing
historical database.

`--evaluation stored` disables automatic reevaluation for a saved-decision
replay. `--rm-evidence evidence.json` retains the existing explicit complete
measurement workflow. There is no need to construct that file for the automatic
path. Neither mode treats an absent annotation as a biological negative.

A below-cutoff R remains a gene in the genome: absence from the retained subset
does not create an orphan M/S annotation. Current retained-context composition
is stored separately from preserved native decisions.

System review coverage is calculated per displayed CDS set. A missing motif for
one family does not make an independently established effector in another family
unassessed. Exact multi-component catalytic cores count their distinct CDSs;
accessories require their own evidence.

## Select and label samples

For native outputs, default sample IDs use the assembly accession declared in
the source GenBank records; labels use their organism and strain qualifiers.
Explicit equivalent strain aliases are preserved in provenance and displayed
with one unambiguous source alias. Directory basenames are the fallback when
the corresponding source metadata is absent. A manifest can choose a cohort,
order rows, supply display labels or highlight strains:

```json
{
  "schema_version": 1,
  "samples": [
    {"sample_id": "isolate-a", "path": "runs/isolate-a", "label": "Example isolate A", "highlight": true},
    {"sample_id": "isolate-b", "path": "runs/isolate-b", "label": "Example isolate B"}
  ]
}
```

```sh
defenseviz from-suite . --manifest cohort.json --out selected-results
```

Manifest paths are relative to the manifest. Optional `cohort` metadata records
primary/external scope without adding automatic cohort labels to the figure.
The values `external` and `external_validation` select the external type-case
header color. Include a short scope label in `label` when readers need that role
explicitly stated; color alone does not explain the sampling design.
Legacy DefenseFinder gene exports without a `type` column can use their exact
`model_fqn` hierarchy for display taxonomy. Native fields remain unchanged and
the projection is recorded. This adds no retention or functional evidence.

Gene worksheets may contain constant HTTP(S) `HYPERLINK` labels; the importer
records each literal-label extraction with its cell and source hash. No formula
is evaluated. All other formulas remain unsupported.

IDs are opaque: the importer does not parse accessions, guess species, or merge
different loci because their proteins have the same sequence.

Optional sample `genus`, `species` and `display` metadata control taxonomy strips
and labels. For example, `display: {"taxon_label": "Example species",
"strain_label": "isolate A", "genus_color": "#537C83", "species_color": "#B29BCB"}`
keeps the taxon italic and the strain upright when these parts match `label`.
Update both `label` and its display parts when adding a suffix such as
`(external)`; inconsistent parts fall back to the literal label. No study-specific target or
transformability key is drawn by default. Sample highlights still use `highlight`
without assigning a biological meaning. An optional layout `highlight_label`
names the pink highlight only when highlighted samples are shown.

For an explicitly supplied binary sample annotation, a layout may set
`sample_marker: {"field": "reviewed", "label": "Reviewed sample", "color": "#526A75"}`.
This marks only samples with `display.reviewed: true`; the marker and its key
appear together. A stored `display.transformable` field has no special rendering
behavior: it requires an explicit `sample_marker` configuration like any other
annotation. Neither phenotype nor taxonomy is inferred from sample names.
Taxonomy and category keys are limited to items shown on the current page.

For nonstandard locations, sample-level `genes` and `fasta` paths can be supplied.
`modules` accepts `defensefinder`, `padloc`, `defensepredictor`, and
`rebasefinder` mappings with explicit artifact paths (`calls`, `id_map`,
`profiles`, `workbook`, `systems`, `fasta`) and a saved `status` when needed.
All such paths are relative to the manifest. Preserve the actual run status;
declaring completion is not a way to fill missing output with zeros.
For multipart DefenseFinder runs, `modules.defensefinder.calls` also accepts a
nonempty list of files sharing the sample's exact ID map. Each file and row
retains its provenance. Replicon-specific systems stay distinct and repeated
memberships do not inflate the unique CDS union.

DefensePredictor's saved `predicted_general_defense` flag is authoritative when
present. For score-only exports, declare `score_threshold` and `score_operator`
(`>` or `>=`), or point `decision_metadata` to the producer's JSON containing
`threshold` and `operator`. Conflicting explicit rules are rejected; a saved
flag/rule disagreement is recorded and rejected under `--strict-inputs`.
Legacy exports with neither flag nor rule use the historical inclusive `>= 4`
fallback, explicitly marked as an assumption in provenance. Missing saved flags
never become negative calls.

## What is preserved

Native DefenseFinder and PADLOC final memberships, general DefensePredictor
predictions, raw REBASE matches and saved DNMB decisions remain distinct.
Repeated system memberships do not inflate gene totals. HMM hits are evidence,
not final system members. Original native source counts never acquire DNMB-only
additions. Native system identifiers remain scoped to their sample and source.

Known source taxonomy names normalize to common display families, including
CAS/CBASS subtypes and numbered PADLOC PDC families. Exact native family names,
model IDs and system memberships remain in the audit. This mapping supplies no
retention or effector evidence. Ambiguous Argonaute exports remain untyped;
unknown family names are preserved. An explicit user `family_map` takes
precedence over built-in display aliases.

The importer supports version-dependent augmented REBASE tables and legacy
workbooks. Explicit saved `curation_keep` and final role/type fields are preserved
when available. Older exports without those fields enter automatic evidence evaluation; their
historical BLAST filter is not substituted for modern DNMB retention. Synthetic context labels never serve as REBASE reference identifiers. A
structured supplemental subject/identity pair can supply the reference match,
while its context role and type remain separate. Raw features,
motif observations, structures and annotation transfers are not interchangeable.

Automatic effector annotations use a versioned catalog of exact component models;
required complete complex cores are counted separately from individual effectors.
Raw InterProScan and EggNOG annotations are preserved for source-bound review.
Generic domains and a broad upstream `rease_motif_verified` flag do not establish
nuclease function: that flag can also describe motor-associated evidence.
Consequently, older outputs without explicit functional reviews do not reproduce
the expanded manuscript review automatically. Supply those saved reviews when
available; unsupported endpoints remain unassessed or observed lower bounds.

Candidate binding and reference-match completeness are separate endpoints.
An unbound no-hit candidate stays quarantined without implying a lost reference
match. An unbound or incomplete reference observation keeps reference coverage
partial. Broader historical candidates remain in the audit; an explicit current
evaluation scope is preserved in full. The displayed R-M input pool combines
original R-M annotations, reference hits and explicitly typed current candidates.
Untyped no-hit nominees remain evaluated/audited without inflating this denominator.

Current exported tables take precedence over cached tables. If live artifacts and
saved run metadata disagree, the conflict is reported rather than merged away.
Use `--strict-inputs` to reject these conflicts. A saved completion sentinel proves
run completion, not that an export was never subsequently edited.

## Restore R-M evaluation for legacy exports

An old REBASE workbook or augmented table may contain candidates but no final
`curation_keep` field. Importing that file alone cannot supply a retained count.
When complete DNMB evidence is available separately, evaluate it in the same
publication workflow:

```sh
defenseviz from-suite runs --manifest cohort.json \
  --rm-evidence rm-evidence.json --out evaluated-results
```

The evidence file binds native input hashes, exact sample/gene identities and
whole-protein SHA-256 values. Each record supplies all required evidence axes;
missing boolean assessments are errors, not false values. The selected policy
and per-CDS decisions are exported in `rm_evaluation.json`. This is fresh
evaluation of stored evidence, with zero detector runs. Original DF/PAD/DP and
REBASE memberships are preserved. Unassessed legacy candidates outside the
evaluated scope remain visible in the exported assessment coverage; partial
coverage remains in the audit; unavailable measurements use an em dash in figures.

Use the versioned evidence contract in [SUITE_CONTRACT.md](SUITE_CONTRACT.md).
It does not infer function from a raw homology match or synthesize missing
motif/context assessments. A contradictory existing saved decision requires an
explicit choice of the active input run (`decision_mode: "replace_scoped"` in
the evidence file). Omitted previous candidates remain unassessed and original
decisions are archived. Later reviewed Type IV components can
be supplied through the separate review supplement and stay distinguishable
from policy-retained CDSs.

## Include saved reviews

Optional reviews can supplement native exports when a later analysis supplied
additional component, context or effector evidence:

```sh
defenseviz from-suite runs --supplement reviews.json --out reviewed-results
```

A supplement has `schema_version: 1`, an `input_sha256` map binding actual imported
artifacts for every reviewed sample, and one or more canonical record arrays:
`annotations`, `effector_evidence`, `groups`, `complexes`. Every gene-specific
record binds its current `protein_sha256` and hashed `source_ref`/`source_refs`.
Group members require their protein hashes. Complex records additionally use
`core_protein_sha256: {"gene-id": "sha256"}` for every required core member.
Paths resolve relative to the supplement. Stale files, mismatched proteins and
contradictory assignments are errors. Supplements cannot add saved DNMB decisions
or native systems: those originate only in imported output. Review groups use
`membership_kind: "saved_context"` with `provider: "Review"`; the importer stamps
their review origin and run identity. Do not supply a run ID or member `retained`
flags. Reviews may contribute explicitly included reviewed annotations, which
remain separate from DNMB-retained CDSs.

See [the saved-result contract](SUITE_CONTRACT.md) for fields and set definitions.
Evidence completeness and activity remain separate. In particular, required
members of a complete restriction complex may count as effector CDSs without
assigning individual catalytic activity to each member.

## Output and interpretation

The Excel workbook is generated automatically by `from-suite`. Classified CDSs are the observed integrated defense union. Retained R-M annotations are the exact current retained subset; priority CDSs follow the exported annotation-based selection routes. These are separate sheets and TSVs because retention and prioritization are different sets. Numbers remain numeric Excel cells, filters and frozen headers are provided, and annotation text is stored literally rather than interpreted as formulas.


The output includes `tables/defense_annotations.xlsx`, `tables/classified_CDSs.tsv`, `tables/retained_RM_annotations.tsv`, `tables/prioritized_annotations.tsv`, an editable publication PDF, SVG pages, layout measurements,
`analysis.json`, original `native_bundle.json`, derived `analysis_bundle.json`,
tabular gene/priority/family/category results and `receipt.json`. The receipt pins inputs, implementation files,
available tool/database metadata, command, elapsed times and output hashes.
Unknown producer versions and absent database hashes are not invented.

`native_bundle.json` and `tables/imported_*` retain the original imported channels.
`analysis_bundle.json` and `tables/analysis_*` contain the evaluated inputs used for
the figures. `rm_evaluation.json` records every candidate's source scope, evidence,
unknown axes and decision; `tables/rm_candidate_evaluation.tsv` and the separate
retained/not-retained/unassessed/out-of-scope tables make those decisions inspectable.
`automatic_system_review.json` records component and effector review coverage.
Original final tool memberships are never augmented by later DNMB decisions.

Automatic nomination includes native DNMB candidates and original final R-M tool
members. The latter have their own `source_scope`; they are downstream integration
candidates, not a new independent detector. Existing assessed DNMB decisions are
preserved with their original provenance. New policy decisions and replayed
native decisions have separate counts and execution-stage labels in the receipt.
`policy_evaluations` counts attempted bounds evaluations; `policy_decisions` counts
newly resolved retention outcomes. The field-level `decision_transitions` and
`stored_vs_auto_summary` also expose assignment, evidence and coverage changes
when a record remains unassessed.
When `--rm-evidence` supplies a complete reviewed evidence file, the automatic
compiler/context/system stages are skipped explicitly to preserve that reviewed
workflow; this is recorded under `evaluation_execution`.

Input references are exported relative to the supplied root when possible.
Execution receipts deliberately retain exact executed command arguments and may
contain machine-specific paths, as may unchanged upstream metadata. This is
execution provenance, not a claim that every provenance field is relocatable.

The default **manuscript** atlas preserves the five-card vector workflow and the
detailed landscape: original R-M source mini-bars, the native-source intersection,
13 R-M role/context columns, retained R/total, noncanonical source memberships,
exclusive-membership symbols, family glyphs and taxonomy legends. The main atlas
orders native columns **PAD, DF, DP**; REB stays a separate reference channel.
Rings on R-M groups describe supplied support, independently of effector glyphs.
The separate reference column counts original reference hits, including untyped
subjects; it does not assign them an R-M type. The intersection excludes the
reference channel and never gains later DNMB-only additions. Unavailable original-context evidence uses an em dash.

In the R-M type tables, **Pool** counts CDSs placed in each type by DNMB or, without a
DNMB type, by a native R-M call. Under **Tool calls in pool (any family)**, DF, PAD and
DP count pool CDSs that each tool called in any defense family: this is CDS inclusion,
not agreement on the R-M type (PADLOC PDC-S or DMS_other calls, for example, are
counted). **DNMB** counts pool CDSs assigned to that type.

Category colors, purple evidence glyphs and optional pink sample highlights match
the manuscript grammar. Comparative tables separate putative nuclease-associated
**Total** and **Effector CDSs**, right-align numbers and order the default native
source columns **DF, PAD, DP**. Native and prioritized overlap views use the same
original source memberships. Larger datasets paginate; counts are never removed
to fit a page. SVG text remains editable; PDF also preserves vector artwork.

Use `--title "Defense comparison" --layout layout.json` to configure presentation:

```json
{
  "sections": ["atlas", "type_benchmarks", "comparison"],
  "font_family": "Arial",
  "font_fallbacks": ["Liberation Sans", "DejaVu Sans"],
  "atlas_style": "manuscript",
  "atlas_source_order": ["PADLOC", "DefenseFinder", "DefensePredictor"],
  "source_order": ["DefenseFinder", "PADLOC", "DefensePredictor"],
  "type_samples": ["isolate-a", "isolate-b"],
  "preview_png": true
}
```

Omit `type_samples` to select highlighted typed samples, or all typed samples
when none are highlighted. Explicitly list every desired comparison sample to
include unhighlighted external samples in that order; use `[]` to omit type tables.
`max_rows_per_page` and `max_families_per_page` limit atlas pages.
`atlas_page_size` sets maximum width/height in millimetres (default 506 × 365).
The manuscript template requires at least 447 × 200; data rows and family blocks
paginate without changing font sizes. Short inventories retain the circular section numbers and titles anchored to
their observed column groups. Measured title collisions move headings to a
second height tier and expand header space; missing families are not invented. Long strain labels expand
the page within its configured limit. Labels too long for that limit give an
actionable error. `atlas_style: "compact"` explicitly selects the earlier compact
view, whose minimum page limit is 270 × 145. `atlas_title` sets panel b's title.
`category_order` controls displayed category
order; an empty list preserves category encounter order while grouping partners.
`atlas_section_titles` overrides compact-view headings. Comparisons remain US Letter portrait and continue onto additional
pages. Additional sources or large intersection numbers use an exact intersection
table when a Venn diagram cannot remain legible. `source_colors`,
`category_colors`, and `aliases` are optional mappings. Unknown layout options
raise an error. The manifest records the actual font files and hashes; an
unavailable Arial installation is never reported as Arial.

Workflow thresholds are supplied presentation provenance, not new evaluation
parameters. `manuscript_workflow.parameters` accepts saved `retain_cutoff`,
`max_divergent_gap_bp` and `structure` values; omitted values use qualified saved
criteria wording. Result-bound values are authoritative; layout may fill a missing
value or repeat it, but conflicting values are rejected. Fresh source-bound R-M
evaluation supplies its actual policy metadata to this schematic.
`manuscript_workflow.pt_pairs` optionally lists displayed modification/restriction
partners; an explicit empty list disables them. Otherwise, observed families with
explicit partner display metadata supply up to two unambiguous examples. Absent
or ambiguous pairs are not invented, and omission reasons are recorded in the
layout. A layout setting never changes classification.

The native pipeline automatically attaches `analysis.manuscript_atlas`, a display
projection built from exact original memberships and saved groups. Python users
who already hold the canonical bundle and analyzed result can call:

```python
from defenseviz.plot.publication.manuscript_model import enrich_manuscript_model
from defenseviz.plot.publication import render_publication

enriched = enrich_manuscript_model(bundle, result)
render_publication(enriched, output_dir, sections=["atlas"])
```

Rendering an older result without that projection remains supported, but channels
that require the original bundle remain unavailable. See [display evidence](SUITE_CONTRACT.md#manuscript-display-projection)
for separately bound R-M support badges and cosmetic metadata.

These are saved annotation/evidence summaries. They do not establish biological
activity, prove a candidate is safe to remove, or supply independently adjudicated
precision/recall. Unreviewed results are not false positives or true negatives.
Family/category rows overlap, so their counts must not be summed as a genome-wide
gene union.
