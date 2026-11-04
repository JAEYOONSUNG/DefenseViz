# Saved DNMBsuite results contract

This additive interface preserves the existing normalized-project evaluation API.
`from-suite` imports saved results and aggregates explicit evidence; it never calls
the historical scoring predicate on incomplete native exports or reruns detectors.

## Entry points

```python
load_suite(root, *, manifest=None, rscript=None, strict_inputs=False) -> dict
analyze_saved(data) -> dict
render_publication(result, output_dir, *, title="DefenseViz", **layout_options) -> list[Path]
```

`defenseviz from-suite ROOT --out OUTPUT` discovers one per-genome output directory
or immediate children of a cohort directory. It also accepts `--manifest FILE`,
`--supplement FILE`, `--rm-evidence FILE`, `--rscript PATH`, and `--strict-inputs`. The manifest has
`schema_version: 1` and `samples` entries with `sample_id`, `path`, `label` and
optional `highlight`, `genus`, `species`, `cohort`, `modules`. Manifest paths are
relative to the manifest. Basenames are literal default sample identifiers;
accessions and species are never guessed from filenames.

## Canonical analysis inputs

All arrays default to empty, except `samples`. No missing array proves a negative.
Identity is `(sample_id, gene_id)`, not a concatenated string or protein hash.
Equal sequences at different loci remain distinct genes. A `source_ref` is an
object containing `path`, `sha256`, and optionally `row`, `sheet`, or `record_id`.
`source_refs` is a list of these objects. Every imported file appears in `inputs`
and `input_sha256`. Optional cross-source reviews must bind the current protein
hash and at least one exact native result artifact for every reviewed sample.
Manifest/map bindings alone are insufficient. Protein hashes use the uppercase
amino-acid string with whitespace removed; terminal/internal `*` are preserved,
matching the imported sequence exactly. Sample IDs are unique, and duplicate
resolved sample directories are rejected.

| Key | Record contract |
| --- | --- |
| `samples` | `sample_id`, `label`; optional `highlight`, `genus`, `species`, `cohort` |
| `genes` | `sample_id`, `gene_id`; optional `protein_sha256`, `protein_length`, `contig`, `start`, `end`, `strand`. Coordinates are 1-based inclusive. |
| `calls` | Existing canonical call fields: `sample_id`, `gene_id`, `caller`, `family`, `category`, `evidence_kind`, `system_id`, `role`, `rm_type`; preserve `native_model`, `native_family`, `run_id`, `source_ref` when available. Only final members enter the native union. Reference and HMM hits remain separate. |
| `coverage` | `sample_id`, `caller`, `status`: `complete`, `partial`, `missing`, or `not_run`; optional endpoint/run metadata. Completion and current artifact consistency are separate. |
| `saved_decisions` | `sample_id`, `gene_id`, `state`: `retained`, `not_retained`, `unassessed`, or `out_of_scope`; optional `protein_sha256`, `family`, `role`, `rm_type`, `source_ref`, `score_total`, `context_state`, `orphan`. Missing scores do not invalidate an explicit saved decision. |
| `annotations` | `sample_id`, `gene_id`, `protein_sha256`, `status`, `assignment_kind`: `reviewed_RM_component` or `noncatalytic_component`; `family`, `role`, optional `rm_type`, `include_in_union`, `source_refs`, `rule_id`, `context_state`, `orphan`. |
| `effector_evidence` | `sample_id`, `gene_id`, `protein_sha256`, `status`, `assignment_kind`: `individual_nuclease_effector`, `nuclease_other_role`, or `uncertain_domain_function`; optional `family`, `motif_state`, `motif_grade`; `source_refs`, `rule_id`. |
| `groups` | `group_id`, `sample_id`, `provider`, `run_id`, `membership_kind`: `native_system` or `saved_context`; `members: [{gene_id, protein_sha256?, retained?}]`, `source_refs`, optional `context_state`, `required_roles`, `family`. No grouping inferred from proximity alone. |
| `complexes` | `complex_id`, `group_id`, `sample_id`, `status`: `complete`, `incomplete`, or `unassessed`; `core_gene_ids`, `source_refs`, `rule_id`. Complete required cores can support a complex-effector endpoint without assigning individual catalysis. |
| `source_metadata` | Source/run objects including provider name, source kind, available software/database versions and provenance. Unknown producer versions remain unknown. |
| `options` | Presentation and explicit source/category mapping options; never thresholds fitted to an evaluation cohort. |
| `suite` | Native observations, module metadata, raw feature/context/motif/structure/profile records, identity maps, conflicts, capabilities and warnings. The importer also supplies the canonical tables above when fields are available. |

The analysis validates contradictory duplicate assignments, mismatched hashes,
unknown genes, invalid subsets and stale reviews. Evidence status is explicit;
accepted/supported records contribute only to their named endpoint. Category
aliases and mechanism stripes change display membership, never catalytic status.
Exact component rules may transfer an annotated role under a documented rule;
family membership alone cannot assign every member as an effector.

`options.assessment_coverage` may declare `effector`, `associated`, and `motif`
review coverage as `complete`, `partial`, or `unassessed`. This describes the
supplied assessment scope; it is not inferred from one positive record and does
not establish a complete biological inventory. Omit it for conservative observed
lower bounds. The declaration and per-endpoint states are exported with results.

Supplement records always carry `origin: "supplement"`. They cannot add
`saved_decisions`, claim native systems/provider/run IDs or set member retention.
Supplement groups use `saved_context` and provider `Review`. Their run ID is the
supplement hash. Explicit reviewed annotations remain separate from native and
DNMB-retained endpoints. Groups are scoped by `(sample_id, group_id)`.
Positive annotation/evidence statuses are `supported`, `accepted`,
`supported_annotation`, and `supported_putative_effector_annotation`. Nonpositive
states are `unassessed`, `unsupported`, `not_supported`, `rejected`, `incomplete`,
`conflicting`, `out_of_scope`, and `not_applicable`. Other spellings are rejected.
Complex `complete` denotes the recorded complete genomic component set, not an
assembled or experimentally active complex.

## Set definitions

- Native union: distinct final members from configured native sources.
- DNMB integrated: native union plus saved retained annotations and explicitly
  included reviewed components. Reference matches alone do not enter this set.
- Individual effectors: supported individual nuclease-effector annotations.
- Effector CDSs: individual effectors plus required cores of complete reviewed
  effector complexes. Export both contributions separately.
- Putative nuclease-associated CDSs: integrated members of groups that directly
  intersect individual effectors or an explicitly complete reviewed complex,
  plus ungrouped effectors. Never expand recursively through overlapping groups.
- Prioritized CDSs: restriction-role annotations retained by DNMB or explicitly
  reviewed, union effector CDSs. This is annotation-based prioritization; activity
  is not inferred. A restriction-only candidate does not become a new effector.
- Motif/site evidence: a graded subset with its original method and completeness;
  it is neither synonymous with effectors nor required for all effectors.

Every category/family number is an intersection of distinct gene sets. Rows may
overlap. A source-on-type-pool count is not native subtype agreement. R/M/S/A
components and unique CDS totals are separate; fused RM counts once in the total.
Orphans and missing retained partners require explicit saved context assessments.
Type IV does not inherit a Type I–III methylation-partner requirement.
R means a restriction-role annotation, M methylation, S specificity and A an
auxiliary component. A context defect does not by itself establish inactivity.

## Source-bound R-M reevaluation

`apply_rm_evidence(data, path)` runs before review supplements and aggregation.
The JSON object requires `schema_version: 1`, `input_sha256`, `selected_policy`,
`records`, and `candidate_scope`. Paths resolve relative to that JSON file.
The full `selected_policy` object must equal the packaged selected preset;
this route does not fit or override thresholds.

Each record requires `sample_id`, `gene_id`, `protein_sha256`, `evidence`, and
hashed `source_refs`. `evidence` supplies every field of the current `Evidence`
dataclass, including explicit optional-axis values. Missing boolean assessments
cannot become false. Exact sample/gene identity performs the join; a matching
protein hash alone never merges paralogs or resolves an identifier.
Optional `context_state`/`orphan` need separate `context_source_refs`.
`old_context` is a scoring axis, not a current retained-system assessment.
Optional `expected_decision` contains `keep`, `score_total`, and hashed
`source_refs`; a mismatch against fresh evaluation is an error. New evidence
without a prior evaluation does not need an expected decision.

Each candidate scope names a bound `sample_id`, distinct `gene_ids`, hashed
`source_refs`, and `status` (`complete` or `partial`). Complete scope cannot hide
unmapped native candidates, unresolved candidate-count conflicts, omitted
decisions, or added candidates outside the original native candidate set.
Partial scope exports the remaining unassessed and added identities explicitly.

`decision_mode` defaults to `evaluate_unassessed`, which rejects a conflict with
an existing assessed decision. Explicit `replace_scoped` selects the fresh
evaluation for the declared samples: all their previous decisions are preserved
in the audit record, and omitted candidates become unassessed. Their original
tool/reference calls are unchanged. Native saved-retention coverage and
old-to-current transitions are also preserved. This selection never edits a
source artifact.

The policy's `integrated_R_component_support` remains separate from retained R
annotation counts. Primary `priority` preserves annotation-based triage: retained
or reviewed R roles union effectors. The separate `priority_support_qualified`
endpoint requires positive policy R support, an independently accepted R review,
or effector evidence. An unmet policy nuclease route is not an R-negative or an
inactivity result. Both membership flags and per-CDS three-state policy support
are exported; `R_only_policy_unmet_CDS` identifies candidates retained in the
broad endpoint without that narrower support. No computational support field
establishes experimental activity. Pre-retention
native contexts do not become current assessments merely because all members
were retained after reevaluation.

Review groups can explicitly name `supersedes_group_ids` for earlier DNMB saved
contexts. They must belong to the same sample, stay within the original member
boundary and carry hashed sources; native detector systems cannot be replaced.
Displaced contexts remain in `suite.superseded_context_groups`. Complex records
cannot reference a displaced active group.

See the [runnable synthetic example](../examples/suite_rm_evidence/README.md).
The figures retain `+` for partial counts. `rm_evaluation.json` records every
Evidence object, policy decision, source binding and actual evaluation count.

Type tables condition on their observed assigned-type pool. Their
`within_pool_counts_status` and `within_pool_role_states` determine the displayed
assessment state. For example, all nine assessed members of a nine-CDS pool are
exact within that pool even if an unrelated historical nominee is unassessed.
An unknown member inside the pool keeps its DNMB subset partial. Broader
`counts_status` and sample retention coverage remain separate exported fields;
a zero inside an empty observed typed pool is not genome-wide absence evidence.

Atlas counts use the same separation of endpoints. `current_saved_decisions_by_sample`
reports completeness within an explicit current candidate scope; absent or
unassessed members inside that scope prevent completeness. Without a declared
current scope, the broader saved-decision state applies. Current retained counts
and component roles also include explicitly retained active decisions outside
that selection; both assessment member keys and outside-selection retained keys
are exported separately. Original-source intersections and
observed input pools use their own source coverage. Fraction numerators and
denominators carry independent states. Historical candidates and superseded
contexts remain auditable; they do not make unrelated current counts partial.

`rm_candidate_scope.member_keys` preserves the complete selected evaluation
scope. `typed_member_keys` and `untyped_member_keys` partition it without dropping
decisions. Only typed current candidates supplement original R-M/reference
memberships in the atlas input denominator. The typed table pool, the broader
atlas input pool, and the full evaluation candidate scope are distinct endpoints.

## Analysis outputs for the renderer

`analyze_saved` returns `schema_version`, `samples`, `sources`, `source_sets`,
`global_summary`, `sample_summary`, `overlap_regions`, `priority_overlap_regions`, `inventory_rows`
(with `categories` and `families`), `cells`, `rm_components`, `type_tables`,
`gene_ledger`, `priority_ledger`, `groups_selected`, `endpoint_definitions`,
`warnings`, `availability`, and `provenance`.

Inventory rows expose `id`, `label`, `category`, `member_keys`, `N`,
`source_counts`, `associated_CDS`, `effector_CDS`, `motif_CDS`, `counts_status`,
and number of genomes; cells add `sample_id` and `family_id`. Source masks carry
an explicit `source_order`, counts and exact member keys; preserve a selected
all-zero mask outside all source sets. Comparisons must use a common explicit
sample/endpoint scope. Missing or incomplete evidence is shown as unavailable or
an observed lower bound, not a biological zero.

Publication views use editable SVG/PDF. Main atlas and type panels are separate
from Letter portrait comparisons. Numeric columns are right-aligned; associated
Total and Effector CDSs use independent columns; DF/PAD/DP is the default table
order. Category colors, purple nested atlas glyphs, pink sample highlights and
uniform comparison tracks are data-independent presentation defaults. Empty
families can be hidden only with an explicit observed-zero inventory; unknown
families and missing evidence remain visible. Large cohorts, long labels and
additional sources require pagination or a counted-intersection fallback.

## Manuscript display projection

`from-suite` adds `manuscript_atlas` after `analyze_saved`. This independent
projection never changes gene retention, effector admission, priority selection
or source membership. Its original R-M intersection uses native sources only;
reference matches are a separate column. I–III component glyphs use saved retained
components; IV also includes explicitly reviewed components already admitted by
the analysis. Exact DNMB saved contexts produce I MS/RMS and III M/MR group
counts. Equal original member sets are deduplicated only within the same known
sample/provider/run/type; adjacent or overlapping groups are not merged.

Optional manifest `options.manuscript_atlas` accepts `source_order`, per-sample
`group_coverage` (`complete`, `partial`, `unassessed`), `absence_families`,
`absence_label` and `taxonomy_legend`. Group coverage must describe the saved
exports; without an explicit complete assessment, observed groups remain partial.
An absence highlight requires both a configured family and a complete zero cell.
It never means genome-wide biological absence. Sample `display` metadata controls
taxonomy colors and taxon/strain typography. There are no default target-strain
or transformability keys. Layout `highlight_label` is nullable text and labels
explicit pink highlights only when present. Optional layout `sample_marker`
requires exactly `field`, `label`, and `color`: a marker is drawn only where
`sample.display[field] is true`, with its key on the same page. Legacy
`display.transformable` is ordinary metadata and does not activate a marker.
Taxonomy/category keys require rendered members; an absence key requires at
least one actually shaded complete-zero cell on that page.
Family display metadata may carry `display_order`, `display_group`,
`partner_group`, `category_tags`, `secondary_category`, `secondary_style` and
`highlight_source`. These affect presentation only; no function follows from a
label or color. Highlight ticks use actual exclusive source memberships.

Historical R-support rings and curated reference-group badges are separate from
nuclease-effector evidence. Optional annotation records use
`evidence_kind: "rm_role_display"` or `"rm_context_display"`,
`assignment_kind: "reviewed_RM_component"`, `status: "unassessed"` and
`include_in_union: false`. Each binds `sample_id`, `gene_id`, current
`protein_sha256`, a `rule_id` and hashed `source_refs`. The explicit fields
`rm_R_support` and (for contexts) `rm_reference_link` accept `supported`,
`not_supported` or `unassessed`. Context records also bind `group_id` and an exact
retained group member. A curated reference badge requires a complete saved R-M
context; raw reference presence does not supply it. Conflicting bindings or
assessments are rejected. The supplement importer checks actual source hashes.

The five workflow cards show the saved method. Optional layout
`manuscript_workflow.parameters` displays supplied saved criteria; it does not
evaluate them. `analysis.manuscript_atlas.provenance`, exact group/member records
and `publication_layout.json` retain the origins and geometry for auditing.

Evidence-only HMM providers are excluded from manuscript native-source columns.
The broader legacy summary coverage can remain conservative when an additional
HMM-only provider is present; this does not create a final call or a zero-hit claim.

Distinct family pools with the same display label retain separate counts and
identities. Rendering falls back to their exact declared IDs to disambiguate
those labels consistently in atlas and comparison views; no pools are merged.
