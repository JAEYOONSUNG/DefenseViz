# Saved input formats

DefenseViz reads local, saved result files. Imports never launch a detector, query a database, or infer sample IDs from filenames. Paths resolve relative to the project JSON, even when the command runs elsewhere. Generated source references are relative to that project when possible (Windows inputs on another drive retain an explicit absolute path); resolved paths and SHA-256 hashes of the bytes actually read remain available for the execution receipt. UTF-8 labels and identifiers are preserved. TSV and CSV are preferred; JSON tables use a list of objects. XLSX requires `pip install 'defenseviz[excel]'` and an explicit `sheet` when a workbook contains multiple sheets. Formula cells must be exported to literal values first.

## Project configuration

```json
{
  "schema_version": 1,
  "samples": "samples.tsv",
  "calls": [
    {"format": "defensefinder", "path": "isolate_a_genes.tsv", "sample_id": "isolate_a", "id_map": "a_id_map.tsv"},
    {"format": "padloc", "path": "isolate_a_padloc.csv", "sample_id": "isolate_a", "id_map": "a_id_map.tsv"}
  ],
  "family_map": "families.tsv",
  "evidence": "evidence.json",
  "options": {"title": "My isolates"}
}
```

`samples` requires unique `sample_id`; optional columns are `label`, `genus`, `species`, and `highlight` (true/false or 1/0). Sample order is display order. Identifiers must be nonempty text without surrounding whitespace. Metadata never changes score decisions.

For any table input, `{ "path": "file.xlsx", "sheet": "Sheet name" }` can replace its filename. Call sources are always objects including `format` and `path`.

## Native caller adapters

| Format | Required native columns | Interpretation |
|---|---|---|
| `defensefinder` gene table | `hit_id`, `sys_id`, `type`; optional `subtype` | One selected system member per row. Display family is subtype when present, otherwise type. |
| `defensefinder` systems table | `sys_id`, `type`, `protein_in_syst`; optional `subtype`, `genes_count` | Expand the comma-separated, explicit protein member list. If given, genes_count must agree. |
| `padloc` per-gene CSV | `system.number`, `seqid`, `system`, `target.name` | Preserve selected genes; system numbers are scoped by contig as well as sample and caller. |
| `defensepredictor` output | `locus_tag`, `mean_log_odds` by default | Require a caller-source `score_threshold` supplied by the user; include scores >= threshold. This is a general-defense context endpoint, without native R-M classes or system IDs. |
| `normalized` | `gene_id`, `family`, `sample_id`, `caller` | `sample_id` and `caller` may instead be fixed in the source descriptor. |

DefenseFinder accepts leading `# ` metadata comments from its native systems exports. Import either a genes or systems file for a caller run, not both: duplicate call identities are rejected. A genes table and its parent system use the same `sys_id`. The adapter does not reinterpret profile names as catalytic roles.

PADLOC `system_id` is a JSON tuple string `[seqid, system.number]`, so equal local numbers on different contigs remain different systems. DefenseFinder native type labels `RM_Type_I` through `RM_Type_IV`, and PADLOC labels `RM_type_I` through `RM_type_IV`, carry their explicit R-M type/category. Other family names are preserved and classified only by the shipped exact family catalog or an explicit user mapping. Unmatched families stay `Unresolved`. The catalog supplies display/category metadata, never catalytic roles or score support.

For DefensePredictor, `score_column` and `id_column` may override the default columns. Threshold choice is explicit and recorded, not fitted by DefenseViz. Negative rows are still validated; they are omitted from included calls but the supplied per-sample file remains a completed result. Missing/invalid scores are errors. The adapter never extracts IDs by splitting `protein_context_id`. Supply that column as `id_column` with an explicit ID map if it is your source identifier.

A supported native export establishes selected-member membership only. Raw HMM evidence must use the normalized schema with `evidence_kind=hmm_hit`; REBASE/reference hits use `reference_hit`. Neither is silently promoted to a selected defense system. No raw HMM parser guesses model-specific qualification thresholds.

### Identifier mapping

All cross-source joins use the explicit `(sample_id, gene_id)` pair. Without an ID map the source ID is used verbatim: the user must ensure callers share the same identifier namespace. An optional `id_map` table has `source_id`, `gene_id`, and `sample_id` (the sample column can be omitted for a per-sample source). The exact exported aliases `query_id`/`locus_tag` are accepted instead. Maps must be one-to-one within each sample; unmapped IDs and collisions fail. Source identifiers containing delimiters are never parsed.

### Normalized calls

Optional columns are `category`, `system_id`, `role`, `rm_type`, `evidence_kind`. Defaults are `category=Unresolved`, `evidence_kind=final_member`, and null for unspecified components. Roles are `R`, `M`, `S`, `RM`, `C`, or `Aux`. R-M types are `Type I` through `Type IV`. A supplied R-M type requires category `R-M`.

`family_map` is an exact mapping table with `family`, `category`, optional `caller`, `display_family`, `rm_type`. Explicit category values in normalized rows (including `Unresolved`) override the shipped catalog. A user `family_map` conflicting with an explicit non-Unresolved row category is rejected. Explicit user mappings override shipped defaults; within the same source, caller-specific entries take precedence over generic entries. `display_family` permits explicit harmonization of native labels; original names remain in `native_family`. Categories are `R-M`, `Modification`, `Nuclease`, `CRISPR-Cas`, `Abortive infection`, `Other mechanisms`, or `Unresolved`. Conflicting category/type claims are rejected. Category membership alone does not assign catalytic roles or add score points.

## Complete, missing and not run

A supplied per-sample final-member source file is marked `complete`, including a header-only zero-call file. A normalized source containing only HMM/reference evidence does not establish native final-member completeness; those records remain available even when native coverage is `missing` or `not_run`. Providers contributing only `reference_hit` records are evidence sources, not native caller columns; only explicit native coverage can add such a provider to the native caller list. HMM-only caller output remains a native caller with missing final coverage. A zero-byte file has no verifiable schema and is rejected. An empty normalized file also needs `caller` in its descriptor.

For multi-sample call tables, provide `coverage.tsv` with `sample_id`, `caller`, `status`. Allowed statuses are `complete`, `missing`, `not_run`. Final-member calls require complete coverage; omitted sample/caller pairs become `missing`, never zero. Complete coverage means a completed *saved annotation output*, not biological absence. An explicit status conflicting with a supplied per-sample completed source is rejected. Coverage can list callers with no observed calls.

## DNMB reevaluation inputs

Canonical `evidence.json` is a list of `{sample_id, gene_id, evidence, provenance}` records. `evidence` contains the complete typed DNMB evidence object described in [evidence.md](evidence.md); old stored `decision` values are ignored. Null evidence may explicitly preserve an unassessed candidate and its provenance. Partial malformed evidence objects do not become valid scores. Duplicate sample/gene records fail.

Existing DNMB `evidence_and_decisions.json` exports with opaque `target_id` are accepted through:

```json
{"evidence": {"format": "dnmb", "path": "evidence_and_decisions.json", "id_map": "evidence_ids.tsv"}}
```

The map requires `target_id`, `sample_id`, `gene_id`, all one-to-one. An optional `sequence_sha256` column or an explicit top-level hash in each DNMB record binds exact sequences for reference/domain enrichment; if both are supplied they must agree. Hashes are validated and copied into the new provenance. Nested historical adjudication fields are not guessed as sequence bindings. DefenseViz never splits a target ID to guess its sample. Upstream adjudication is retained as provenance; the stored decision is not reused as a fresh calculation.

### Native feature tables and older DNMB exports

```json
{"evidence": {"format": "dnmb_features", "path": "native_features.tsv", "sample_id": "isolate_a", "id_column": "locus_tag"}}
```

This adapter invokes the original `from_feature` allowlist after verifying that required fields were actually supplied. It accepts TSV, CSV, JSON-row tables and XLSX. Identifiers may come from an explicit `id_column` or `gene_id`, `locus_tag`, `query_id`, `target_id` in that order, with an optional one-to-one `id_map`. Exact native aliases `rm_type`/`subunit` can supply `family_id`/`enzyme_role`. A supplied `sequence_sha256` column is validated and retained in provenance, allowing the optional reference/domain extension to verify an exact sequence join; no hash is invented when it is absent.

Complete reevaluation needs native factor flags (`factor_strong`, `factor_moderate`, `factor_profile`, `factor_annotation`, `factor_motif`, `factor_context`, `factor_alternative_annotation`, `factor_conflict`, `factor_partial`), sequence length, reference and curated-reference identity/coverage/E-value/compatibility fields, role/motif qualification flags, explicit structural admissibility and divergent-context flags, plus a `profiles` JSON list containing individual coherent Profile records. Empty `profiles=[]` is valid; separate per-column HMM maxima cannot establish a coherent profile and are not substituted. `blast_family_compatible` can be unknown, and numeric observations can be missing, but an unassessed boolean is not converted into an evaluated false predicate.

Older raw REBASE tables often lack these factors. Their candidate identities and supplied fields are preserved with `evidence=null`, `status=unassessed`, and a `missing_fields` list. They are not assigned a zero score, and absent coordinates or structural measurements are never fabricated. This is an honest import path for old DNMB outputs, not a claim that old workbooks contain all current evidence. Supply complete native feature exports or typed evidence from the upstream analysis for full reevaluation.

## Optional coordinates

`genes` requires `sample_id`, `gene_id`, `contig`, `start`, `end`, `strand`; optional `sequence_sha256` must be a lowercase SHA-256. Coordinates are 1-based inclusive, with `start <= end` and strand `+` or `-`. Duplicate sample/gene records fail. Coordinate availability alone does not grant context or activity support.

## Migration from DefenseFinderViz

The earlier DNMB `DefenseFinder_Heatmap` in `R/DefenseFinder_plot.R` read `*_finder_systems.tsv` from the working directory, derived sample names from filenames, and consulted `Genome_summary_df` in R's global environment. DefenseViz replaces those implicit dependencies with a versioned project, explicit samples, paths and joins. It accepts the same native systems-member lists, while distinguishing their system counts from distinct CDS counts and preserving multi-caller provenance. No R global state or manuscript-specific genome list is needed.

## Optional reference/domain extension inputs

`evidence_extension` supplies the generic reference/domain assessment stage:

```json
{"evidence_extension": {
  "registry": "registry.json",
  "candidates": "candidate_evidence.json",
  "admission_rules": "admission_rules.json",
  "source_root": "."
}}
```

All three JSON files are read and hashed by the loader. `registry_dir` resolves to the registry file's parent, while `source_root` (default project directory) locates separately bound evidence files. The analysis stage validates reference sequences, exact model identities, source hashes and domain/site correspondence before admitting evidence; loading JSON does not itself grant support. See [evidence.md](evidence.md) for the registry and admission contract. No database downloads or new searches happen during import.
