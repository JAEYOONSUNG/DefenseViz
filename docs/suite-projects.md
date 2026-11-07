# Repeatable DNMBsuite projects

Put a `defenseviz_project.json` file beside the inputs and settings for a comparison:

```json
{
  "schema_version": 1,
  "root": "DNMBsuite_outputs",
  "engine_config": "engine.json",
  "evaluation": "full",
  "layout": "layout.json"
}
```

Run from any working directory:

```sh
defenseviz from-suite /path/to/project --out results
```

The project may also be supplied as an explicit JSON path. Its paths resolve against that file. An explicit command-line option overrides the corresponding project setting. A folder without `defenseviz_project.json` uses ordinary recursive suite discovery.

Supported fields are `root`, `manifest`, `supplement`, `rm_evidence`, `layout`, `engine_config`, `reference_cache`, `system_review_rules`, `system_review_evidence`, `evaluation`, `strict_inputs`, `rscript`, `title`, `external_type_projects` and an optional `description`. Unknown fields are rejected. The project file and all consumed inputs are hashed in the receipt.

## Saved evidence

For DNMB outputs with sufficient saved evidence, use `evaluation: "auto"` without an engine configuration. A source-bound `rm_evidence` file requests policy evaluation of supplied measurements; a `supplement` supplies reviewed annotations. Those files contain CDS identities and evidence, not expected aggregate counts. Their use is recorded in the receipt.

`stored` replays saved decisions and skips new family-scope inference. It cannot be combined with `rm_evidence`, `reference_cache` or `engine_config`. To switch an existing project to replay, remove incompatible settings or create a separate replay project. External projects have their own evaluation settings; the main command's overrides do not silently change them.

Saved DNMB curation is imported separately from downstream decisions. Source rows, protein identities, native fields and unavailable fields remain explicit. Legacy InterProScan 13-column TSV files are accepted when only the optional GO and pathway columns are absent. Malformed rows or mismatched sequence bindings are not repaired by guessing.

## Letter pages

Example `layout.json`:

```json
{
  "preview_png": true,
  "page_composition": "letter_composite",
  "type_samples": ["sample-a", "sample-b", "external-case"],
  "rm_fraction_scope": "retained_I_III"
}
```

`letter_composite` puts the landscape atlas and selected R–M type tables together on a Letter landscape page when they fit. Source overlap and family comparisons use Letter portrait. Larger cohorts, more families or additional case tables continue onto further pages; no rows are dropped to force two pages. The layout ledger records native geometry, final physical geometry and the export scale.

The default `page_composition: "separate"` keeps the atlas, type tables and comparison views separate. The default fraction scope is `assigned_I_IV`. `retained_I_III` explicitly displays the policy-retained types I–III subset, calculated from CDS identities. This display option does not remove Type IV assignments from exported tables or other panels.

## Independent case tables

An external example can be displayed beside primary cases without entering cohort totals:

```json
{
  "schema_version": 1,
  "root": "primary_outputs",
  "evaluation": "auto",
  "layout": "layout.json",
  "external_type_projects": [
    {"project": "external/defenseviz_project.json", "sample_ids": ["external-case"]}
  ]
}
```

Each external project goes through the same input validation and analysis. Its selected sample IDs must be distinct from the primary cohort and other external examples. Nested external projects are rejected. External cases affect only the type tables: primary source sets, Venn counts, family counts and workbook totals stay within the primary cohort.

`external_type_cases.json` binds the displayed cases to separate analyzed JSON artifacts. The receipt records each external project's inputs, execution mode and search counts. Top-level search and scoring run totals include external execution; policy decisions and biological totals remain separately scoped. `inspect-suite` also validates external projects rather than ignoring them.
