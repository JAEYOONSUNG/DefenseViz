# Use saved DNMBsuite results

DefenseViz reads saved DNMBsuite outputs. It does not install or rerun DNMB,
DefenseFinder, PADLOC, DefensePredictor, REBASE searches or structure prediction.
The public analysis pipeline name is **DNMB**; historical source paths and module
names remain in provenance.

## Install and run

```sh
python -m pip install '.[suite]'
defenseviz inspect-suite /path/to/DNMBsuite_outputs
defenseviz from-suite /path/to/DNMBsuite_outputs --out defenseviz-results
```

The input can be one genome's result directory or a parent containing per-genome
directories. Discovery is limited to immediate children. Original files are
read-only, and the output directory must be new. `inspect-suite` validates and
reports available endpoints without generating figures.

R is used only to read saved RDS completion/version metadata through a fixed
packaged helper. Set `--rscript /path/to/Rscript` when it is not on `PATH`.
No original DNMB R code is sourced. Missing completion evidence must remain
explicit; a directory without results is not a successful zero-hit run.

## Select and label samples

Default sample IDs and labels are literal directory basenames. To choose a cohort,
order rows, supply concise display labels or highlight strains, use a manifest:

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
IDs are opaque: the importer does not parse accessions, guess species, or merge
different loci because their proteins have the same sequence.

For nonstandard locations, sample-level `genes` and `fasta` paths can be supplied.
`modules` accepts `defensefinder`, `padloc`, `defensepredictor`, and
`rebasefinder` mappings with explicit artifact paths (`calls`, `id_map`,
`profiles`, `workbook`, `systems`, `fasta`) and a saved `status` when needed.
All such paths are relative to the manifest. Preserve the actual run status;
declaring completion is not a way to fill missing output with zeros.

## What is preserved

Native DefenseFinder and PADLOC final memberships, general DefensePredictor
predictions, raw REBASE matches and saved DNMB decisions remain distinct.
Repeated system memberships do not inflate gene totals. HMM hits are evidence,
not final system members. Original native source counts never acquire DNMB-only
additions. Native system identifiers remain scoped to their sample and source.

The importer supports version-dependent augmented REBASE tables and legacy
workbooks. Explicit saved `curation_keep` and final role/type fields are preserved
when available. Older exports without those fields remain unassessed; their
historical BLAST filter is not substituted for modern DNMB retention. Raw features,
motif observations, structures and annotation transfers are not interchangeable.

Automatic effector annotations use a versioned catalog of exact component models;
required complete complex cores are counted separately from individual effectors.
Raw InterProScan and EggNOG annotations are preserved for source-bound review.
Generic domains and a broad upstream `rease_motif_verified` flag do not establish
nuclease function: that flag can also describe motor-associated evidence.
Consequently, older outputs without explicit functional reviews do not reproduce
the expanded manuscript review automatically. Supply those saved reviews when
available; unsupported endpoints remain unassessed or observed lower bounds.

Current exported tables take precedence over cached tables. If live artifacts and
saved run metadata disagree, the conflict is reported rather than merged away.
Use `--strict-inputs` to reject these conflicts. A saved completion sentinel proves
run completion, not that an export was never subsequently edited.

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

The output includes an editable publication PDF, SVG pages, layout measurements,
`analysis.json`, `native_bundle.json`, tabular gene/priority/family/category
results and `receipt.json`. The receipt pins inputs, implementation files,
available tool/database metadata, command, elapsed times and output hashes.
Unknown producer versions and absent database hashes are not invented.

Input references are exported relative to the supplied root when possible.
Execution receipts deliberately retain exact executed command arguments and may
contain machine-specific paths, as may unchanged upstream metadata. This is
execution provenance, not a claim that every provenance field is relocatable.

The main atlas uses category colors and purple evidence glyphs, with optional pink
sample highlights. Comparative tables separate putative nuclease-associated
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
  "source_order": ["DefenseFinder", "PADLOC", "DefensePredictor"],
  "type_samples": ["isolate-a", "isolate-b"],
  "preview_png": true
}
```

Omit `type_samples` for the default selection; use `[]` to omit type tables.
`max_rows_per_page` and `max_families_per_page` limit atlas pages.
`atlas_page_size` sets maximum width/height in millimetres (at least 270 × 145);
smaller datasets use less space. `category_order` controls displayed category
order; an empty list preserves category encounter order while grouping partners.
`atlas_section_titles` overrides section headings. Comparisons remain US Letter portrait and continue onto additional
pages. Additional sources or large intersection numbers use an exact intersection
table when a Venn diagram cannot remain legible. `source_colors`,
`category_colors`, and `aliases` are optional mappings. Unknown layout options
raise an error. The manifest records the actual font files and hashes; an
unavailable Arial installation is never reported as Arial.

These are saved annotation/evidence summaries. They do not establish biological
activity, prove a candidate is safe to remove, or supply independently adjudicated
precision/recall. Unreviewed results are not false positives or true negatives.
Family/category rows overlap, so their counts must not be summed as a genome-wide
gene union.
