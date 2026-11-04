# Changelog

## 2.1.4

- Correct atlas uncertainty states: original-source intersections, observed R-M pools and current retained counts use their own evidence scopes. Historical unassessed nominees remain in the audit without marking unrelated current counts as partial.
- Preserve genuine incomplete current evidence, unknown component roles and unassessed Type IV contexts.
- Restore circular section numbers and group-anchored titles on short inventories; measured title collisions determine vertical spacing.

## 2.1.3

- Connect complete, source-bound R-M evidence to the DNMBsuite workflow with `--rm-evidence`; export fresh policy decisions and their provenance alongside the figures.
- Preserve original tool/reference calls and unassessed legacy candidates; reject incomplete evidence and stale protein/source bindings.
- Reconcile retained context completeness after filtering and count reciprocal saved contexts once.
- Preserve multipart DefenseFinder systems, saved DefensePredictor decisions and incomplete REBASE observations without conflating source coverage with candidate coverage.
- Normalize known native family names for display while retaining source identities and explicit user mappings.
- Export annotation-based prioritization separately from its stricter policy-support subset.
- Keep the full manuscript atlas, type tables and Letter portrait comparisons for small cohorts.

## 2.0.0

- Standalone Python API and CLI with explicit file/config inputs instead of working-directory and R global-environment dependencies.
- Native caller adapters and identity mapping, typed DNMB evidence reevaluation, separated gene/system and final-member/HMM/reference endpoints.
- General comparative Figure 2 a/b rendering with editable SVG, pink sample highlights and consistent bar geometry.
- Transactional output, table exports, input/software hashes and reproducibility receipts.
- Synthetic examples, adversarial tests and clean-install validation.

This is a new interface, not a drop-in replacement for DefenseViz 1.x R functions. Native detector and recognition-specificity execution is not included.
