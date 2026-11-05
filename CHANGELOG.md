# Changelog

## 2.4.0

- Make source-bound R-M feature collection, frozen-policy evaluation and current-context review the default DNMBsuite path.
- Evaluate exact protein motif architecture and preserve unknown measurements through conservative score bounds; export every candidate disposition.
- Review all observed systems and calculate functional count coverage over each actual displayed CDS set.
- Support explicit current-reference homology evaluation when legacy exports lack coverage; preserve original reference and detector memberships.
- Keep stored-decision replay and the explicit reviewed-evidence path available.
- Preserve incomplete and quarantined R-M scopes; never infer orphanhood from a retention cutoff.
- Merge reference HSPs under native measurement rules and separate subject-specific negatives from missing global evidence.
- Audit field-level assignment/coverage changes and distinguish attempted evaluations from resolved decisions.

## 2.3.1

- Read legacy source taxonomy, gene worksheet labels and supplemental reference identities without changing native memberships.
- Qualify incomplete priorities and improve short-cohort label spacing.

## 2.3.0

- Add source-pinned offline system information from the official DefenseFinder wiki, DefenseFinder definitions and CasFinder definitions.
- Export a searchable catalogue and automatic per-observation metadata beside DNMBsuite and normalized-project analyses.
- Keep wiki knowledge, model identity, version compatibility and query evidence distinct; no changes to accepted scientific memberships or figure categories.

## 2.2.0

- Review every observed native system and ungrouped detection, including unknown families, with distinct evidence axes, reason codes and a review queue.
- Add exact configurable review rules and source-verified domain/reference evidence; preserve accepted figure memberships while exposing new review proposals.
- Scope functional assessment declarations and conservative defaults by sample, preventing unrelated samples from changing each other's count status.
- Preserve InterPro evidence for all observed/reviewed identities with sequence checks, explicit exclusions and quarantine accounting.
- Export review JSON/TSV and input/rule provenance alongside the existing publication figures.

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
