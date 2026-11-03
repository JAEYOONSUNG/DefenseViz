# Evidence and counting contract

DefenseViz evaluates saved annotations and explicit DNMB evidence. It does not
run a defense detector, infer biological activity, or estimate independently
validated sensitivity, specificity, precision, recall or AUC. Empty or missing
observations are never evidence of biological absence.

## Units and identity

Every gene identity is the pair `(sample_id, gene_id)`. IDs are opaque strings;
punctuation, case and delimiters are preserved. Identical sequences in different
genes or samples are not merged. Duplicate evidence records and duplicate call
identities raise errors. Multiple distinct family calls for one gene are kept.

The family matrix (`cell_scope: native_final_members`) counts distinct genes in
native `final_member` records. `hmm_hit` and `reference_hit` records remain in
`evidence_calls` and never inflate native caller counts. Explicit categories are
honored; an unspecified category becomes `Unresolved`. Conflicting categories
for the same family require an explicit consistent mapping before integration.

Reference-only providers are reported in `evidence_sources`, not automatically
as native callers; this uses evidence kind rather than a provider-name exception.
An explicit native coverage row can opt such a provider into caller scope.
HMM-only detector sources remain native caller names with missing final output
and null counts. `evidence_sources` lists providers with auxiliary records but
no supplied native final-member records; a provider can therefore be in both
lists when its native caller status is known but its final output is missing.

The sample summary distinguishes:

| Field | Meaning |
| --- | --- |
| `native_union_cds` | Union of distinct native final-member genes across callers |
| `retained_dnmb_only_cds` | Retained DNMB evidence genes outside that native union |
| `union_cds` | Native union plus retained DNMB-only genes |
| `rm_pool_cds` | R-M native/reference candidates plus supplied R-M policy candidates |
| `scored_cds` | R-M pool candidates with complete, in-scope policy evidence and a valid sequence |
| `retained_cds` | Scored candidates satisfying the selected evidence policy |
| `unassessed_cds` | R-M pool candidates without complete usable policy evidence |
| `out_of_scope_cds` | R-M pool candidates explicitly outside scored roles/types, including control/auxiliary roles |

`rm_pool_cds = scored_cds + unassessed_cds + out_of_scope_cds`. Families and
categories can overlap, so neither their sum nor a sum of role rows is the
integrated gene union. Exact caller overlap groups use native final-member
genes only. A group's `all_callers_complete` flag records whether absent caller
membership can be interpreted as absence from all supplied final outputs.

Native system identity is `(sample_id, caller, system_id)`. A system cannot be
assigned conflicting families within that identity. `system_count` is the sum
of explicitly caller-specific native systems, not de-duplicated biological
systems. If any member lacks a system ID, its aggregate system count is null
with `unassessed_missing_system_ids`; `known_count`/`known_native_system_count`
preserve the observed lower bound. A complete native export with no genes has
zero observed systems. Pooled components never establish system completeness.

## Coverage and missing values

Coverage describes **native final-member export completeness**, not HMM search
coverage. Only `complete` allows a caller's displayed zero. `missing` and
`not_run` yield null caller counts. Native final members conflict with incomplete
coverage and are rejected; HMM/reference records may coexist with incomplete
native coverage. A raw-HMM-only file therefore cannot imply zero native calls.

For direct Python calls with `coverage=None`, observed final-member sample/caller
groups are provisionally treated as supplied complete outputs, and a warning
records this assumption. Every omitted or HMM/reference-only group remains
missing. Prefer explicit coverage for partial exports. Explicit empty native
files can establish complete coverage with zero saved calls. The observed
integer count in a cell with no complete callers is zero, but `count_status`
is `unassessed` and the renderer must display NA.

## DNMB evidence policy

The numerical predicate is derived from the historical DNMB evidence policy.
`analysis.policy.ORIGINAL_SOURCE` and `ORIGINAL_SOURCE_SHA256` pin its source.
The packaged selected preset uses `r_identity_min = 0.60`, rather than the old
class definition's 0.70 default. Evidence gates were selected using historical
manuscript-cohort calibration. Weights and the retention cutoff were fixed.
This provenance does not establish independent validation, a unique optimal
threshold, or predictive performance on a new dataset. No sample identifiers,
locus lists, desired cohort counts or experimental labels occur in the scoring
predicate or admission rules.

Any parameter override is identified separately as `custom_parameter_override`,
with an exact `parameter_overrides` difference map and `is_selected_preset: false`.
The custom configuration is not described as having undergone the historical
selection or independent validation.

The main sequence evidence axis is capped at the maximum of strong reference
(35), moderate reference (25), profile (25), and annotation (10). Motif, context
and strict structure can contribute 20, 20 and 25 points respectively. Alternative
annotation and conflict each carry a 15-point penalty, with conflict suppressed
when alternative annotation is already penalized; partial evidence subtracts
10. Scores are floored at zero. Retention requires at least 50 points and the
configured capped-axis count. These axes are evidence categories, not proven
statistical independence.

`parse_evidence` requires all non-default `Evidence` fields. Booleans must be
actual JSON booleans; `false`, null, zero and a missing field are not silently
interchangeable. Numeric measurements may be null, but finite fraction values
must be in [0,1] and E-values nonnegative. A fully missing record is represented
as `evidence: null` and remains unassessed. A malformed purportedly complete
record raises `ValueError`. Imported native feature exports must explicitly
document which fields are missing before constructing evidence; the historical
`from_feature` helper alone does not establish completeness.

The scored roles are R, M, S and RM. C and Aux remain in the R-M candidate pool
but are out of scope. The policy's known R-M families are Type I–IV; the historical
NA/Unknown/Unresolved family state preserves existing evidence while preventing
new type-specific motif criteria from being invented. New non-R-M families can
appear in native matrices but are outside this R-M scoring preset. An invalid or
missing sequence makes a candidate unassessed even if legacy support flags were
present. The low-level historical `evaluate(Evidence, Policy)` predicate is
preserved for exact numerical replay; use `integrate` for the strict public
validation and scope boundary.

All HMM numeric gates must hold for **one profile**. Coverage and E-value cannot
be stitched between separate hits. New R-component admission requires typed,
compatible reference and profile evidence and the registered motif rules.
New motif score credit additionally requires context by default. A raw PLD-HKD
class string cannot replace the two-occurrence/domain rule. Type IV does not
receive a newly invented generic motor/nuclease motif score. Reference presence,
conserved residues, model confidence and active-site geometry remain distinct.

`rm_components` describes explicit candidate type/role assignments. The same
gene may have conflicting or multiple native assignments, retained as separate
rows; support is only credited to the type/role matching its evaluated Evidence.
For R/RM, `supported_count` uses `integrated_R_component_support`; for M/S it uses
retention. This count is neither complete systems nor biological activity.

## Context and structure

Coordinates are validated as 1-based inclusive intervals with an explicit
contig and strand. A sequence hash is optional for geometry-only input, but when
present it must be a valid SHA-256. If both geometry and policy provenance specify
a hash for the same sample/gene, conflicting hashes are rejected. Missing hashes
cannot supply new sequence-bound evidence. DefenseViz does not infer context from arbitrary
adjacency, does not pool components across samples/contigs, and does not infer
system completeness from pooled R/M/S counts. The current policy consumes
previously adjudicated `old_context` / `divergent_context` facts; it does not
recompute these from imported coordinates. Their producing methods and bounded
geometry belong in record provenance. Changing the configured maximum gap does
not retrospectively recompute imported context facts. Structure is likewise an
explicit precomputed evidence axis, never manufactured from domain or motif
presence.

## Optional source-bound registry assessment

`analysis.registry.assess_registry(registry, candidates, registry_dir=..., source_root=...)`
evaluates user-supplied model/reference/site/module evidence offline. It does not
run HMMER, BLAST, alignments or structure prediction. No raw biological reference
database is bundled. The preserved generic registry core is source-pinned by
its module constants.

Registries use schema `1.0` with family entries containing explicit substrate,
component, exact model accession/source/hash, and reference contracts. Reference
sequences are inline or local FASTA paths relative to `registry_dir`; each is
bound to its SHA-256 and full-sequence 1-based numbering. Catalytic sites specify
residue, role, primary evidence source, and whether they belong to the required
tested subset. Candidate records require explicit `candidate_id`, `sample_id`,
`query_id`, exact uppercase sequence and SHA-256. Neither class/domain presence
nor a motor/modification domain is silently transferred to nuclease evidence.

The default API checks actual local files against every cited candidate source
hash. `source_root` resolves relative evidence paths; network URLs are rejected.
Exact model identities/hashes must agree with the supplied registry. This checks
provenance integrity and evidence-contract consistency; it does not independently
rerun or authenticate the biological assertions inside the source files.
`verify_sources=False` is an explicit trusted-fixture mode and is reported as
unverified; score admission never uses that mode.

Required sites need agreement between the registered independent mapping
methods, coherent order, valid exact-sequence positions and a qualified domain.
Sites cannot be assembled from different domain instances. Generic nuclease
domains do not establish a reference subfamily without separate qualified
reference evidence. RNA/unknown substrate classifications do not become DNA
nuclease evidence. Source-bound modules require explicit native caller/system
identity and distinct member identities; module support never transfers
catalytic assignment to every component. Local structural confidence is not
active-site geometry. Missing evidence remains unassessed, distinct from
incomplete or conflicting evidence and from biological inactivity.

## Optional registry-to-policy admission

`analysis.admission.apply_reference_anchors(evidence_records, candidates, registry,
admission_rules, registry_dir=..., source_root=..., policy=...)` returns updated
records and an audit without mutating inputs. Joining requires exact sample,
gene and `provenance.sequence_sha256` agreement. Candidate `gene_id` may be used
as an alias for `query_id` only when both agree. Candidates outside the policy
pool are reported, not injected.

Every admission rule has exactly these fields:

```json
{
  "family_id": "Registered nuclease family",
  "RM_families": ["Type II"],
  "roles": ["R"],
  "entry_sha256": "64 lowercase hex characters",
  "model_accession": "Exact accession.version",
  "model_sha256": "64 lowercase hex characters",
  "reference_accession": "Exact registered reference",
  "required_anchor_subset": [2, 4],
  "source": "local_rule_provenance.txt",
  "source_sha256": "64 lowercase hex characters"
}
```

Use `analysis.admission.entry_sha256(entry)` for the canonical entry digest.
The example positions are illustrative, not a biological rule. Sample/locus
selectors and additional hidden fields are rejected. The family, substrate,
role, exact qualified model, reference and entire required anchor subset must
all match. A coherent source-verified domain yields its own linked Profile with
coverage computed from exact intervals and sequence length. The new anchor
cannot borrow an unrelated passing profile. The normal DNMB policy still decides
whether motif/context/retention gates are satisfied; no score, weight, threshold,
structure result or activity label is invented by this bridge.

Registries, every rule and every supplied candidate are validated even when the
candidate batch is empty or none join the policy pool. Actual external reference
FASTA file hashes are recorded in addition to sequence hashes. Reapplying the same
extension is idempotent: existing bound profile indices are verified and reused.
Changed or conflicting previously admitted sources require original unextended
policy records; duplicate profiles are never silently appended by the batch bridge.
