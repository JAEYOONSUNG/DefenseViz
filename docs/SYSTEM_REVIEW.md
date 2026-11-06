# Review every observed system

`from-suite` automatically reviews every saved detection, including families
outside the bundled component catalogue. It preserves the original calls and
publication counts. The review reports what supports each interpretation and
what evidence is still needed.

```sh
defenseviz from-suite RESULTS --out reviewed-results
```

The default auto mode uses a configured DNMB engine when available; otherwise it
evaluates stored evidence. Use `--evaluation stored` to replay saved decisions.
Full mode performs DNMB R–M sequence searches before system review. Original
DefenseFinder, PADLOC and DefensePredictor calls remain saved inputs, and
recognition-specificity algorithms are unchanged.

## Outputs

`system_review.json` contains the complete review. The same records are available
as TSV files in `tables/`:

| File prefix `system_review_` | Contents |
| --- | --- |
| `summary.tsv` | Flat per-unit counts, assessment states and reasons; start here |
| `units.tsv` | Native systems and explicitly ungrouped observations |
| `gene_reviews.tsv` | Evidence axes for each observed coding sequence |
| `rule_evaluations.tsv` | Exact rules evaluated, outcomes, and reasons |
| `pending_queue.tsv` | Unresolved questions and the required evidence |
| `detection_accounting.tsv` | Original detection-to-review-unit mapping |
| `sample_coverage.tsv` | Sample-specific input and assessment scope |
| `quarantine.tsv` | Raw records that could not be bound safely |

In the flat summary, `members_with_*_evidence` counts evidence about the member
proteins, including evidence from other systems. The adjacent unit `*_status`
columns describe only that unit's own bound contract. These scopes are separate;
unit counts must not be summed as unique CDS counts. `raw_observation_count`
describes the observations preserved by the importer, not a new census of every
line in the source files.

A native system is scoped by sample, provider, run, replicon and system ID.
Matching family names or neighboring genes cannot merge systems. Partial
supplied group lists cannot hide other detections. Predictions without a system
ID remain ungrouped; HMM and reference matches remain separate from final calls.

## Read the outcomes

`review_execution: completed` means that the supplied observations were examined.
It does not mean that every function was supported or that input data were
exhaustive. An unfamiliar family receives a completed review with explicit
unassessed axes and a queue entry; it is not treated as a negative.

The review separates native detection, family annotation, component role,
nuclease domain, motif/site evidence, system architecture and complex core.
Evidence attached to a member protein is reported separately from an assessment
of its particular native system. A complete core in one system cannot complete
another system that happens to contain the same protein. Exact anti-defense
namespaces retain their own scope; unknown scope stays unresolved.
Computational support does not establish biological activity. A system can be
reported by a detector while its effector or full component set remains
unconfirmed. Different detectors may use different component quorums.

Existing accepted assignments retain their provenance. New extension-rule
conclusions are review proposals and do not silently change the accepted
effector, associated or prioritized sets. An adjudicated change can enter the
existing source-bound `--supplement` workflow. The figure remains a view of the
accepted evidence, with the new review ledger beside it.

## Add rules without changing package code

```sh
defenseviz from-suite RESULTS --system-review-rules rules.json --out reviewed-results
```

A rule registry supplies its schema version, ID, version and exact rules.
Native-component rules require an exact provider and model, current protein and
source bindings, and an allowed model/database release or model checksum. A bare family
name, description or substring cannot establish a function. Unknown fields and
unsupported predicates are rejected rather than ignored.

`match.model_versions`, `database_versions` and `software_versions` refer to
different provenance fields and never substitute for one another. A rule can
declare `defense_scope` as `host_defense`, `anti_defense` or `unresolved`; an
explicit incompatible native scope prevents support. Protein identity imported
through the suite's ID maps is labelled `importer_attached`, distinct from an
independently producer-reported sequence binding.

Component-set rules require an explicit source-pinned architecture contract and
the required components in the same native unit. They assess that declared
component set; they do not replay the full detector grammar or establish
individual catalysis for every member of a complex. Missing partners, ambiguous
multiplicity and conflicting evidence remain visible.

See [the synthetic example](../examples/system_review/README.md) for runnable
unknown-family, exact-rule and missing-binding cases. Its sequences, models and
claims are invented solely to demonstrate the interface.

## Use qualified domain and site evidence

Raw InterProScan signatures and eggNOG annotations are included as context.
InterProScan sequence identity, length and coordinates are checked where the
saved sequence permits it; unbound or conflicting rows are quarantined. The
import scope covers observed detections and supplied review members, with
excluded rows counted explicitly. An absent InterPro row is not evidence of
absence. An eggNOG-transferred Pfam label is not a query-domain alignment.

For exact domain/reference-site evaluation, supply an evidence pack:

```json
{
  "schema_version": 1,
  "registry": "registry.json",
  "candidates": "candidates.json",
  "source_root": "."
}
```

```sh
defenseviz from-suite RESULTS --system-review-evidence evidence-pack.json --out reviewed-results
```

Paths are relative to the pack. `source_root` is optional and otherwise uses the
registry directory. The registry and candidates use the existing
[`assess` evidence contract](evidence.md). Every candidate's sample/query and
whole-protein SHA-256 must match the imported suite gene. Actual evidence files
and external reference sequences are checked and included in the receipt.
Qualified domain and site conclusions stay separate from effector-role claims.

## Declare assessment scope per sample

A cohort manifest can retain legacy `options.assessment_coverage` and override
individual samples:

```json
{
  "assessment_coverage_by_sample": {
    "sample-a": {"effector": "complete", "motif": "complete", "associated": "complete"},
    "sample-b": {"effector": "partial", "motif": "unassessed", "associated": "partial"}
  }
}
```

Place this object under `options`. These are declarations of the supplied
assessment scope, not claims of exhaustive biological knowledge. Explicit sample
values take precedence over the legacy declaration; otherwise conservative
defaults use that sample's evidence. External-cohort labels do not change a
scientific outcome. Completing the automatic review never upgrades these
declarations or removes a justified partial-count marker.
