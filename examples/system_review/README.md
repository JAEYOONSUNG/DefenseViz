# Reviewing a previously unknown defense family

This example creates three **entirely synthetic** DNMBsuite result folders. The
protein sequences, model names, and rule are invented software fixtures. They
are not a biological rule set or evidence of activity.

Generate and verify the example from the repository root:

```bash
python examples/system_review/make_example.py /tmp/defenseviz-system-review-inputs --verify
```

The output directory must be new. The verification prints
`SYSTEM_REVIEW_EXAMPLE_PASS` and writes `review-before.json`, `review-after.json`,
and `verification.json` with exact input hashes.

| Input | Without an extension rule | With this exact rule |
|---|---|---|
| Exact model and model-package version | Role unassessed | Component role supported by the synthetic annotation rule |
| Same model, unavailable model-package version | Role unassessed | Role unassessed; missing version queued |
| An unregistered model | Role unassessed | Role unassessed; review required |

All three source detections remain accounted for. Motif and activity support
remain unassessed. The accepted effector, associated, motif, and priority counts
remain unchanged: this extension adds a review result, not a new biological
adjudication or a new detector run.

To render these native saved inputs through the regular pipeline:

```bash
defenseviz from-suite /tmp/defenseviz-system-review-inputs/runs \
  --manifest /tmp/defenseviz-system-review-inputs/cohort.json \
  --system-review-rules /tmp/defenseviz-system-review-inputs/review-rules.json \
  --out /tmp/defenseviz-system-review-figures
```

The source table explicitly records `model_version` when available. It must
match `match.model_versions` exactly; a software version or database release
cannot stand in for a model-package version. Native artifact and current
protein hashes bind the observation. The rule explicitly declares host-defense
scope; an anti-defense namespace requires its own scoped rule.

For a real new family, replace the invented rule with a reviewed rule that has
primary literature, an exact producer/model scope, and a claim limited to the
available evidence. Sequence-domain and catalytic-site conclusions need their
separate qualified evidence contracts. A descriptive name, transferred eggNOG
PFAM label, or neighboring nuclease does not establish those conclusions.
