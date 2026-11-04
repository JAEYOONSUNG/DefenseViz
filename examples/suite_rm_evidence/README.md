# Reevaluate stored R-M evidence

This example creates three entirely invented DNMBsuite output folders. Their
legacy REBASE exports have no saved retention column. A separate, fully bound
Evidence file supplies every required axis to the unchanged DNMB policy.

```sh
python examples/suite_rm_evidence/make_example.py demo-inputs
defenseviz from-suite demo-inputs/runs --manifest demo-inputs/cohort.json \
  --rm-evidence demo-inputs/rm-evidence.json --out demo-results
```

The workflow renders the manuscript atlas, type tables and comparison pages,
with editable SVG and policy decisions in `rm_evaluation.json`. Original tool
calls remain unchanged. This synthetic fixture is a software demonstration;
its invented evidence is not biological validation. Do not fill missing real
assessments with the example's boolean values.
