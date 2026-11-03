# Synthetic source-bound reference extension

All sequences, model identifiers and evidence here are artificial test data. They demonstrate contracts; they do not assert a real enzyme or biological activity.

From the repository root:

```sh
defenseviz validate examples/reference_extension/project.json
defenseviz run examples/reference_extension/project.json --out extension-output
defenseviz assess --registry examples/reference_extension/registry.json --candidates examples/reference_extension/candidates.json --out assessment-output
```

The extension admits a coherent sequence-bound domain and two independently recorded site mappings into the existing motif-score axis. Model/reference/entry/source hashes must match. Sequence, profile and context requirements still apply; a mapped residue alone does not assign motif credit. Alter the reference or source file without updating its declared evidence provenance to see validation fail.

For real data, provide independently justified family rules, official model identities and qualified domain hits, exact reference sequences/site provenance, and source-backed coordinate mappings. This command evaluates saved evidence; it does not generate those searches or maps. No real databases or manuscript samples are bundled.
