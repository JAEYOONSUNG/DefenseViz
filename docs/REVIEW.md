# Independent critique and resolution log

The analysis, IO and rendering components were implemented with separate ownership and then cross-reviewed by developers who did not write the reviewed components. Root reruns the integrated tests and verifies real stored data independently.

## Scientific and data-contract findings resolved

- Missing caller output could be mistaken for zero selected calls: final-member coverage now explicit; HMM/reference records do not imply final output completeness.
- Reference-only providers could appear as native callers: auxiliary evidence providers are separate, without name-based exceptions.
- Native genes without system IDs could yield false zero system counts: unknown counts are null and explicitly scoped.
- Control C/Aux roles could enter R/M/S scoring: preserved in the R-M pool but outside score scope.
- Known gene/evidence sequence hashes could conflict: incompatible bindings now fail; incomplete evidence remains unassessed.
- Empty or unmatched extensions could bypass registry/rule validation: all configured contracts and source bindings are checked before filtering candidates.
- Reference FASTA raw input hashes were omitted: recorded and stability-checked alongside normalized sequence identities.
- Repeated anchor enrichment could duplicate profiles: exact source-bound reapplication is idempotent.
- Native feature sequence hashes did not reach admission: existing hashes are lifted through explicit adapters, with disagreement rejected.
- Modified policy parameters retained preset labeling: custom overrides now carry parameter differences and distinct status.
- Default display categories could overwrite explicit input classification: explicit input/user-map precedence is enforced.

## Engineering and figure findings resolved

- Staging absolute paths broke relocation equivalence: generated provenance uses config-relative presentation pointers, with actual file hashes retained in receipts.
- Excel extra documentation differed from packaging: standardized on `.[excel]`.
- Caller headers, long titles and rotated family labels could overlap: measured wrapping and dynamic layout bounds replace fixed assumptions.
- Dollar signs could activate math parsing in identifiers: text is literal, including malformed math-like strings.

## Operational scope

Use a fresh output target for each run. Do not let unrelated programs concurrently write the same destination directory; such external writers are outside the output transaction contract. Clean installation tests use a new virtual environment and a built wheel. No source detector, structure prediction or recognition-specificity code is run by this package. Saved structural/context facts remain upstream inputs; their raw computation is not reconstructed from a bare caller table.
