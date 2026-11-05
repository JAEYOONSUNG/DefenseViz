# System reference information

DefenseViz includes an offline, source-pinned catalogue of the official
[DefenseFinder wiki](https://defensefinder.mdmlab.fr/wiki/defense-systems),
[DefenseFinder model definitions](https://github.com/mdmparis/defense-finder-models),
and the separately maintained [CasFinder definitions](https://github.com/macsy-models/CasFinder).

```sh
defenseviz catalog --out all-system-information
defenseviz catalog Wadjet --out wadjet-information
```

Without `--out`, the command prints matching reference records as JSON. Search is
for browsing; a search hit does not assign a detected gene to that system. The
command works outside a source checkout and does not contact the internet.

The export includes system information, native model definitions, model
components, source links, and structured JSON. Wiki information includes the
reported sensor, activation mode, effector mechanism, PFAM identifiers and paper
references, with concise source-based mechanism summaries where available.
An explicit `Unknown` in the source is distinct from an absent field. Article
abstracts, full wiki prose, HMM profiles and executable model XML are not bundled.

A wiki page describes a system or family. A model definition can describe a
subtype, detection alternative, reference component collection or anti-defense
system. These record counts are reported separately. They must not be added to
produce a supposed number of unique biological systems.

## Information for detected systems

`from-suite` and normalized-project `run` automatically write:

- `system_information.json`: structured per-observation and per-unit links.
- `tables/system_information.tsv`: the system/unit information index.
- `tables/system_information_*.tsv`: detailed observation and source tables.

`inspect-suite` reports the information-matching coverage. Unmatched names and
ambiguous matches remain in the output. A saved model no longer present in the
current catalogue is not discarded from the analysis.

Identity resolution uses original provider, native family/subtype and explicit
model FQN information. A component HMM name is not a system-model identity.
Display aliases, lowercasing, punctuation changes and search similarity do not
silently resolve native identities. Exact source-model identity, family-level
context and package-version compatibility are separate outcomes.

An exact cross-provider label match to a wiki title or slug supplies descriptive
context only. It is explicitly labelled as a name coincidence between independent
vocabularies, not a reviewed cross-tool equivalence. Native namespace scope and
the matched catalogue definition's scope are reported separately. Display-family
labels in figures may pool several native systems; join information through
original detection/unit IDs, not through display names.

The catalogue describes the pinned source snapshot. An older or unversioned run
can receive current descriptive context, labelled with its version relationship;
that does not establish that the current definition produced the historical
call. A matching family name across DefenseFinder and PADLOC does not assign
DefenseFinder's detection rules to PADLOC.

`same_declared_package_version` means package-version agreement, not identity of
an executed XML definition. A DefenseFinder `models_version` field describes its
own model package; it cannot bind the separately versioned CasFinder package.
`native_model_version` is the first applicable report from canonical
`model_version`, `source_version.models_version`, then
`native_fields.model_version`; all reports and conflicts remain in
`native_model_versions_reported`. Current
component-name comparisons under an unknown or different version remain unbound
to the historical run. Empty producer fields remain missing, and malformed
optional metadata is retained with diagnostics rather than blocking publication.

In TSV files, read `sensor`, `activator`, `effector` and `pfam` with their `_state`
columns. `source_mechanism_statement_state` is an editorial description of the
source's statement, not a confidence score for biological activity. Matching
counts are also reported by provider; model and wiki links can overlap.

## Relationship to evidence review

Reference information is additive context. It does not change original calls,
display families, accepted CDS sets, effector counts, motifs, prioritization,
system-review evidence states, or publication figures. In particular, a model's
mandatory component is not automatically a nuclease, and a wiki statement about
a studied system does not validate the user's query protein.

Use the separate [source-bound review rules](SYSTEM_REVIEW.md) when there is
evidence supporting an actual query-level claim. Source entries, explicit
unknowns, snapshot commits and file hashes remain inspectable in the catalogue.
Updates are reviewed release inputs; analysis never silently downloads a newer
database.

## Source snapshot

The bundled snapshot contains 152 official wiki pages, 666 DefenseFinder/ADF/R-M
model definitions, 44 separate CasFinder definitions, and 6,970 declared model
component entries. It retains every page and definition, including unknown or
undocumented mechanisms. Source genomic examples are examples, not mandatory
architectures for every member of a family.

The wiki is pinned to `ee7647d8da3a6dd1750d233987f3ed8feb0fc0eb`,
DefenseFinder models to `afb0e5a8b466be53586b13266f5d38d98c3ac268`, and
CasFinder to `c7b02a90900e38a43a4f5e4f669f13c50dcf6491`. Every entry records its
original path, content hash and source link. File-tree and source-content
reconciliation checks coverage; it does not independently validate every source
claim. Citations remain the references listed on the original pages, with
source-quality issues retained where identified.

See [third-party attribution](../THIRD_PARTY_NOTICES.md) for original source
terms and the upstream model-license discrepancy. Attribution is also installed
under `defenseviz/data/notices/`. Updating these reviewed snapshots is an explicit
release change, not an automatic network operation during an analysis.
