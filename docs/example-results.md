# Reviewed benchmark figures

The README shows the reviewed 52-genome Geobacillus benchmark. These exports include study-specific context review and sample annotations. The B. coagulans DSM 1 case table is an external comparison and does not contribute to the 52-genome totals.

| Endpoint | CDSs |
|---|---:|
| Integrated annotations | 2,929 |
| Assigned R–M, types I–III | 342 |
| Assigned type IV | 151 |
| Assigned R–M, including type IV | 493 |
| Effector annotations | 497 |
| Prioritized annotations | 500 |

These sets overlap and should not be added together. In the figures, **assigned** means a CDS-level annotation; orphan components and incomplete contexts remain identified. It does not imply a complete or experimentally active system. Historical decision fields and the `retained_RM_annotations.tsv` filename retain their existing names for compatibility.

In the priority panel, “Active R” denotes DNMB-selected or reviewed restriction-role annotations; it is not a measurement of enzyme activity. The `retained_RM_annotations.tsv` export contains policy-selected CDSs, while reviewed additions remain in the integrated annotation list.

The snapshot includes four methylase CDSs reassigned during context review to BREX, DISARM and Hma. Their methylase evidence is preserved. These supplied corrections are part of the reviewed figures, not evidence that the current automatic engine reproduces every reviewed assignment. The synthetic demo is a separate dataset.

- [Two-page PDF](assets/benchmark52/publication.pdf)
- [Atlas SVG](assets/benchmark52/atlas.svg)
- [Comparison SVG](assets/benchmark52/comparison.svg)
- [R–M assignment flow SVG](assets/benchmark52/rm-assignments.svg)
- [Figure hashes and snapshot metadata](assets/benchmark52/manifest.json)

Only rendered figures and their metadata are distributed here. Reproducing a reviewed analysis requires its DNMBsuite outputs, runtime, databases, configuration and review annotations.
