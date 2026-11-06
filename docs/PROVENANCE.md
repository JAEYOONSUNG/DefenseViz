# Source lineage and scientific scope

DefenseViz is a Python implementation of defense annotation integration and plotting, informed by two historical tools by Jae-Yoon Sung:

- DNMB `R/DefenseFinder_plot.R`: the original `DefenseFinder_Heatmap()`/DefenseFinderViz visualization entry point (MIT).
- The earlier DefenseViz R restriction–modification analysis package (MIT). The current Python CLI has a separate interface and release history.

The numerical DNMB policy derives from the source-bound revision027 policy and the selected revision026 parameter preset. Original module names and SHA-256 identities are recorded in the analysis modules and release validation. Biological rule defaults are versioned configuration; dataset identifiers, manuscript membership and desired output counts do not enter the decision predicate.

The preset's historical calibration provenance is disclosed. Its parameters are not claimed to have independent population-level validation. This migration does not optimize thresholds against the manuscript cohort. A parity check establishes software behavior reproduction, not biological correctness or classifier accuracy.

The evidence catalog evaluator keeps family, domain, catalytic-site, module and structural outcomes separate. External catalog/reference/model files remain user-supplied and checksum-bound; raw REBASE/Pfam databases and unpublished manuscript data are not redistributed. Source reference correspondence is evidence of annotation, not proof of activity. Old module names may appear in provenance without representing separate public DNMB products.

Full mode executes DNMB R–M HMM and REBASE searches with the configured runtime and databases. Stored mode replays saved decisions. DefenseFinder, PADLOC and DefensePredictor annotations remain saved inputs; structure prediction is not rerun. A supplied policy feature such as prior structural support records an upstream assessment. Coverage gaps and unsupported input fields remain visible.
