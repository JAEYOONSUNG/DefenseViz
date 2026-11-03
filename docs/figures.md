# Figure 2 panels a and b

`defenseviz run project.json --out output` writes a vector PDF, editable SVG pages and `figure2_layout.json`. A single page is `figure2_ab.svg`; multi-page exports use `figure2_ab_page001.svg` and so on. PDF pages follow the same order. SVG text stays text; no raster images are embedded. Text can be edited in Illustrator or Inkscape. The source font must be installed when editing; the layout JSON records its family name.

Panel **a** shows the actual saved-input/evidence/evaluation/reporting stages. Support counts come from evaluated records. The displayed threshold and required capped axes come from the selected policy. Parameter overrides are labeled `DNMB policy (custom)` when the analysis metadata identifies a modified preset. A zero supported-evidence count does not mean that a fresh search or structure prediction was performed and found zero results. No detector or structural engine runs in this renderer.

Panel **b** contains:

- Native caller bars: distinct CDS identifiers in each caller's final selected output. Raw HMM hits do not become native calls.
- **Union**: the analysis layer's integrated union, including retained DNMB-only CDSs.
- **RM pool**: explicit R-M candidate loci from native/reference evidence and the evaluated evidence input.
- **Retained**: CDSs passing the selected DNMB evidence policy. **Unass.**: R-M candidates without an evaluable decision.
- R-M component columns: candidate CDS counts for explicit type/role pairs. A ring means at least one of those CDSs has qualified support in the evaluated evidence; it does not assert a complete system or a biological activity assay. Exact supported counts are in the result tables and layout JSON.
- Family matrix: distinct CDSs in **native final caller output**, colored by the supplied category mapping. Novel families remain visible. A locus may appear in more than one family; these columns cannot be summed to reconstruct the integrated union.

Native system counts remain separately available in the machine-readable analysis tables. The figure never treats a pooled R/M/S combination as a complete biological system.

A small gray dot is zero observed native calls in the available outputs. `NA` is missing/not-run caller coverage; unknown output is not converted to absence. Partial caller coverage can therefore yield an observed partial count. The analysis result preserves coverage and count status. These displays describe annotation, not genome-wide biological truth.

## Layout and editing

Samples keep the configured order. `highlight: true` produces pink shading; there are no built-in strain identities. Caller and summary bars use one physical thickness, fixed padding and a shared full-cohort maximum within each section. Their scales do not change between pages. Numeric bar lengths are never forced equal.

Rows split at `max_rows_per_page` (default 40) or a physical row-height limit. Long sample labels and figure titles wrap according to actual font metrics. Multi-line titles expand the reserved header area, keeping the page number, workflow, table and labels separate. Family columns split at `max_families_per_page` (default 24). Extremely many caller or R-M type/role columns use additional column pages. Each sample/family combination is retained, and a/b headings repeat for context. The PDF's physical page dimensions adapt to the required columns; exports are not forced into an unreadably reduced Letter page.

The layout JSON records page dimensions, measured text bounds, row/cell bounds, page membership, counts, missing-data states and exact bar geometry. No page-bound clipping is accepted. Unicode labels are retained. If no installed font supports a label, the renderer reports an actionable error instead of silently emitting missing-glyph boxes; install a suitable Noto Sans/CJK or Unicode font.

The renderer refuses to overwrite existing figure outputs. The CLI runs rendering in an atomic staging directory and publishes only successful results.
