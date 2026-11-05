# Preserve a relocatable DNMB runtime

DefenseViz source code and local analysis assets have separate lifecycles. Keep a
single runtime/database bundle, original DNMBsuite results, and new analysis
outputs in separate directories:

```text
project/
  runtime/engine/     # preserved DBs, R packages, source packages, engine.json
  inputs/            # original DNMBsuite results; read-only
  results/run-001/    # new figures, tables and execution evidence
```

The exporter below creates the first directory without copying sample data or
changing the originals. It requires an installed DefenseViz and a checkout of
this repository for the script. The bundle is a local preservation artifact;
database and third-party package licenses still apply. Do not add it to Git.

## Export the configured environment

Provide the actual source packages corresponding to the installed DNMB and
legacy R DefenseViz runtime, including any local patches. These R source
packages are distinct from this Python repository. The exporter records their
content and Git metadata without assuming that a version number identifies all
local changes.

```sh
python scripts/export_runtime_bundle.py export \
  --engine-config /existing/engine.json \
  --dnmb-source /preserved/source/DNMB \
  --legacy-source /preserved/source/DefenseViz \
  --evidence /preserved/source-patches-and-receipts \
  --evidence /preserved/database-provenance \
  --out /project/runtime/engine
```

`--evidence` may be repeated for files or directories containing patches, build
receipts, database versions and license provenance. The destination must be new.
An interrupted/failed export remains diagnostic evidence and must not be used as
a completed bundle.

The export contains:

- One canonical REBASE metadata/FASTA/index set, plus the original raw reference
  and reference-state file when present.
- One configured R-M HMM database and its exact profile-role mapping.
- The measured R namespaces and their recursive `Depends`, `Imports` and
  `LinkingTo` dependencies, excluding base/recommended packages supplied by R.
- Exact installable DNMB and legacy R DefenseViz source files and an executable
  rebuild script. Installed local package code is retained as a fallback.
  Historical sample workbooks, FASTA/GenBank examples, validation outputs and
  root run artifacts are excluded from both source and fallback packages; their
  omission hashes are recorded. Audited DNMB structure references/templates
  remain package assets. Other third-party package resources retain their
  original package content and license.
- A relative CATH wrapper, its configured binary and runtime libraries. Shared
  library alias filenames are preserved.
- A relative-resource `engine.json`, exact file hashes, original source bindings,
  package licenses, platform information and system executable hashes.

R, BLAST+, `blastdbcmd` and HMMER remain explicit system dependencies; their
resolved executable paths/hashes are recorded. The copied R packages and
CATH/ICU binaries are platform-specific. Moving the bundle on the same compatible
machine is supported; a different operating system, architecture or R version
requires rebuilding dependencies and a new validated runtime configuration.

## Verify after relocation

Move the complete bundle without editing its files, then verify from its new
location. Write the verification report outside the preserved bundle.

```sh
python scripts/export_runtime_bundle.py verify \
  --bundle /new-location/runtime/engine \
  --report /new-location/results/runtime-check-001
```

Verification checks every file against the bundle manifest, loads the actual R
runtime, rejects packages loaded from an unpreserved user/workspace library,
executes CATH and checks every BLAST index entry against the exact reference
FASTA. This is a runtime/resource check, not a full analysis.

Then run the full pipeline against original inputs:

```sh
defenseviz from-suite /project/inputs \
  --engine-config /new-location/runtime/engine/engine.json \
  --evaluation full \
  --out /new-location/results/run-001
```

For a manifest-selected cohort, also provide `--manifest cohort.json`. The
analysis performs fresh R-M HMM/reference searches and downstream evaluation;
saved DF/PAD/DP and independent annotation inputs remain the original observations.
Compare gene/system counts, complete candidate/search inventories, tables and
figures to the accepted run before considering removal of old intermediate
directories. A runtime preflight alone is insufficient for that decision.

## Rebuild the exact R source packages

The installed package copy is not the only recovery route. The bundle includes
an executable installer for its preserved source trees:

```sh
python /project/runtime/engine/setup/rebuild_r.py \
  --library /project/runtime/rebuilt-r-library
```

The new library must be outside the preserved bundle and must not already exist.
The script verifies bundle hashes, copies the preserved third-party dependency
packages, and executes `R CMD INSTALL` for legacy R DefenseViz followed by DNMB.
It writes exact commands, exits, logs and elapsed times to
`rebuild_receipt.json`. Supply `--rscript /absolute/path/to/Rscript` only when
intentionally selecting a different installation.

Use the rebuilt library in a new engine configuration. Resolve the original
configuration's relative database/tool paths against the bundle first; do not
copy relative strings to an unrelated directory. Retain the original bundle and
manifest. Run the full analysis and compare outcomes before adopting the rebuilt
environment. This script rebuilds the two preserved source packages; it does not
download or rebuild every third-party dependency for a different platform.

## What a GitHub clone does and does not provide

The Python code, fixed native evaluation modules and exporter can be versioned
in GitHub. Exact sample annotations, local source revisions, reference snapshots
and externally distributed HMMs are separate inputs. The REBASE snapshot used in
an earlier analysis may have no immutable public release identifier, so fetching
the current database is not necessarily equivalent. An HMM database may be
ignored by its upstream code repository. Preserve their exact files and supplied
license/version records until a clean recovery reproduces the hashes.

The bundle manifest proves preserved content identity, not biological activity,
independent accuracy or permission to redistribute third-party data. Historical
manual reviews, immutable evaluation evidence and final publication outputs
also require their own preservation. This exporter deletes nothing.

## Check resources through the installed CLI

Once `engine.json` points to the prepared resources, the installed package can
check them without a source checkout:

```sh
defenseviz check-resources --engine-config /project/runtime/engine/engine.json \
  --out /project/results/resource-check
```

This checks reference metadata and every indexed protein, the HMM gathering
thresholds and exact profile-role mapping, frozen evaluation modules, executable
versions and loaded R packages. It writes `resource_check.json` and resource
hashes; failure is explicit. It does not download data or perform an analysis.

Run `from-suite --evaluation full` after the check. The pipeline writes the
classified Excel workbook, retained R-M and prioritized annotation lists, and
publication PDF/SVG pages in that single invocation. R-M retention and broader
annotation-based prioritization remain separate output sheets.
