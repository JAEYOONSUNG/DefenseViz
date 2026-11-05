# Install the native runtime from public source

DefenseViz can install its required R namespaces directly from pinned GitHub
source. A previous benchmark workspace or local patched source archive is not
required by this installer. Databases, compatible third-party R dependencies and
original DNMBsuite outputs remain separate inputs. For a new installation,
prepare a source-bound database directory using [resource-preparation.md](resource-preparation.md).
The [DNMBsuite container](dnmbsuite-container.md) already includes the compatible
Python, R and native tool environment described below.

Install the Python package first, then run the setup script from its source
checkout:

```sh
python -m pip install '.[suite]'
python scripts/install_native_runtime.py \
  --rscript /absolute/path/to/Rscript \
  --dependency-library /existing/DNMBsuite/R/library \
  --out /project/runtime/native-r
```

Repeat `--dependency-library` for each compatible dependency library needed by
your DNMBsuite environment. The destination must be new and outside those
libraries. The installer never changes existing packages or automatically
downloads an unpinned dependency collection. Missing R dependencies are named in
the failure message and recorded logs. Use matching R/platform dependency builds.

The fixed public source pins are:

| R namespace | Repository | Commit |
|---|---|---|
| DNMB | `JAEYOONSUNG/DNMB` | `976489e9967e9e4dfacc61f7e83e66b1da0f4f2b` |
| Legacy R DefenseViz | `JAEYOONSUNG/DefenseViz` | `53a1dc4b8c6595945762e21354fbdc86d2b3731a` |

Both carry MIT license text. The legacy R namespace is a pinned dependency of
the current Python package; it is not a second public analysis product. The
current fixed annotation/context/motif/scoring modules are supplied by the
Python package and retain their own source hashes.

The installer uses sparse Git checkouts of installable code, documentation,
licenses and audited reference assets. It excludes historical sample workbooks,
GenBank/FASTA examples and validation/output directories. It verifies the actual
commit, source working-tree state and selected file hashes, then installs the
two packages into `native-r/library`. Dependency origins and versions, source
licenses, commands, logs, exits and elapsed times are retained in the new
directory. The 19 inherited DNMB helper functions are compared structurally with
the pinned source, including function arguments and bodies; loading an older
DNMB namespace from a dependency directory is rejected.

`runtime_installation.json` supplies the exact `r_libraries` list to use in a new
engine configuration: the newly installed library first, followed by the
explicit dependency libraries. Preserve the resource paths from the database
configuration, resolving relative paths against that configuration's original
directory when writing the new configuration elsewhere. See
[engine.example.json](engine.example.json) for the supported resource fields.

Validate the resulting configuration, then analyze original DNMBsuite results:

```sh
defenseviz check-resources --engine-config /project/runtime/engine.json \
  --out /project/results/resource-check-001
defenseviz from-suite /project/inputs/DNMBsuite-results \
  --engine-config /project/runtime/engine.json --evaluation full \
  --out /project/results/run-001
```

Successful installation and resource checking do not claim a completed analysis.
The full command performs fresh R-M HMM/reference searches and produces the
evaluated annotation lists, tables and figures. Original DF/PAD/DP annotations
remain saved observations; this does not rerun those tools or InterPro/EggNOG.

REBASE and HMM databases are not downloaded by this source installer. Use the
[resource preparation commands](resource-preparation.md) to prepare the exact
databases and profile metadata in your resource configuration,
with their license/version provenance. A source clone does not reproduce a
historical, unversioned reference snapshot. The separate
[runtime bundle exporter](runtime-bundle.md) preserves exact local resources and
an installed-package fallback, without requiring those databases in GitHub.
