# Install once, analyze result parents

Run the commands in the [README](../README.md#first-install) from a new source
checkout. `defenseviz setup` is an installation/resource operation; it does not
analyze genomes. `defenseviz from-suite PARENT --out RESULTS` performs the full
R–M analysis and all exports once setup has completed.

## What is installed

The Docker image supplies the compatible Python, R libraries, BLAST, HMMER and
CATH runtime. Its DNMBsuite base image is pinned by digest. The setup receipt
records the actual source files, Git commit and working-tree state, final image
ID, commands, elapsed times and exit codes. A changed source checkout does not
silently reuse an image built from different bytes. Reinstall the Python package
from the same checkout before setting up an updated runtime.

Only the downstream references needed by DefenseViz are downloaded:

- REBASE's complete `protein_seqs.txt`, from the [official sequence download
  page](https://rebase.neb.com/rebase/rebase.seqs.html). The provider currently
  distributes it via anonymous FTP. Setup records its acquisition date and
  SHA-256; its release remains unknown unless the provider supplies one. FTP
  does not authenticate the transfer, and a newly obtained snapshot is not
  claimed to match an older database. Provider terms remain applicable.
- [MicrobeMod database v3](https://doi.org/10.6084/m9.figshare.24456481.v3), through
  the official Figshare API. Its fixed file ID, size and provider MD5 are checked;
  the receipt also records SHA-256 and the provider's license. Only the R–M HMM,
  profile metadata and exclusion list are extracted. DefenseViz does not run
  MicrobeMod or change these models and cutoffs.

InterProScan, eggNOG, DefenseFinder, PADLOC and DefensePredictor observations come
from the completed DNMB input folders. Their large upstream databases are not
downloaded by DefenseViz setup. Run DNMB on the genomes first when those results
do not yet exist.

## Persistent cache

The default is `~/.cache/defenseviz`, or `defenseviz` within `XDG_CACHE_HOME`.
For a different disk, set the same variable before setup and every later analysis:

```sh
export DEFENSEVIZ_CACHE=/path/to/persistent/defenseviz-cache
defenseviz setup --source /path/to/DefenseViz
defenseviz from-suite /path/to/completed-DNMB-parent --out /path/to/new-results
```

`--cache` is also available on setup, but subsequent commands must use the same
`DEFENSEVIZ_CACHE` location. The ready runtime profile is published only after
provider integrity, database parsing, indexes and native runtime checks pass.
Incomplete downloads and failed setup logs remain separate from ready snapshots.
Interrupted setup can be retried; it does not activate a partial database.

Repeated setup with the same source validates cached references and the image
without provider requests or rebuilding indexes. Analysis mounts databases and
inputs read-only and has no container network access. It runs new searches for
the requested samples; cached databases do not mean cached sample results.

An explicit `--refresh-databases` downloads new provider snapshots and prepares a
new resource directory. A failed refresh leaves the previous ready profile
available. Existing snapshots are not silently overwritten or removed.

## Analysis and logs

```text
completed-DNMB-parent/
  genome-one/ ...native DNMB outputs...
  genome-two/ ...native DNMB outputs...
```

Use `defenseviz from-suite completed-DNMB-parent --out new-results`. Folder names
are arbitrary; normal native artifacts and assembly metadata identify samples.
When a native total workbook exists it supplies the complete gene annotations,
including optional modules; the smaller REBASE input table is used only when no
total workbook is available. Ambiguous total workbooks require an explicit gene
input selection, and malformed selected inputs fail rather than being guessed.

The standard output directory contains `tables/`, `figures/publication.pdf`,
three-page exports when those panels are applicable, and `receipt.json`.
A sibling hidden `.RESULTS.defenseviz-*` directory retains the exact Docker
command, stdout/stderr, runtime profile and execution receipt, including on
failure. Setup logs live under the cache's `setups/` directory. Data and scientific
results are never published by setup.

The managed command currently supports Linux/macOS POSIX input paths; use WSL
paths on Windows. Spaces and Unicode are supported. Docker's mount syntax does
not support commas/newlines in these paths. Keep outputs outside the input parent.
Simple project JSON files and explicit input selections are supported; projects
with external native execution settings require an explicit native engine or
the documented direct container workflow. Explicit CLI, project, environment
and user-config engine settings take precedence over the managed profile.

Docker must be installed and running. A network or provider failure leaves setup
incomplete with the failing command and logs; retry after restoring access.
Firewall restrictions on official FTP must be resolved for automatic REBASE
acquisition. Alternatively, an authorized local provider snapshot can be prepared
with [prepare-resources](resource-preparation.md) and an explicit native/container
configuration. TLS verification is never disabled.

The default renderer uses a font available in the managed image. It does not copy
or redistribute proprietary workstation fonts, so layout remains reproducible
within the recorded runtime but is not promised pixel-identical to an earlier
Arial installation. Analysis output remains annotation evidence rather than
proof of biological activity.
