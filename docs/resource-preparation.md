# Prepare the analysis databases

For a completely new installation, use [`defenseviz setup`](first-install.md).
It downloads the supported public provider files, prepares these resources in a
persistent user cache and configures the managed runtime automatically. The
commands below are the advanced path for explicit local snapshots or existing
native/container installations.

A source checkout and DNMBsuite result folders do not include the analysis
databases. Use an existing validated [resource bundle](runtime-bundle.md), or
use the installed `defenseviz prepare-resources` command with authorized provider
files as described below. Keep the
provider files, versions, licenses and hashes beside the prepared resources.
Preparation is a separate, one-time step; analysis then accepts one result
parent directory and produces the tables and figures in one command.

```text
validated-databases/
  rebase/rebase_data.rds
  rebase/rebase_db.fasta
  rebase/rebase_db.fasta.*       # complete BLAST index
  rm-hmm/RM_HMMs.hmm
  rm-hmm/profile_roles.tsv
  rm-hmm/preparation.json
```

## R-M HMM library and role metadata

The supported provider format is the **MicrobeMod database**, with these three
matching source files: `HMMs/RM_HMMs.hmm`, `restriction_metadata.csv` and
`HMMs/off_target.txt`. The v3 provider release is
[MicrobeMod Database v3, DOI 10.6084/m9.figshare.24456481.v3](https://doi.org/10.6084/m9.figshare.24456481.v3).
Its `microbemod_db_v3.tar.gz` file has Figshare ID `55209713`, size `195132093`
bytes and provider MD5 `279f19f6d204fc4e1b750aea9828151d`. Select that version
when reproducing that resource set; downloading a later release is a database
change. The preparation script does not download anything or execute MicrobeMod.

From a DefenseViz source checkout, point to the directory containing the
extracted `HMMs` directory and CSV:

```sh
HMM_SOURCE=/absolute/path/to/provider-db
RESOURCES=/absolute/path/to/validated-databases
python scripts/prepare_hmm_resources.py \
  --hmm "$HMM_SOURCE/HMMs/RM_HMMs.hmm" \
  --metadata "$HMM_SOURCE/restriction_metadata.csv" \
  --off-target "$HMM_SOURCE/HMMs/off_target.txt" \
  --out "$RESOURCES/rm-hmm"
```

The output directory must be new. This standard-library Python script copies
the HMM byte-for-byte, retains its native GA cutoffs, and writes a paired TSV
plus an input/output SHA-256 receipt. The receipt lists unused provider metadata
and exclusion names, so a library subset is explicit. No sample data is needed.

`profile_roles.tsv` contains all six fields consumed by the analysis:

| Field | Source |
|---|---|
| `profile_name` | Exact HMM `NAME`, joined to CSV `Name` |
| `roles` | CSV `Enzyme_type`: `MT` → `M`, `RE` → `R`, `SP` → `S`, `IIG` → `M;R` |
| `type_hints` | Explicit CSV `System`, removing only `RM_Type_` or `Type_` |
| `off_target` | Exact membership in the supplied `off_target.txt` |
| `profile_length` | Exact HMM `LENG` |
| `library_sha256` | SHA-256 of the complete unchanged HMM file |

The supported type values are I, II, III, IV, IIG and IIM. Off-target profiles
receive no R/M/S assignment. Missing, unknown or conflicting metadata for a
nonexcluded profile fails preparation; a model's name is not used to guess its
role. IIG is an explicit family-level dual-role hint, not evidence that a query
protein is biologically active. Extra provenance fields retain source row
numbers, original enzyme/system labels and native GA values.

For a container-only installation, use the
[DNMBsuite image recipe](dnmbsuite-container.md) and run the same preparer:

```sh
mkdir -p "$RESOURCES"
docker run --rm --platform linux/amd64 \
  --mount "type=bind,src=$HMM_SOURCE,dst=/source,readonly" \
  --mount "type=bind,src=$RESOURCES,dst=/resources" \
  --entrypoint /opt/defenseviz/venv/bin/python \
  dnmbsuite-defenseviz:local /opt/defenseviz/source/scripts/prepare_hmm_resources.py \
  --hmm /source/HMMs/RM_HMMs.hmm \
  --metadata /source/restriction_metadata.csv \
  --off-target /source/HMMs/off_target.txt --out /resources/rm-hmm
```

The resource destination is writable during preparation. Mount it read-only
for subsequent analysis.

## REBASE metadata, FASTA and index

A DNMBsuite shared cache normally resides at
`$DNMB_CACHE_ROOT/db_modules/rebasefinder/cache`. In the standard suite container,
`DNMB_CACHE_ROOT` defaults to `/opt/dnmb-cache`; use its corresponding host mount.
The per-genome result folder is not this shared reference cache.

The installed `defenseviz prepare-resources` command prepares both databases and
an engine configuration together. It needs the native runtime described in
[runtime-install.md](runtime-install.md), including the pinned R helper namespace,
BLAST and HMMER. The HMMER installation must provide `hmmstat` beside the
configured `hmmsearch`; setup uses it to parse every complete model without
query sequences. A Python installation alone does not supply that runtime.

To preserve the exact databases from an existing engine configuration:

```sh
defenseviz prepare-resources --engine-config existing-engine.json \
  --out /project/validated-databases
```

The destination must be new. The command copies the complete REBASE metadata,
FASTA and index set, preserves available raw/reference-state files, and copies
the paired HMM and role metadata. It verifies reference metadata against the
FASTA, the BLAST index against the FASTA, and HMM names, lengths and library hash
against the role table. HMMER must also successfully parse every model body;
valid headers alone are insufficient. Existing source resources remain unchanged. It writes
`engine.json`, `preparation.json` and validation logs in the new directory. Runtime
paths are resolved from the original configuration, so moving the new config
does not reinterpret its relative source paths.

For a new resource set, obtain an authorized REBASE `protein_seqs.txt` snapshot
and the matching HMM provider files described above. Use a runtime configuration
with the native tool and R library paths; its database paths will be replaced:

```sh
defenseviz prepare-resources --engine-config runtime-engine.json \
  --rebase-source /source/REBASE_protein_seqs.txt \
  --hmm /source/HMMs/RM_HMMs.hmm \
  --hmm-metadata /source/restriction_metadata.csv \
  --off-target /source/HMMs/off_target.txt \
  --source-label 'Provider snapshots and dates recorded in the accompanying provenance' \
  --out /project/validated-databases
```

Use `--makeblastdb /path/to/makeblastdb` if it is not available beside the
configured `blastp` executable. `--reference-cache` selects a complete existing REBASE cache;
`--hmm` together with `--profile-metadata` selects an already paired HMM set.
The default uses the database paths in the source engine configuration. Do not
mix REBASE metadata and FASTA from different snapshots.

Raw REBASE preparation uses the pinned legacy R helper's native parser and FASTA
writer. The historical R namespace is an actual dependency API, not another
analysis product. The command records input/output hashes, source label, commands,
versions, exit status and elapsed time. It does not download references or run
searches against genome queries. An unspecified provider release stays unknown;
file modification times are not evidence of a provider release date. A newly
acquired snapshot does not reproduce an older unversioned reference.

For the DNMBsuite container, mount provider files read-only, mount a new resource
parent writable, and invoke `/opt/defenseviz/venv/bin/defenseviz prepare-resources`
with `/opt/defenseviz/container/engine.json` and the corresponding mounted paths.
The generated config references the container runtime. Mount the prepared
resource directory read-only for subsequent analysis.

## Validate, then analyze

Use the generated `/project/validated-databases/engine.json` for analysis.
Preparation already validates the assembled resources. To check them again later,
use `defenseviz check-resources --engine-config /project/validated-databases/engine.json
--out resource-check001`.
Preserve its report and input hashes. Then run
`defenseviz from-suite /path/to/suite-parent
--engine-config /project/validated-databases/engine.json --evaluation full --out results/run001`, or use the complete
[Docker command](dnmbsuite-container.md#run-the-full-analysis-with-one-parent-folder).

Each check/analysis destination must be new. A setup or resource-check receipt
is not an analysis result; only the full command generates the evaluated CDS
lists, workbook and figures.
