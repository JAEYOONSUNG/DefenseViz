# Prepare the analysis databases

A source checkout and DNMBsuite result folders do not include the analysis
databases. Use an existing validated [resource bundle](runtime-bundle.md), or
prepare a new resource set from authorized provider files below. Keep the
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

If the original DNMBsuite run prepared REBASE, its cache normally resides at
`$DNMB_CACHE_ROOT/db_modules/rebasefinder/cache`. In the standard suite container,
`DNMB_CACHE_ROOT` defaults to `/opt/dnmb-cache`; find the corresponding host cache
mount. The per-genome result folder is not this shared reference cache.

Preserve a complete existing set together:

```sh
REBASE_CACHE=/absolute/path/to/suite-cache/db_modules/rebasefinder/cache
mkdir -p "$RESOURCES"
mkdir "$RESOURCES/rebase"
cp -p "$REBASE_CACHE/rebase_data.rds" "$RESOURCES/rebase/"
cp -p "$REBASE_CACHE"/rebase_db.fasta* "$RESOURCES/rebase/"
```

Also preserve `REBASE_protein_seqs.txt` and `rebase_reference_state.rds` when
available. Do not combine metadata from one snapshot with another FASTA/index.
The full pipeline checks encoded reference IDs and exact protein sequences;
`check-resources` verifies the BLAST index against its FASTA.

If no complete cache exists, obtain an authorized REBASE `protein_seqs.txt`
snapshot and prepare a new directory. The pinned legacy R helper namespace
installed by [runtime-install.md](runtime-install.md) exports
`parse_rebase_sequences(rebase_file)` and
`write_fasta_for_blast(data, output_file, id_col, seq_col)`. The following is a
source-compatible preparation script; save it as `prepare_rebase.R`:

```r
a <- commandArgs(TRUE)
stopifnot(length(a) == 2L)
raw <- normalizePath(a[1], mustWork = TRUE)
out <- a[2]
stopifnot(!file.exists(out), !dir.exists(out))
stopifnot(requireNamespace("DefenseViz", quietly = TRUE))
suppressPackageStartupMessages(library(dplyr))
makeblastdb <- Sys.which("makeblastdb")
stopifnot(nzchar(makeblastdb))
reference <- DefenseViz::parse_rebase_sequences(raw)
required <- c("enzyme_name", "sequence", "enz_type", "rm_type", "subunit",
              "seq_length", "rec_seq", "is_gold_standard")
stopifnot(nrow(reference) > 0L, all(required %in% names(reference)))
dir.create(out, recursive = TRUE)
saveRDS(reference, file.path(out, "rebase_data.rds"))
fasta <- file.path(out, "rebase_db.fasta")
DefenseViz::write_fasta_for_blast(reference, fasta,
                                id_col = "enzyme_name", seq_col = "sequence")
status <- system2(makeblastdb, c("-in", shQuote(fasta), "-dbtype", "prot",
                               "-out", shQuote(fasta), "-title", "REBASE"))
stopifnot(status == 0L)
```

With the configured R dependency libraries available, execute:

```sh
Rscript --vanilla prepare_rebase.R /absolute/path/to/protein_seqs.txt \
  "$RESOURCES/rebase"
```

Set `R_LIBS` to the installed helper/dependency library paths if they are outside
R's default libraries. The historical R namespace name is preserved here because
it is the actual dependency API. These helpers retain the native FASTA filtering
and identifier encoding; do not substitute a generic FASTA writer. Record the
raw file hash, provider date/version, helper source pin, R/BLAST versions and
command output with the new set. A newly downloaded reference does not reproduce
an earlier unversioned snapshot. The setup script does not search genome queries
or change the analysis policy.

For the DNMBsuite image, mount the raw-file/script directory read-only at
`/source` and the new resource parent writable at `/resources`. Select the
actual R executable through `--entrypoint /usr/bin/env` with arguments
`-u R_HOME /opt/biotools/bin/Rscript --vanilla /source/prepare_rebase.R
/source/protein_seqs.txt /resources/rebase`, and set
`--env R_LIBS=/opt/defenseviz/native-r/library:/opt/biotools/lib/R/library`.
Clearing inherited `R_HOME` avoids mixing the base image's two R installations.

## Validate, then analyze

Set the engine's `reference_cache`, `hmm_database` and `profile_metadata` to these
prepared paths; keep its tool and R library settings from the runtime setup.
The [container engine](../docker/dnmbsuite/engine.json) already expects the layout
above at `/resources`. Validate the assembled resources using
`defenseviz check-resources --engine-config engine.json --out resource-check001`.
Preserve its report and input hashes. Then run
`defenseviz from-suite /path/to/suite-parent --engine-config engine.json
--evaluation full --out results/run001`, or use the complete
[Docker command](dnmbsuite-container.md#run-the-full-analysis-with-one-parent-folder).

Each check/analysis destination must be new. A setup or resource-check receipt
is not an analysis result; only the full command generates the evaluated CDS
lists, workbook and figures.
