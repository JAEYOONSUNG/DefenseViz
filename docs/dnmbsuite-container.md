# DNMBsuite container with DefenseViz

`Dockerfile.dnmbsuite` adds the complete DefenseViz analysis and export stage to a pinned DNMBsuite image. It includes Python, the native R helper packages, BLAST, HMMER and CATH. Databases, inputs and optional licensed fonts are supplied separately. The original DNMBsuite entrypoint and tools remain available in the same image.

This recipe targets **Linux amd64**. Use `--platform linux/amd64` on Apple Silicon; execution then uses emulation. A native installation is also supported; Docker is optional. The smaller default `Dockerfile` is a Python-only package image and does not provide this full native engine.

## Build from GitHub

From a clean checkout at the revision you intend to use:

```sh
git clone https://github.com/JAEYOONSUNG/DefenseViz.git
cd DefenseViz
# Optionally: git checkout <published-tag-or-commit>
docker build --platform linux/amd64 \
  --build-arg DEFENSEVIZ_REVISION="$(git rev-parse HEAD)" \
  -f Dockerfile.dnmbsuite -t dnmbsuite-defenseviz:local .
```

The Dockerfile-specific ignore file sends only package code, its installer and container files. It excludes databases, genomes, output folders and local fonts. The base is pinned by digest; CATH is pinned by version and SHA256; the R helper installer verifies exact public source commits. Resolved Python packages and R source/install receipts remain under `/opt/defenseviz`. Record the built image ID with `docker image inspect dnmbsuite-defenseviz:local --format '{{.Id}}'` when sharing results. Python dependencies are resolved during the build, so preserving the final image digest is necessary to reproduce that exact environment.

The base uses `/opt/biotools/bin/Rscript` (measured R4.3.1), not the R version advertised by the Rocker base environment. The dedicated DefenseViz executable clears the inherited conflicting `R_HOME` for its own process. It does not replace the suite's global PATH, Python or R library.

## Prepare the resources

Mount the database directory with this layout:

```text
validated-databases/
  rebase/
    rebase_data.rds
    rebase_db.fasta
    rebase_db.fasta.phr
    rebase_db.fasta.pin
    rebase_db.fasta.psq
    ...other BLAST index files...
  rm-hmm/
    RM_HMMs.hmm
    profile_roles.tsv
```

Use the validated resource export described in [runtime-bundle.md](runtime-bundle.md), or an equivalent source-traceable DNMB resource set. The paired HMM and role metadata must describe the same profiles. Copying a Mac runtime bundle wholesale into Linux is not supported; mount only its `databases` directory. Database distribution terms and provenance remain with the resources.

Set absolute host paths; create a fresh writable output parent:

```sh
RESOURCES=/absolute/path/to/validated-databases
SUITE_RESULTS=/absolute/path/to/completed-suite-parent
RESULTS=/absolute/path/to/new-defenseviz-results
mkdir -p "$RESULTS"
```

Validate resources once inside the image:

```sh
docker run --rm --platform linux/amd64 \
  --mount "type=bind,src=$RESOURCES,dst=/resources,readonly" \
  --mount "type=bind,src=$RESULTS,dst=/output" \
  --entrypoint /opt/defenseviz/bin/defenseviz \
  dnmbsuite-defenseviz:local check-resources \
  --engine-config /opt/defenseviz/container/engine.json \
  --out /output/resource-check
```

This verifies the actual R/tool environment, profile metadata and reference index/FASTA binding. It does not run the biological analysis or certify biological activity.

## Run the full analysis with one parent folder

```sh
docker run --rm --platform linux/amd64 \
  --mount "type=bind,src=$RESOURCES,dst=/resources,readonly" \
  --mount "type=bind,src=$SUITE_RESULTS,dst=/input,readonly" \
  --mount "type=bind,src=$RESULTS,dst=/output" \
  --entrypoint /opt/defenseviz/bin/defenseviz \
  dnmbsuite-defenseviz:local from-suite /input \
  --evaluation full --out /output/run001
```

The parent may contain nested per-genome DNMBsuite results. The CLI discovers completed native roots, performs fresh full-profile HMM and nominated-candidate reference searches, evaluates R-M/system/effector evidence, and writes classified Excel, filtered annotation lists, comparison tables, publication figures and provenance. No manual cohort or evidence manifest is needed. Original DF/PAD/DP observations are imported from the suite outputs; this stage does not claim to rerun those upstream programs.

Inputs and databases are read-only. Temporary GenBank/FASTA/search files are isolated in the fresh result directory. Choose a new `run001` name for another run. The bundled engine uses four CPUs. For a different CPU setting or resource path, copy `docker/dnmbsuite/engine.json`, edit the relevant fields, mount that file read-only at `/engine.json`, and add `--engine-config /engine.json`.

On Linux, `--user "$(id -u):$(id -g)"` may be added so output files belong to the invoking user; the output directory must be writable by that UID. The app wrapper uses writable temporary font/cache directories. Do not add that option to the original DNMBsuite workflow without checking its own UID/cache behavior.

## Arial and figure appearance

Arial is not redistributed in this image. To retain Arial, prepare a local directory containing licensed `Arial.ttf`, `Arial Bold.ttf`, `Arial Italic.ttf` and `Arial Bold Italic.ttf`, then add this mount to the DefenseViz command:

```sh
--mount "type=bind,src=/absolute/local/arial-fonts,dst=/usr/local/share/fonts/defenseviz,readonly"
```

The wrapper refreshes Fontconfig in an app-local cache. On Docker Desktop, protected system font directories may not be shared; a private byte-identical local copy in a shared folder works. Do not commit that directory or COPY it into a public image. Without Arial, the renderer records its available fallback (for this base, DejaVu Sans); inspect figure typography and layout rather than assuming identical pixels across fonts or platforms.

## Original DNMBsuite commands in the same image

The inherited entrypoint, default command and `/data` working directory are unchanged. For example, with a writable working copy of GenBank inputs and a separate writable suite cache:

```sh
docker run --rm --platform linux/amd64 \
  --mount type=bind,src=/absolute/suite-work,dst=/data \
  --mount type=bind,src=/absolute/suite-cache,dst=/opt/dnmb-cache \
  --env DNMB_AUTO_UPDATE=0 --env DNMB_INTERPROSCAN_UPDATE=never \
  dnmbsuite-defenseviz:local run /data
```

After the suite finishes, use the parent-folder command above on its completed outputs. The frozen analysis resources and writable upstream suite cache are different directories; do not bind the read-only reference resources as the suite cache. See [DNMBsuite](https://github.com/JAEYOONSUNG/DNMBsuite) for upstream module options.

Historical R package names and third-party licenses are documented in [the container notices](../docker/dnmbsuite/THIRD_PARTY.md). The scientific pipeline remains DNMB; the Python DefenseViz entrypoint runs its evidence review and publication export stage.
