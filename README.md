# DefenseViz

[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-%3E%3D4.0-blue)]()

## Overview

**DefenseViz** is an R package for comprehensive analysis of bacterial Restriction-Modification (R-M) systems. It automatically detects R-M system components from genome annotation data and performs comparative analysis against the REBASE database to classify enzyme types, identify recognition sequences, and visualize system distributions.

## Key Features

- **Methyltransferase (MTase) Detection**: Identifies MTases using PFAM domain signatures and keyword-based annotation mining
- **Restriction Enzyme (REase) Detection**: Detects REases through PFAM domains, product annotations, and catalytic motif patterns
- **REBASE Integration**: Automated BLAST comparison against the REBASE Gold Standard database for enzyme classification
- **Operon Analysis**: Genomic context analysis to identify co-localized R-M system components
- **Catalytic Motif Detection**: Searches for conserved functional motifs (PD-(D/E)xK, DPPY, Walker A/B, DEAD/DEAH helicase motifs)
- **Recognition Sequence Extraction**: Retrieves DNA recognition sequences from REBASE matches
- **Visualization**: Generates publication-ready dotplots showing R-M system type and subunit distributions

## Installation

### From GitHub

```r
# Install devtools if not already installed
install.packages("devtools")

# Install DefenseViz from GitHub
devtools::install_github("your-username/DefenseViz")
```

### From Source

```r
# Clone the repository, then:
devtools::install("path/to/DefenseViz")

# Or using R CMD:
# R CMD build DefenseViz
# R CMD INSTALL DefenseViz_*.tar.gz
```

### Dependencies

#### Required R Packages

```r
# CRAN packages
install.packages(c(
  "dplyr",
  "tidyr",
  "stringr",
  "rlang",
  "readxl",
  "ggplot2",
  "httr",
  "jsonlite"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
```

#### Optional Packages

```r
# For Excel export (at least one recommended)
install.packages("xlsx")      # Full Excel support (requires Java)
install.packages("writexl")   # Lightweight alternative
```

#### External Tools

- **NCBI BLAST+**: Required for local BLAST searches against REBASE
  - Download: https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html
  - Ensure `blastp` is available in your system PATH

## Quick Start

```r
library(DefenseViz)

# Run the complete analysis pipeline
results <- rmscan_pipeline(
  input_file = "your_annotation.xlsx",
  output_dir = "."
)

# Print summary statistics
print_rmscan_summary(results)
```

## Pipeline Architecture

The DefenseViz pipeline consists of 12 sequential analysis steps:

```
┌─────────────────────────────────────────────────────────────────────────┐
│                    DefenseViz Analysis Pipeline                         │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  ┌─────────────────────────────────────────────────────────────────┐   │
│  │ STEP 1: Data Loading                                            │   │
│  │   • Load genome annotation file (Excel/CSV)                     │   │
│  │   • Detect available annotation sources (PFAM, InterPro, etc.)  │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│                              ↓                                          │
│  ┌─────────────────────────────────────────────────────────────────┐   │
│  │ STEP 2-3: Initial Candidate Detection                           │   │
│  │   • MTase detection (PFAM domains + keyword search)             │   │
│  │   • REase detection (PFAM domains + keyword search)             │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│                              ↓                                          │
│  ┌─────────────────────────────────────────────────────────────────┐   │
│  │ STEP 4: Operon Analysis                                         │   │
│  │   • Identify gene clusters based on genomic proximity           │   │
│  │   • Analyze strand orientation patterns                         │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│                              ↓                                          │
│  ┌─────────────────────────────────────────────────────────────────┐   │
│  │ STEP 5-6: REBASE Comparison                                     │   │
│  │   • BLAST search against REBASE Gold Standard database          │   │
│  │   • Filter results (identity ≥10%, alignment length ≥50)        │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│                              ↓                                          │
│  ┌─────────────────────────────────────────────────────────────────┐   │
│  │ STEP 7: Annotation & Motif Detection                            │   │
│  │   • Create annotated tables with BLAST results                  │   │
│  │   • Search for catalytic motifs in protein sequences            │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│                              ↓                                          │
│  ┌─────────────────────────────────────────────────────────────────┐   │
│  │ STEP 8-9: TRD & Recognition Sequence Extraction                 │   │
│  │   • Extract Target Recognition Domain (TRD) regions             │   │
│  │   • Retrieve recognition sequences from REBASE matches          │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│                              ↓                                          │
│  ┌─────────────────────────────────────────────────────────────────┐   │
│  │ STEP 10-12: Results Compilation & Export                        │   │
│  │   • Generate comprehensive summary table                        │   │
│  │   • Create R-M system distribution dotplot                      │   │
│  │   • Export results to Excel workbook                            │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

## Output Files

The pipeline generates two main output files:

| File | Description |
|------|-------------|
| `R-M_REBASE_analysis.xlsx` | Comprehensive analysis results (3 worksheets) |
| `RM_system_dotplot.pdf/png` | R-M system type and subunit distribution visualization |

### Excel Workbook Structure

| Worksheet | Contents |
|-----------|----------|
| **RM_Comprehensive** | All R-M candidate genes with complete annotation data |
| **High_Identity_50pct** | High-confidence matches (BLAST identity ≥50%) |
| **Type_Summary** | Count summary by R-M Type × Subunit |

### Output Column Descriptions

| Column | Description |
|--------|-------------|
| `locus_tag` | Gene identifier |
| `product` | Gene product annotation |
| `passed_blast_filter` | Boolean: passed BLAST quality filter |
| `rm_type` | R-M system type (Type I, II, III, or IV) |
| `subunit` | Functional subunit (M, R, or S) |
| `blast_match` | Best matching REBASE enzyme name |
| `rec_seq` | DNA recognition sequence |
| `blast_identity` | Sequence identity (0-1 scale) |
| `blast_evalue` | BLAST E-value |
| `operon_id` | Assigned operon cluster ID |
| `motif_confidence` | Catalytic motif detection confidence |
| `rease_motif_positions` | Detected REase motif positions |
| `mtase_motif_positions` | Detected MTase motif positions |

## R-M System Classification

### Type Classification

R-M system types are extracted from REBASE `enz_type` annotations:

| Type | Characteristics | Typical Components |
|------|-----------------|-------------------|
| **Type I** | Multi-subunit, ATP-dependent, cleave away from recognition site | HsdR (R), HsdM (M), HsdS (S) |
| **Type II** | Single or homodimeric, cleave within/near recognition site | REase (R), MTase (M) |
| **Type III** | Two-subunit, ATP-dependent, cleave ~25bp from recognition site | Mod (M), Res (R) |
| **Type IV** | Methylation-dependent restriction, no cognate MTase | McrA, McrBC, Mrr (R only) |

### Subunit Classification

Subunit types are determined from REBASE enzyme name prefixes:

| Prefix | Subunit | Function |
|--------|---------|----------|
| `M.XXX` | **M** (Methyltransferase) | DNA methylation |
| `R.XXX` | **R** (REase) | DNA restriction/cleavage |
| `S.XXX` | **S** (Specificity) | DNA recognition specificity (Type I only) |
| `XXX` (no prefix) | **R** | Default assignment for unprefixed enzymes |

### Fixed System Categories for Visualization

The dotplot displays all possible R-M system combinations:

- Type I: M, R, S subunits
- Type II: M, R subunits
- Type III: M, R subunits
- Type IV: R subunit only

## Catalytic Motif Detection

### REase Catalytic Motifs

| Motif | Regex Pattern | Biological Function |
|-------|---------------|---------------------|
| **PD-(D/E)xK** | `PD[DE].K` | Type II endonuclease catalytic core |
| **Walker A** | `G.{4}GK[TS]` | ATP-binding (P-loop) in helicases |
| **Walker B** | `[LIVMF]{4}DE` | ATP hydrolysis in helicases |
| **DEAD/DEAH** | `DEA[DH]` | Helicase superfamily 2 core motif |
| **SAT/SAH** | `S[ATV][THPSAQK]` | Helicase motif III |
| **TAN** | `T[AT][NT]` | Helicase motif III variant |
| **QxxR** | `Q..R` | Helicase motif VI |
| **ARGID** | `[AR]G[IL]D` | Helicase motif V |

### MTase Catalytic Motifs

| Motif | Regex Pattern | Biological Function |
|-------|---------------|---------------------|
| **DPPY** | `[DNSH]PP[YFW]` | Catalytic motif IV (base flipping) |
| **NPPY** | `NPP[YF]` | Motif IV variant |
| **FxGxG** | `[FY].G.G` | SAM-binding motif I |

## Advanced Usage

### Custom Pipeline Parameters

```r
results <- rmscan_pipeline(
  input_file = "annotation.xlsx",
  output_dir = "results",


  # BLAST Filtering Thresholds
  blast_min_identity = 0.30,    # Minimum identity (default: 0.10 = 10%)
  blast_min_length = 100,       # Minimum alignment length (default: 50)

  # Operon Detection Parameters
  max_operon_gap = 3000,        # Maximum intergenic distance (default: 5000 bp)
  max_intervening = 2,          # Maximum intervening genes (default: 1)

  # Pipeline Options
  compare_rebase = TRUE,        # Perform REBASE comparison (default: TRUE)
  search_motifs = TRUE,         # Search for TRD motifs (default: TRUE)
  verbose = TRUE,               # Print progress messages (default: TRUE)
  save_intermediates = TRUE     # Save intermediate results to global env (default: TRUE)
)
```

### Accessing Intermediate Results

When `save_intermediates = TRUE`, intermediate data objects are saved to the global environment:

```r
# Available after pipeline execution:
rmscan_dnmb                  # Original input data
rmscan_mtase                 # MTase candidates
rmscan_rease                 # REase candidates
rmscan_operon                # Operon analysis results
rmscan_blast_filtered        # Filtered BLAST results
rmscan_rebase_db             # Cached REBASE database
rmscan_rm_comprehensive      # Final comprehensive table
rmscan_catalytic_summary     # Catalytic motif summary
```

### Using Individual Functions

```r
library(DefenseViz)

# Step 1: Load annotation data
dnmb_data <- load_dnmb("annotation.xlsx")

# Step 2: Detect MTases
mtase_candidates <- detect_methyltransferases(dnmb_data)

# Step 3: Detect REases
rease_candidates <- detect_restriction_enzymes(dnmb_data)

# Step 4: Identify operons
operons <- identify_rm_operons(mtase_candidates, rease_candidates, dnmb_data)

# Step 5: Get REBASE data
rebase_db <- get_rebase_data()

# Step 6: Detect catalytic motifs
motifs <- detect_catalytic_motifs(candidates, seq_col = "translation")
motif_summary <- summarize_catalytic_motifs(motifs)

# Step 7: Generate dotplot
generate_rm_type_heatmap(rm_data, output_dir = ".", min_identity = 0.5)
```

## Input Data Format

### Required File Format

The input file should be an Excel (.xlsx) or CSV file containing genome annotation data.

### Required Columns

| Column | Description | Example |
|--------|-------------|---------|
| `locus_tag` or `protein_id` | Unique gene identifier | `ECOLI_00123` |
| `start` | Start coordinate | `12345` |
| `end` | End coordinate | `13456` |
| `strand` or `direction` | Coding strand | `+` or `-` |
| `product` | Gene product description | `DNA methyltransferase` |
| `translation` | Protein sequence | `MKFLILLFNILC...` |

### Optional Columns

| Column | Description |
|--------|-------------|
| `PFAMs` or `Signature accession_Pfam` | PFAM domain annotations |
| `seqid` or `contig` | Chromosome/contig identifier |
| `InterPro` or `Signature accession_CDD` | Additional domain annotations |

## System Requirements

| Requirement | Specification |
|-------------|---------------|
| **R Version** | ≥ 4.0 |
| **Memory** | 4 GB RAM minimum (8 GB recommended for large genomes) |
| **Disk Space** | ~100 MB for REBASE cache |
| **BLAST+** | Required for REBASE comparison (optional if skipping BLAST) |

## Troubleshooting

### Common Issues

**BLAST not found**
```r
# Check if blastp is in PATH
system("which blastp")

# Or specify full path in environment
Sys.setenv(PATH = paste("/path/to/blast/bin", Sys.getenv("PATH"), sep = ":"))
```

**Memory issues with large genomes**
```r
# Process in chunks or increase memory limit
options(future.globals.maxSize = 2000 * 1024^2)  # 2GB
```

**Excel export fails**
```r
# Try alternative package
install.packages("writexl")  # Doesn't require Java
```

## Citation

If you use DefenseViz in your research, please cite:

```
DefenseViz: An R Package for Automated Analysis of Bacterial
Restriction-Modification Systems
https://github.com/your-username/DefenseViz
```

## References

1. Roberts RJ, Vincze T, Posfai J, Macelis D. (2023) REBASE - a database for DNA restriction and modification: enzymes, genes and genomes. *Nucleic Acids Research*, 51(D1):D629-D635. doi: 10.1093/nar/gkac975

2. Loenen WAM, Dryden DTF, Raleigh EA, Wilson GG. (2014) Type I restriction enzymes and their relatives. *Nucleic Acids Research*, 42(1):20-44.

3. Pingoud A, Wilson GG, Wende W. (2014) Type II restriction endonucleases - a historical perspective and more. *Nucleic Acids Research*, 42(12):7489-7527.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

**Jae-Yoon Sung**
Email: sungjaeyoon92@gmail.com

---

*DefenseViz - Comprehensive R-M System Analysis for Bacterial Genomics*
