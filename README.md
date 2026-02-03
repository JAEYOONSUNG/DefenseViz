# DefenseViz

[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Restriction-Modification (R-M) System Analysis Pipeline

**DefenseViz**는 세균 게놈의 DNMB 어노테이션 데이터에서 R-M 시스템을 자동으로 탐지하고 REBASE 데이터베이스와 비교 분석하는 R 패키지입니다.

## Features

- **MTase Detection**: PFAM 도메인 및 키워드 기반 Methyltransferase 탐지
- **REase Detection**: PFAM 도메인, 키워드, catalytic motif 기반 Restriction Enzyme 탐지
- **REBASE BLAST**: REBASE 데이터베이스와 BLAST 비교로 R-M 타입 및 인식서열 확인
- **Operon Analysis**: R-M 시스템 오페론 구조 분석
- **Catalytic Motif Detection**: PDxDxK, DPPY, Walker A/B, DEAD/DEAH 등 핵심 모티프 탐지
- **Visualization**: R-M 시스템 Type/Subunit 분포 dotplot

## Installation

```r
# GitHub에서 설치
devtools::install_github("your-username/DefenseViz")

# 또는 로컬에서 빌드 및 설치
devtools::install("path/to/DefenseViz")

# 또는
R CMD build DefenseViz
R CMD INSTALL DefenseViz_0.1.0.tar.gz
```

### Dependencies

**Required:**
```r
install.packages(c("dplyr", "tidyr", "stringr", "rlang", "readxl", "ggplot2", "httr", "jsonlite"))

# Bioconductor
BiocManager::install("Biostrings")
```

**Optional:**
```r
install.packages(c("xlsx", "writexl"))  # Excel 출력
```

**External:**
- BLAST+ (NCBI): 로컬 BLAST 실행에 필요

## Quick Start

```r
library(DefenseViz)

# 전체 파이프라인 실행
results <- rmscan_pipeline(
  input_file = "your_dnmb_annotation.xlsx",
  output_dir = "rmscan_results"
)

# 결과 요약 출력
print_rmscan_summary(results)
```

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────────┐
│                      DefenseViz Pipeline (12 Steps)                 │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  [1] Load DNMB Data ─────────────────────────────────────────────  │
│       ↓                                                             │
│  [2] Detect MTase (PFAM + Keyword) ──────────────────────────────  │
│       ↓                                                             │
│  [3] Detect REase (PFAM + Keyword) ──────────────────────────────  │
│       ↓                                                             │
│  [4] Operon Analysis (Genome Context) ───────────────────────────  │
│       ↓                                                             │
│  [5] REBASE BLAST Comparison ────────────────────────────────────  │
│       ↓                                                             │
│  [6] Filter BLAST Results (identity ≥10%, length ≥50) ───────────  │
│       ↓                                                             │
│  [7] Create Annotated Tables (BLAST + Operon info) ──────────────  │
│       ↓                                                             │
│  [7b] Detect Catalytic Motifs (PDxDxK, DPPY, Walker, DEAD...) ───  │
│       ↓                                                             │
│  [8] TRD Extraction ─────────────────────────────────────────────  │
│       ↓                                                             │
│  [9] Extract Recognition Sequences from REBASE ──────────────────  │
│       ↓                                                             │
│  [10] Create Comprehensive Table ────────────────────────────────  │
│       ↓                                                             │
│  [11] Generate R-M Type Dotplot (identity ≥50%) ─────────────────  │
│       ↓                                                             │
│  [12] Export Results (xlsx) ─────────────────────────────────────  │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

## Output Files

| File | Description |
|------|-------------|
| `R-M_REBASE_analysis.xlsx` | 종합 분석 결과 Excel (3 sheets) |
| `RM_system_dotplot.pdf` | R-M 시스템 Type/Subunit 분포 시각화 |
| `RM_system_dotplot.png` | Dotplot (PNG 형식) |

### Excel Sheets

| Sheet | Description |
|-------|-------------|
| **RM_Comprehensive** | 전체 R-M 후보 유전자 (모든 데이터) |
| **High_Identity_50pct** | BLAST identity ≥50% 고신뢰 결과 |
| **Type_Summary** | R-M Type × Subunit 카운트 요약 |

## Output Columns

| Column | Description |
|--------|-------------|
| `locus_tag` | 유전자 ID |
| `product` | 유전자 기능 설명 |
| `passed_blast_filter` | BLAST 필터 통과 여부 |
| `rm_type` | R-M 시스템 타입 (Type I ~ IV) |
| `subunit` | 서브유닛 (M, R, S) |
| `blast_match` | REBASE 매칭 효소명 |
| `rec_seq` | 인식 서열 (Recognition Sequence) |
| `blast_identity` | BLAST identity (0~1) |
| `operon_id` | 오페론 ID |
| `motif_confidence` | 모티프 신뢰도 |
| `rease_motif_positions` | REase 모티프 위치 |
| `mtase_motif_positions` | MTase 모티프 위치 |

## R-M System Classification

### Type (from REBASE enz_type)

| Type | Description |
|------|-------------|
| Type I | Type I R-M system (HsdR, HsdM, HsdS) |
| Type II | Type II (including IIG, IIS subtypes) |
| Type III | Type III (Mod, Res) |
| Type IV | Methylation-dependent restriction |

### Subunit (from enzyme_name prefix)

| Prefix | Subunit | Function |
|--------|---------|----------|
| `M.XXX` | **M** | Methyltransferase |
| `R.XXX` | **R** | Restriction endonuclease |
| `S.XXX` | **S** | Specificity subunit |
| `XXX` (no prefix) | **R** | Default (REase) |

## Catalytic Motifs Detected

### REase Motifs

| Motif | Pattern | Function |
|-------|---------|----------|
| PD-(D/E)xK | `PD[DE].K` | Type II nuclease catalytic |
| Walker A | `G.{4}GK[TS]` | ATP binding (helicase) |
| Walker B | `[LIVMF]{4}DE` | ATP hydrolysis |
| DEAD/DEAH | `DEA[DH]` | Helicase core |
| SAT/SAH | `S[ATV][THPSAQK]` | Helicase motif III |
| TAN | `T[AT][NT]` | Helicase motif III variant |
| QxxR | `Q..R` | Helicase motif VI |
| ARGID | `[AR]G[IL]D` | Helicase motif V |

### MTase Motifs

| Motif | Pattern | Function |
|-------|---------|----------|
| DPPY | `[DNSH]PP[YFW]` | Catalytic motif IV |
| NPPY | `NPP[YF]` | Motif IV variant |
| FxGxG | `[FY].G.G` | SAM-binding motif I |

## Advanced Usage

### Custom Parameters

```r
results <- rmscan_pipeline(
  input_file = "annotation.xlsx",
  output_dir = "output",

  # BLAST filtering
  blast_min_identity = 0.30,   # Default: 0.10 (10%)
  blast_min_length = 100,      # Default: 50

  # Operon detection
  max_operon_gap = 3000,       # Default: 5000 bp
  max_intervening = 2,         # Default: 1 gene

  # Options
  compare_rebase = TRUE,       # Run REBASE comparison
  search_motifs = TRUE,        # Search TRD motifs
  verbose = TRUE,              # Print progress
  save_intermediates = TRUE    # Save to global environment
)
```

### Access Intermediate Results

```r
# 파이프라인 실행 후 (save_intermediates = TRUE 필요)
rmscan_rm_comprehensive   # 종합 테이블
rmscan_mtase_annotated    # MTase 어노테이션
rmscan_rease_annotated    # REase 어노테이션
rmscan_blast_filtered     # 필터링된 BLAST 결과
rmscan_rebase_db          # REBASE 데이터베이스 캐시
rmscan_operon             # 오페론 분석 결과
rmscan_catalytic_summary  # 모티프 요약
```

### Individual Functions

```r
# 1. 데이터 로드
dnmb_data <- load_dnmb("annotation.xlsx")

# 2. MTase 탐지
mtase <- detect_methyltransferases(dnmb_data)

# 3. REase 탐지
rease <- detect_restriction_enzymes(dnmb_data)

# 4. REBASE 데이터 가져오기
rebase_db <- get_rebase_data()

# 5. Catalytic motif 탐지
motifs <- detect_catalytic_motifs(candidates, seq_col = "translation")
motif_summary <- summarize_catalytic_motifs(motifs)

# 6. 오페론 분석
operons <- identify_rm_operons(mtase, rease, dnmb_data)
```

## Input Format

DNMB 어노테이션 Excel/CSV 파일:

**Required columns:**
- `locus_tag` 또는 `protein_id`: 유전자 식별자
- `start`, `end`: 좌표
- `strand` 또는 `direction`: 방향 (+/-)
- `product`: 기능 설명
- `translation`: 단백질 서열

**Optional columns:**
- `PFAMs` 또는 `Signature accession_Pfam`: PFAM 도메인
- `seqid` 또는 `contig`: 염색체/콘티그 ID

## System Requirements

- **R**: ≥ 4.0
- **Memory**: 4GB+ 권장
- **Disk**: ~100MB (REBASE 캐시)
- **BLAST+**: NCBI BLAST+ (optional, 로컬 BLAST용)

## Citation

```
DefenseViz: R-M System Analysis Pipeline for Bacterial Genomes
https://github.com/your-username/DefenseViz
```

## References

- Roberts RJ, et al. (2023) REBASE - a database for DNA restriction and modification.
  Nucleic Acids Res. 51(D1):D629-D635.

## License

MIT License

## Author

- **Jae-Yoon Sung**
- Email: o3wodbs@gmail.com
