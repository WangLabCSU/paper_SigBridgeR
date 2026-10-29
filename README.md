# Paper SigBridgeR 

Code for the paper 'SigBridgeR: An R Toolkit for Standardized Identification and Benchmarking of Phenotype-Associated Cell Subpopulations in Single-Cell Transcriptomics'

```r
setwd(usethis::proj_path())
# or
usethis::proj_activate(".")
```

## Project Structure

This repository is a sequential analysis pipeline, and the numeric prefixes of the
directories encode their order in that pipeline. Each stage reads the artifacts of the
previous one and writes its own results to the **external data directory**.

```
0_bench_data_preprocess   Raw data acquisition and harmonization
        |
        v
1_bench_screen            Screening with default parameters (exposes method quality)
        |
        v
2_method_acc              Hyper-parameter search + accuracy metrics
        |
        +--> 3_negative_ctrl / 3_neg_ctrl_ucell    Negative controls (random markers)
        +--> 4_positive_ctrl / 4_pos_ctrl_ucell    Positive controls (real DEG markers)
        +--> 5_GSEA                                Hallmark pathway enrichment
        +--> 9_method_time_cost                    Runtime / memory benchmark
        +--> paper_plot                            Manuscript figures
```

### 0_bench_data_preprocess — Data acquisition and harmonization

Downloads and standardizes the public cohorts used throughout the paper into a uniform
"expression matrix + phenotype label + survival information" form.

### 1_bench_screen — Screening with default parameters

Runs `SigBridgeR::Screen` with default settings for the eight benchmarked methods
(Scissor, scAB, SCIPAC, scPAS, scPP, DEGAS, LP-SGL, PIPET) and inspects whether the
resulting subpopulations are biologically plausible on UMAP.

### 2_method_acc — Method accuracy and hyper-parameter search

Benchmarks the accuracy of each method, running a two-round hyper-parameter search and
scoring every configuration against ground-truth labels.

### 3_negative_ctrl — Negative control (random markers, ssGSEA)

Tests whether the identified subpopulations can be reproduced by chance, by scoring each
cell population with 100 replicates of 20 randomly drawn markers.

### 3_neg_ctrl_ucell — Negative control, UCell implementation

The UCell counterpart of the negative control, run beside `3_negative_ctrl` rather than
after it.

### 4_positive_ctrl — Positive control (real DEG markers, ssGSEA)

Validates that SigBridgeR recovers subpopulations defined by genuine differentially
expressed markers.

### 4_pos_ctrl_ucell — Positive control, UCell implementation

UCell replication of the positive control, reusing the marker tables produced by
`4_positive_ctrl`.

### 5_GSEA — Hallmark pathway enrichment

Functional interpretation of the selected subpopulations through MSigDB Hallmark gene sets.

### 9_method_time_cost — Runtime and memory benchmark

Measures computational cost as a function of dataset size.

### paper_plot — Manuscript figures

Figures that are independent of the analysis results and used only in the manuscript,
together with the code that generates them.

- `Figure1A/` — conceptual panels: `GSEA/`, `boxplot/` (accuracy and F1),
  `bulk/` (bulk PCA), `gene_marker/` (marker volcano), `half_box_half_violine/`
  (negative-control mini panel), `pheno/` (phenotype heatmap), `random_gene_marker/`,
  `ssGSEA/` (score heatmap).
- `Figure5A/` — `R/` (`DEG_mock.R`, `benchmarkdata_sankey.R`, `cell_trajectory.R`) and
  `plot/` outputs (cohort-matching Sankey diagram, volcano plot, cell trajectory,
  DEG heatmap).

### deprecated — Archived approaches

Previous implementations kept for provenance: the early `8_ssGSEA/` marker-library
routes (BRCA markers, CellMarker2, Cell Taxonomy), the earlier per-dataset `GSEA/`
downstream and visualization scripts, and superseded `0_bench_data_preprocess/`,
`5_GSEA/`, `1_bench_screen/` variants.

## About

relevant tools: [SigBridgeR](https://github.com/WangLabCSU/SigBridgeR)

## License

- **Code** (`*.R`, `*.py`, `*.cpp`, etc): [MIT](LICENSE.md)
- **Data** (`*.csv`, `*.png`, `*.xlsx`, `*.qs2`, etc): [CC0 1.0](LICENSE-data.md)
