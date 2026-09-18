# Transcriptomic Data Analysis & Genomic Pipelines

## Overview
This repository contains a curated collection of custom R scripts, automation pipelines, and interactive web tools developed to process, filter, cross-validate, and interpret large-scale transcriptomic datasets (bulk RNA-seq). These workflows bridge raw statistical matrices with downstream biological systems interpretation, utilizing strict genomic control metrics to eliminate statistical bias.

## Repository Directory & Toolsets

### 1. Interactive Web Applications (R Shiny Framework)
*   **Directory:** `/shiny-venn-diagrams/app.R`  
    *   **Purpose:** Multi-contrast differential gene expression (DEG) alignment GUI. Automatically extracts up/down-regulated genes from divergent contrasts, interactively scales and adjusts rendering layouts, and exports intersecting coordinate sheets alongside vector files.
*   **Directory:** `/shiny-volcano-plots/app.R`  
    *   **Purpose:** Automated Volcano Plot configuration suite. Supports automated parsing across multi-format outputs (`.csv`, `.tsv`, `.txt`, `.xlsx`), dynamic expression pattern labeling via `ggrepel`, contrast flipping logic, and exports 300-DPI publication graphics.

### 2. Functional Enrichment Engines
*   **Filename:** `go_ora_pipeline.R`  
    *   **Purpose:** Hypergeometric enrichment profiling across BP, MF, and CC domains. Implements a stringently filtered gene background universe constraint to eliminate false-positive pathway inflation.
*   **Filename:** `kegg_enrichment_pipeline.R`  
    *   **Purpose:** Programmatic annotation mapping over Kyoto Encyclopedia of Genes and Genomes databases tailored specifically for murine models (`organism = "mmu"`).
*   **Filename:** `msigdb_hallmark_pipeline.R`  
    *   **Purpose:** Incorporates the Molecular Signatures Database (MSigDB) Hallmark collections via generalized computational `enricher` pipelines.

### 3. Matrix Manipulation Utilities
*   **Filename:** `log2fc_contrast_flipper.R`  
    *   **Purpose:** Inverts log2FoldChange contrast configurations via cross-platform graphic UI (`tcltk`) paths to switch control-versus-treatment variables without re-processing.

---

## Technical Environment & Dependencies
*   `shiny` — Reactive cross-platform user interface pipelines.
*   `ggplot2` & `ggrepel` — Spatial geometry rendering and collision resolution.
*   `dplyr`, `readxl`, `readr`, `tcltk` — Data frame manipulation and data input streams.
*   `clusterProfiler`, `msigdbr`, `org.Mm.eg.db` — Annotation mapping and database queries.
*   `openxlsx`, `grid` — Multilayer cell formatting and workbook layout compilation.

## How to Run
1. Clone this repository locally.
2. Verify package attachments via `install.packages()`.
3. Launch web apps by invoking `shiny::runApp("directory_name")` in your console.

## License
This repository is licensed under the permissive MIT License. See the `LICENSE` file for full text details.
