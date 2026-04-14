# Gliosis in Degeneration Models


This repository contains a memory-optimized computational pipeline for analyzing RNA sequencing data from retinal degeneration mouse models. 

## Biological Goal
To identify differentially expressed genes associated with gliosis (specifically the reactivity of astrocytes, Müller glia, and microglia) in retinal degeneration models (such as rd1 and rd10) when compared against wild-type controls.

## Research Question
*Is it feasible to identify viable candidate targets for immunohistochemistry (IHC) using public RNA-seq data processed on resource-limited hardware (e.g., 16 GB of RAM)?*

## Analytical Methodology
1. **Data Acquisition**: Retrieval of single-cell datasets from public repositories.
2. **Quality Control (QC)**: Robust filtering of ambient RNA contamination and low-quality cells.
3. **Pseudo-Bulk Aggregation**: Aggregation of single cells into pseudo-bulk expression matrices to prepare for robust statistical testing.
4. **Differential Expression Analysis**: Utilizing `DESeq2` with applied shrinkage estimators (e.g., apeglm) to compute reliable log2 fold changes.
5. **Functional Enrichment Analysis**: Assessing the overrepresentation of specific biological pathways (such as gliogenesis) using tools like `clusterProfiler`.

## Infrastructure Requirements & Reproducibility
The computational workflow is completely self-contained, auditable, and reproducible. 
It requires the use of isolated containers (via Docker). To guarantee stability across local environments and circumvent the resource-intensive compilation of R packages from source, the environment installs dependencies using native operating system binaries.

### Required R Libraries

- `GEOquery`
- `Matrix`
- `Seurat`
- `DESeq2`
- `clusterProfiler`
- `org.Mm.eg.db`
- `EnhancedVolcano`
- `patchwork`

## Repository Layout

* **`Dockerfile`**: Defines the reproducible container environment. It handles the installation of all necessary system dependencies and R packages via pre-compiled binaries.
* **`analisis.R`**: The core, fully-automated R script encompassing the full analysis pipeline. It employs strict memory management (`rm()`, `gc()`) to ensure successful execution within the 16 GB RAM constraint.
* **`Report.Rmd` / `Report.pdf`**: The comprehensive analytical report detailing the methodology, rationale, and interpretation of the results. 
* **`Rplots.pdf`**: A collection of generated visual outputs, including UMAP projections and Volcano plots highlighting differential expression.
* **`resultados/`**: The output directory containing all exported datasets, result matrices, and significant DEGs identified by the pipeline.
* **`GSM*_data/`**: Directories containing the raw input matrices from public repositories:
  * `GSM3854512_data/`
  * `GSM8728964_data/`
  * `GSM8728965_data/`
  * `GSM8728966_data/`
  * `GSM8728967_data/`
  * `GSM8728968_data/`

*(Note: "Other Drafts" directory contains deprecated intermediate files and is deliberately excluded from this workflow).*

## Data Sources & Bibliography

* **Simon CJ, Khabou H, Chaffiol A, Rucli M et al.** *Reactivating the phototransduction cascade with a mutation agnostic gene therapy preserves vision in rod-cone dystrophies.* iScience 2025. Data from rd10 and rd1 models are obtained from this study.
* **Heng JS, Hackett SF, Stein-O'Brien GL, Winer BL et al.** *Comprehensive analysis of a mouse model of spontaneous uveoretinitis using single-cell RNA sequencing.* Proc Natl Acad Sci U S A 2019. Data from the wild-type (C57) model are obtained from this study.
