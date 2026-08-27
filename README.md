# Coudray2026

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![R](https://img.shields.io/badge/R-%3E%3D%204.0-276DC3.svg)
[![ClimbTheCliff](https://img.shields.io/badge/package-ClimbTheCliff-blue.svg)](https://github.com/alexdray86/ClimbTheCliff)
[![Preprint](https://img.shields.io/badge/bioRxiv-10.1101%2F2023.05.10.540140-b31b1b.svg)](https://www.biorxiv.org/content/10.1101/2023.05.10.540140v3)

Reproducibility repository for `Computational deconvolution of drug sensitivity via single-cell–to–bulk transcriptome mapping`, Coudray et al. — [preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2023.05.10.540140v3) (posted under an earlier title, will be updated on publication). Contains the R notebooks that generate every main and supplemental figure. The underlying CLIMB and CLIFF deconvolution methods live in the companion R package [ClimbTheCliff](https://github.com/alexdray86/ClimbTheCliff), including full method documentation and a small standalone example dataset — see [alexdray86.github.io/ClimbTheCliff](https://alexdray86.github.io/ClimbTheCliff/).

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Data](#data)
- [Notebook Descriptions](#notebook-descriptions)
- [Shared Helper Functions](#shared-helper-functions)
- [Python Environment for TAPE and Scaden](#python-environment-for-tape-and-scaden)
- [Citation](#citation)
- [License](#license)
- [Issues and Support](#issues-and-support)

---

## Requirements

- R > 4.0, run through the [IRkernel](https://irkernel.github.io/) for Jupyter
- [ClimbTheCliff](https://github.com/alexdray86/ClimbTheCliff), pinned to the frozen release used throughout this repository (`1.1.0`) — see [Installation](#installation)
- ClimbTheCliff's own dependencies, attached alongside it: [`Biobase`](https://bioconductor.org/packages/Biobase/), [`glmnet`](https://cran.r-project.org/package=glmnet), [`reshape2`](https://cran.r-project.org/package=reshape2)
- Additional packages used by individual notebooks for comparison methods and preprocessing (e.g. `Seurat`, `MuSiC`, `BisqueRNA`, `matrixStats`) — each notebook lists what it needs in its first cell; install those before running it
- A Python environment for the TAPE and Scaden comparison methods — see [Python Environment for TAPE and Scaden](#python-environment-for-tape-and-scaden)

## Installation

Clone this repository, then install the exact frozen `ClimbTheCliff` release used to produce these results:

```r
# Using devtools
devtools::install_github("alexdray86/ClimbTheCliff", ref = "1.1.0", force = T)

# Or using remotes
remotes::install_github("alexdray86/ClimbTheCliff", ref = "1.1.0", force = T)
```

```bash
git clone https://github.com/alexdray86/Coudray2026.git
cd Coudray2026
```

Then download the data bundle — see [Data](#data) below — before running any notebook.

---

## Data

All data used across the notebooks are stored in a single, publicly available Google Drive bundle:

**Download link:** https://drive.google.com/file/d/1mL-ZxlkETq9EzZ85TP6Jmtx54NNv1cx1/view?usp=sharing

Unzip `data.zip` and place the resulting `data/` folder at the root of this repository.

| Folder | Contents |
|---|---|
| `invitro_experiment/` | Raw single-cell and bulk RNA-seq (`ExpressionSet`) objects for the 4-cell-line in-vitro mixtures, plus their measured drug response (Fig. 2, Fig. 5) |
| `pseudobulks_climb/` | Paired single-cell references and pseudo-bulk mixtures for each cross-dataset pair (AML, CRC, MEL, BREAST, GBM) used in the deconvolution benchmark (Fig. 3) |
| `aml_3cohorts/` | BeatAML, TCGA-LAML and Leucegene bulk RNA-seq, LSC17 scores, and precomputed deconvolution results for the three-cohort AML analysis (Fig. 3.4) |
| `beataml_drug/` | BeatAML drug screening (AUC) data and precomputed CLIMB expression output, for cell-type-specific drug sensitivity inference (Fig. 6) |
| `simulation_aml/` | Simulated pseudo-bulks and drug sensitivity ground truth used to benchmark CLIFF against alternative methods (Fig. 6 simulation) |
| `cliff_sc/` | Van Galen single-cell drug screening data and single-cell-level CLIFF-SC inputs/outputs (Fig. 7) |
| `GeneLength.txt` | Gene length reference required by TAPE/Scaden |
| `review_litterature_drug_sensitivity.csv` | Literature-curated drug sensitivity values used for cross-checks |

Some notebooks also write intermediate results back into `data/` so that later steps (or later re-runs) don't have to repeat expensive deconvolution methods.

---

## Notebook Descriptions

Each notebook corresponds to a specific figure in the manuscript. Below is a summary of each notebook's contents.

### Fig2_invitro_experiment.ipynb

In-vitro experiments with cell line mixtures created through precise cell counting. Cell line proportions were documented and mixtures were subjected to bulk RNA-seq. This notebook applies various deconvolution techniques to infer cell line fractions using a single-cell RNA-seq reference dataset of the same four cell lines.

### Fig3_1_preproc_integration.ipynb

Integration strategy for pairs of single-cell datasets:
- Extract cell and gene subsets
- Apply Seurat integration pipeline to merge datasets
- Transfer cell subtype labels between datasets
- Create pseudo-bulks with known cell subtype proportions
- Store datasets as RDS objects for downstream deconvolution

### Fig3_2_deconvolution_panel.ipynb

Applies deconvolution techniques in a cross-dataset framework using processed datasets from Fig3_1. Results are saved for analysis in Fig3_3.

### Fig3_3_analysis.ipynb

Compares deconvolution outcomes against ground truth cell subtype proportions using four metrics. Generates boxplots and ranking plots for all twelve analyses.

### Fig3_4_amlCohorts.ipynb

Executes CLIMB deconvolution across three AML patient cohorts and visualizes results as heatmaps. Compares CLIMB results with Van Galen signature scores and LSC17 signature scores.

### Fig3_5_CLIMB_simulation.ipynb

Simulation framework to validate CLIMB's performance:
- Creates pseudo-bulk samples from single-cell reference with known cell-type proportions
- Samples cells according to specified proportions to generate realistic bulk mixtures
- Compares CLIMB deconvolution against other methods using simulated ground truth
- Evaluates correlation between pseudo-bulks and signature matrices

### Fig5_CLIFF_invitro.ipynb

CLIFF deconvolution of drug sensitivity data from in-vitro experiments. Ground truth based on drug sensitivity screening of individual pure cell lines.

### Fig6_CLIFF_beataml.ipynb

CLIFF deconvolution on BeatAML drug screening data using CLIMB-deconvoluted cell subtype proportions and expression as reference. Generates figures visualizing the results.

### Fig6_CLIFF_simulation.ipynb

Simulation of drug sensitivity screening data at BeatAML scale:
- Uses 200 simulated pseudo-bulks to simulate drug screening data with diverse patterns
- Performs CLIFF deconvolution on simulated datasets
- Compares CLIFF results with alternative methods (Bottomly-like approach, Zeng, and linear regression using CLIMB proportions)

### Fig7_CLIFF-SC.ipynb

CLIFF-SC: Single-cell level drug sensitivity prediction (adaptation of CLIFF). Applies to BeatAML dataset to deconvolute Van Galen's single-cell drug screening data using single-cell to bulk coefficients from CLIMB-SC.

### SFig10_PLOSone_CompBio.ipynb

Round-2 revision notebook for the PLOS Computational Biology submission. Builds the 7-panel supplemental figure (Supplemental Fig. 10) addressing all reviewer comments: CLIMB component ablation (A), reference coverage sweep (B), CLIFF identifiability (C), EM initialization sensitivity (D-E), runtime scaling (F), and uncertainty propagation from CLIMB into CLIFF (G).

---

## Shared Helper Functions

[`shared_functions.R`](shared_functions.R) collects a handful of small utility functions (metric computation, cell-type label reformatting, TAPE/Scaden result loading, etc.) that are used identically across most notebooks. Each notebook `source()`s it instead of redefining these functions inline. Functions that genuinely differ between notebooks (e.g. `Fig3_2_deconvolution_panel.ipynb`'s own `reformat_celltypes`) are kept local to that notebook.

---

## Python Environment for TAPE and Scaden

To run TAPE and Scaden deconvolution methods, you need to set up a Python environment with the dependencies listed in `requirements_tape_scaden.txt`.

### Setup Instructions

1. Create a new Python environment:
   ```bash
   python3 -m venv tape_scaden_env
   ```

2. Activate the environment:
   ```bash
   source tape_scaden_env/bin/activate
   ```

3. Install required libraries:
   ```bash
   python3 -m pip install -r requirements_tape_scaden.txt
   ```

### Implementation

The R notebooks trigger bash scripts that:
1. Load the Python environment with required libraries
2. Execute Python scripts for TAPE and Scaden methods
3. Save results to CSV files that are read back into R

**Note:** TAPE and Scaden generate figures during execution that may pause the script. To avoid interruption, comment out the relevant plotting sections in the source code.

---

## Citation

If you use this code or the underlying methods in your research, please cite:

> Coudray A, Forey R, Bejar Haro B, Martins F, Carlevaro-Fita J, Sheppard S, Offner SE, La Manno G, Obozinski G, Trono D. *Computational deconvolution of drug sensitivity via single-cell–to–bulk transcriptome mapping.* Manuscript in revision, PLOS Computational Biology.

A preprint is available on bioRxiv: [doi.org/10.1101/2023.05.10.540140](https://www.biorxiv.org/content/10.1101/2023.05.10.540140v3) (posted under an earlier title, *Profiling drug sensitivity of leukemic stem cells via bulk-to-single-cell deconvolution* — the title above reflects the current manuscript).

*(Full citation with volume/DOI will be updated once the manuscript is published.)*

---

## License

This project is licensed under the MIT License - see the [LICENSE](license) file for details.

---

## Issues and Support

For bug reports, feature requests, or questions, please open an issue on the GitHub repository:

https://github.com/alexdray86/Coudray2026/issues
