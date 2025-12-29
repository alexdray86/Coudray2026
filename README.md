# Computational deconvolution of drug sensitivity via single-cell–to–bulk transcriptome mapping

This repository contains code to reproduce all figures from Coudray et al. 2026.

## ClimbTheCliff R Package

The main methods developed in this work (CLIMB and CLIFF) are available as an R package:

**Repository:** https://github.com/alexdray86/ClimbTheCliff

### Quick Installation

```r
# Using devtools
devtools::install_github("alexdray86/ClimbTheCliff", force=T)

# Or using remotes
remotes::install_github("alexdray86/ClimbTheCliff", force=T)
```

### Methods Overview

- **CLIMB** (Cell-type abundance and expression deconvolution using Inverse Modeling of Bulk RNA-seq): Deconvolutes bulk RNA-seq data using single-cell RNA-seq references to predict cell-type proportions and high-resolution cell-type-specific gene expression for each sample.

- **CLIFF** (Cell-type-specific drug sensitivity from deconvoluted expression profiles): Integrates CLIMB deconvolution results with drug sensitivity data (AUC values) to identify cell-type-specific drug responses, optionally incorporating somatic mutation data.

See the [ClimbTheCliff README](https://github.com/alexdray86/ClimbTheCliff) for detailed documentation and examples.

## Setup 

All notebooks are in R (requires R > v4.0) and were run with the IRkernel for Jupyter notebook. Libraries needed for each notebook are listed in the first cell and should be installed before running. Installation commands are provided in extra cells within the notebooks.

## Data: IMPORTANT 

All required data are stored in a publicly available Google Drive: 

**Download link:** https://drive.google.com/file/d/1mL-ZxlkETq9EzZ85TP6Jmtx54NNv1cx1/view?usp=sharing

Unzip `data.zip` and place the `data/` folder at the root of this directory. 

The `data/` folder contains:
- Raw single-cell datasets
- Pre-processed data (e.g., deconvolution results)
- Intermediate results for analysis without re-running all deconvolution methods
 

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
 

## Notebook Descriptions

Each notebook corresponds to a specific figure in the manuscript. Below is a summary of each notebook's contents:

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

### Fig4_expression_deconvolution.ipynb

Two deconvolution scenarios:
1. Cell subtype expression deconvolution from in-vitro experiment (bulk RNA-seq ground truth) and cross-dataset context (Van Galen → Naldini)
2. Over-expression procedure: 20 genes over-expressed in HSC-like cell subtype, guided by Van Galen → Naldini analysis

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

