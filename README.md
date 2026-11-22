# Bioinformatics Pipeline for Gene Expression Analysis

A Python-based bioinformatics pipeline developed as part of a master's thesis on gene expression changes in ovarian cancer in the context of chemoresistance.

The tool automates the statistical analysis of gene expression data obtained from the **Gene Expression Omnibus (GEO)** database, enabling comparisons between two groups of samples.

## Features

- **Single-gene analysis** – compares gene expression levels individually using the Mann–Whitney U test or Welch's t-test.
- **Gene signature analysis** – evaluates multivariate differences between sample groups using PERMANOVA or MANOVA.
- **Automatic statistical test selection** – selects the appropriate test based on sample size and normality assessment.
- **Multiple testing correction** – applies Benjamini–Hochberg FDR and Bonferroni correction in single-gene analysis.
- **Effect size estimation** – calculates measures such as log2 fold change, Hedges' g, eta-squared and rank-biserial correlation.
- **Automated result interpretation** – classifies statistical significance, effect size and direction of expression changes.
- **CSV reporting** – exports standardized results and, when applicable, a separate normality assessment report.
- **Command-line interface** – allows users to configure and run analyses directly from the terminal.

## Tech Stack

- **Language:** Python
- **Data processing:** pandas, NumPy
- **Statistical analysis:** SciPy, statsmodels, Pingouin
- **Multivariate analysis:** scikit-bio
- **Command-line interface:** argparse

## Analysis Workflow

1. **Input:** Loads gene expression data from a structured CSV file.
2. **Preprocessing:** Removes empty rows and columns and optionally applies log2(x + 1) transformation.
3. **Group definition:** Divides samples into two groups based on the specified number of samples in the first group.
4. **Statistical analysis:** Automatically selects and performs the appropriate statistical tests.
5. **Result interpretation:** Calculates statistical measures, applies multiple testing correction where applicable, and generates descriptive interpretations.
6. **Output:** Saves the results in CSV format.

## Statistical Methods

| Analysis mode | Statistical tests | Additional methods |
|---|---|---|
| Single-gene | Mann–Whitney U, Welch's t-test | Shapiro–Wilk, FDR, Bonferroni |
| Gene signature | PERMANOVA, MANOVA | PERMDISP, Box's M |

## Input Data Format

The pipeline accepts CSV files with a semicolon-separated structure and decimal commas.

Required columns:

- `sonda` – probe identifier
- `gen` – gene symbol
- `GSM_ID` columns – gene expression values for individual samples

Example:

```csv
sonda;gen;GSM_ID_1;GSM_ID_2;GSM_ID_3
33824_at;KRT8;5,21;4,87;6,12
201596_x_at;KRT18;7,34;6,92;7,01
```

The input data should be pre-normalized and log-transformed where appropriate. An optional log2(x + 1) transformation is also available.

## Usage

Run the script from the terminal:

```bash
python program.py \
  --input data.csv \
  --output results.csv \
  --group1_size 10 \
  --transform none \
  --mode single_gene
```

### Command-line arguments

| Argument | Description |
|---|---|
| `-i`, `--input` | Input CSV file |
| `-o`, `--output` | Output CSV file |
| `--group1_size` | Number of samples assigned to the first group |
| `--transform` | Data transformation: `none` or `log2` |
| `--mode` | Analysis mode: `single_gene` or `signature` |

## Project Context

This pipeline was developed as part of the master's thesis:

*Analysis of the Expression of Selected Keratin and EMT-Related Genes in Ovarian Cancer in the Context of Chemoresistance Using Data from the Gene Expression Omnibus Database and an In-house Developed Bioinformatics Pipeline.*

The tool was used to investigate gene expression patterns associated with epithelial–mesenchymal transition (EMT) and chemoresistance in ovarian cancer datasets.

**Project type:** Master's thesis / Bioinformatics / Statistical Data Analysis
