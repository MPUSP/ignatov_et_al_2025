# Proteomics analysis of the yebC mutant
Description and analysis of source data to replicate the primary results of the mass spectrometry (MS) proteomics data analysis for the yebC mutant.

## 1. Input files
Input files are located in the folder `input_files`:

**File:** `24C003_DVI_df_ProteinLevelData.tsv`

---
**File:** `24C003_DVI_df_comparisonResult.tsv`

---
**File:** `Proteome annotation.tsv`

---
**File:** `rnaseq_results.tsv`

## 2. Data analysis

**Script:** `01_combine_data.R`

**Description:** Reads all input files, identifies differentially expressed genes, and creates a combined table of MS results: `delta0316_ms.tsv.`

---
**Script:** `02_de_graphics.R`

**Description:** Plots MA plots for differentially expressed proteins and saves them as `ma_plot_log.pdf` and `ma_plot_stat.pdf`.

---
**Script:** `03_pausing_effect_on_expression.R`

**Description:** Uses `delta0316_ms.tsv` and `rnaseq_results.tsv` to compare the expression of genes with yebC-dependent pausing sites to the rest of the genes at the protein and RNA levels. Generates boxplots (`boxplot_prot_DE.pdf`, `boxplot_rna_DE.pdf`, and `boxplot_translation_DE.pdf`) and calculates statistics for the difference in expression levels between these two gene groups.