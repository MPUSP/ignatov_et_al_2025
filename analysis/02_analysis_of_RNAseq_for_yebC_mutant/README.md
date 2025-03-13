# RNA-seq analysis of the yebC mutant

Description and analysis of source data to replicate the primary results of the RNA-seq data analysis for the yebC mutant.

## 1. Input files
Input files are located in the folder `input_files`:

**File:** `gene_go_kegg_meta_info.tsv`

**Description:** Information about *S. pyogenes* genes, compiled from the NCBI's RefSeq database and the eggNOG Firmicutes database (version 5.0.2).

---
**File:** `all_tpm_counts.tsv`

**Description:** TPM (transcripts per million) gene expression values for each RNA-seq sample.

---

**File:** `log_d-vs-log_w.tsv`

**Description:** Differential expression (DE) analysis results of strain Δspy_0316 vs. WT (wild-type) in Log growth phase.

---

**File:** `log_d-vs-log_c.tsv`

**Description:** DE analysis results of strains Δspy_0316 vs. complement in Log growth phase.

---
**File:** `log_d-vs-log_wc.tsv`
**Description:** DE analysis results of strain Δspy_0316 vs. WT and complement in Log growth phase.

---

**File:** `stat_d-vs-stat_w.tsv`

**Description:** DE analysis results of strain Δspy_0316 vs. WT in Stat growth phase.

---
**File:** `stat_d-vs-stat_c.tsv`

**Description:** DE analysis results of strains Δspy_0316 vs. comp in Stat growth phase.

---
**File:** `stat_d-vs-stat_wc.tsv`


**Description:** DE analysis results of strain Δspy_0316 vs. WT and comp in Stat growth phase.

## 2. Data analysis

**Script:** `01-rna-seq-analysis.R`

**Description:** Uses all input files to prepare a single table with all RNA-seq DE results `rnaseq_results.tsv` and generates several figures:

1. Merges all the data into a single table containing TPM values and results of DE analysis, and identifies differentially expressed genes. Saves the table as `rnaseq_results.tsv`.
2. Generates correlation plots of differences in expression for the comparisons: `cor_log.pdf` and `cor_stat.pdf`.
3. Generates MA plots for DE analysis of Δspy_0316 vs. WT in Log and Stat growth phases. Highlights the DE genes with the largest expression difference or only the speB gene in: `ma_log.pdf`, `ma_stat.pdf`, and `ma_stat_speB.pdf`.
4. Generates Venn diagrams for the DE genes identified in the comparisons Δspy_0316 vs. WT and Δspy_0316 vs. comp in Log and Stat growth phases.
5. Saves the graphics to the `rnaseq_main.pdf` and `rnaseq_main_2.pdf` files.
