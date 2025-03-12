# Ribo-seq analysis of the yebC mutant

Description and analysis of source data to replicate the primary results of the ribo-seq data analysis for the yebC mutant. The input files have to be downloaded from [DRYAD](https://doi.org/10.5061/dryad.j0zpc86rg) as explained in the projects main `README.md` file.

## 1. Input files

Input files are located in the downloaded folder from DRYAD `data/genome` and `data/riboseq`:

**File:** `NC_002737.2.fa`

**Description:** *S. pyogenes* genome sequence file obtained from the NCBI database.

---
**File:** `220215_NC_002737.2_refseq_gene.gff`

**Description:** *S. pyogenes* genome annotation file downloaded from NCBI RefSeq on 25.02.2022.

---
**File:** `riboseq-wt_1.bam`

**Description:** Mapping of Ribo-seq reads for sample `WT_rep1`.

---

**File:** `riboseq-wt_2.bam`

**Description:** Mapping of Ribo-seq reads for sample `WT_rep2`.

---
**File:** `riboseq-wt_3.bam`

**Description:** Mapping of Ribo-seq reads for sample `WT_rep3`.

---
**File:** `riboseq-d0316_1.bam`

**Description:** Mapping of Ribo-seq reads for sample `ΔyebC_rep1`.

---
**File:** `riboseq-d0316_2.bam`

**Description:** Mapping of Ribo-seq reads for sample `ΔyebC_rep2`.

---
**File:** `riboseq-d0316_3.bam`

**Description:** Mapping of Ribo-seq reads for sample `ΔyebC_rep3`.

---
**File:** `riboseq-comp_1.bam`

**Description:** Mapping of Ribo-seq reads for sample `ΔyebC / yebC+_rep1`.

---
**File:** `riboseq-comp_2.bam`

**Description:** Mapping of Ribo-seq reads for sample `ΔyebC / yebC+_rep2`.

---
**File:** `riboseq-comp_3.bam`

**Description:** Mapping of Ribo-seq reads for sample `ΔyebC / yebC+_rep3`.

---
**File:** `riboseq-m6_1.bam`

**Description:** Mapping of Ribo-seq reads for sample `ΔyebC / yebC+_M2+_rep1`.

---
**File:** `riboseq-m6_1.bam`

**Description:** Mapping of Ribo-seq reads for sample `ΔyebC / yebC+_M2+_rep1`.

---
**File:** `riboseq-m6_1.bam`

**Description:** Mapping of Ribo-seq reads for sample `ΔyebC / yebC+_M2+_rep1`.

## 2. Data analysis

**Script:** `01-read-bam.R`

**Description:**
1. Processes Ribo-seq data and calculates the ribosome P-site occupancy for each nucleotide of coding sequences (CDSs), and saves the results to `riboseq.csv`.
2. Creates IGV files for visualization in the Integrative Genomics Viewer (IGV) and saves them to `igv` folder.

---
**Script:** `02-read-bam_A_site.R`

**Description:** Processes Ribo-seq data and calculates the ribosome A-site occupancy for each nucleotide of CDSs, and saves the results to `riboseq_A_site.csv`.

---
**Script:** `03-read-bam_E_site.R`

**Description:** Processes Ribo-seq data and calculates the ribosome E-site occupancy for each nucleotide of CDSs, and saves the results to `riboseq_E_site.csv`.

---
**Script:** `04-plot-coverage-start-stop.R`

**Description:**
1. Uses `riboseq.csv` to plot ribosome occupancy around start and stop codons.
2. Plots the ribosome occupancy around the start codon and saves it as `coverage_start_codon.pdf`.

---
**Script:** `05-calculate-pausing.R`

**Description:** Uses `riboseq.csv` to calculate pausing scores and identify amino acids with statistically significant difference in pausing between samples with functional YebC protein and those without it.

1. Calculates the pausing scores of amino acids in CDSs with sufficient sequencing coverage.
2. Uses the edgeR software package to identify amino acids with statistically significant differences in pausing between samples with functional YebC protein (WT and ΔyebC/yebC+) and those without it (ΔyebC and ΔyebC/yebC_M2+).
3. Applies the following criteria to identify amino acids with significant changes in pausing: (i) FDR-adjusted p value < 1e-03, (ii) fold change in pausing score > 4, and (iii) average pausing score in the samples with increased pausing > 2.
4. Saves the amino acid positions with changes in pausing to `riboseq_cds_diff.csv` and `riboseq_cds_diff.tsv`.

---
**Script:** `06-plot-pausing-sites.R`

**Description:** Uses `riboseq_cds.csv` to plot pausing scores for loci with significant difference in pausing scores. The plots are saved to the folder `gene_pausing_figures`.

---
**Script:** `07-plot-pausing-consensus.R`

**Description:** Uses `riboseq_cds.csv` to analyse the data and plot several figures:
1. Extracts sequences around the differential pausing sites and saves the results to `pausing_seq.fa`. Plots the consensus around the differential pausing sites as WebLogo plots, saving them as `weblogo_bits.pdf` and `weblogo_prob.pdf`.
2. Clusters the neighbouring pausing sites and saves the results to `pausing_clusters.tsv`. Plots the number of different motifs in the clustered pausing sites and saves the plot as `pausing_motifs_pie_chart.pdf`.
3. Plots the distribution of pausing scores for each of the 20 amino acids and saves the plot as `scores_amino_acids.pdf`.
4. Calculates and plots the pausing scores around the pausing motifs `DIP`, `PIP` and `PPG` and saves the plot as `DIP_PIP_PPG_scores.pdf`.
