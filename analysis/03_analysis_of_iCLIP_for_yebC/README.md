# Analysis of iCLIP data for yebC

Description and analysis of source data to replicate the primary results of the iCLIP data analysis for yebC. The input files have to be downloaded from [DRYAD](https://doi.org/10.5061/dryad.j0zpc86rg) as explained in the projects main `README.md` file.

## 1. Input files

Input files are located in the downloaded folder from DRYAD `data/genome` and `data/iclip`:

**File:** `NC_002737.2.fa`

**Description:** *S. pyogenes* genome sequence file obtained from the NCBI database.

---
**File:** `220215_NC_002737.2_refseq_gene.gff`

**Description:** *S. pyogenes* genome annotation file downloaded from NCBI RefSeq on 25.02.2022.

---
**File:** `iCLIP-0316-1124_yebc_uv_1.bam`

**Description:** Mapping of iCLIP reads for sample `yebC:3FLAG_UV+_rep1`.

---
**File:** `iCLIP-0316-1124_yebc_uv_2.bam`
**Description:** Mapping of iCLIP reads for sample `yebC:3FLAG_UV+_rep2`.

---
**File:** `iCLIP-0316-1124_smi_1.bam`

**Description:** Mapping of iCLIP reads for sample `SMI_rep1`.

---
**File:** `iCLIP-0316-1124_smi_2.bam`

**Description:** Mapping of iCLIP reads for sample `SMI_rep2`.

---
**File:** `iCLIP-0316-1124_wt_uv.bam`

**Description:** Mapping of iCLIP reads for sample `WT_UV+`.

---
**File:** `iCLIP-0316-1124_yebC_no_uv.bam`

**Description:** Mapping of iCLIP reads for sample `yebC:3FLAG_UV-`.

## 2. Data analysis

**Script:** `01_read_bam.R`

**Description:** Uses the input files to assign the cross-linked (CL) sites to each nucleotide of the *S. pyogenes* genome:

1. Assigns CL sites to each genomic position.
2. Collapses CL sites of rRNAs to the first operon, accounting for the differences among the six rRNA operons.
3. Combines data for CL sites and coverage information, and normalizes to the total library size.
4. Writes coverage data on plus and minus strands to WIG files.
5. Writes cross-link data for plus and minus strands to WIG files.
6. Adds the data to the final data frame and saves it as `iclip.csv` and `iclip.tsv` files.

---

**Script:** `02_clustering_and_graphics.R`

**Description:** Uses the file `iclip.csv` as input:

1. Plots the coverage information for 23S rRNA and saves it as `23S_coverage.pdf`.
2. Clusters the CL sites located close to each other and saves the clusters to `cl_clusters.csv` and `cl_clusters.tsv`.
3. Plots the distribution of clusters in biotypes and the top cluster for each sample, and saves it as `barplot_cl_stat.pdf`.
