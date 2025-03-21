# Analysis of OOPS and RBS-ID data
Description and analysis of source data to reproduce the main results for the analysis of OOPS and RBS-ID datasets.

## 1. Input files
Input file are located in the folder `Raw data`:

**File:** `mainAnnot.streptococcus_pyogenes_serotype_m1.txt`

**Description:** Annotation of *S. pyogenes* proteome. File was downloaded from http://annotations.perseus-framework.org on 24.08.2021. This file should serve as the genome annotation in Perseus program, however there is an error while reading this file.

---
**File:** `Interpro.STRP1.tsv`

**Description:** Domain annotation from InterPro database.

---
**File:** `Domains_PNK_assay.csv”`

**Description:** Manual domain annotation for candidate RBPs selected for PNK assay.

---
**File:** `200519_20C011_results.xlsx`

**Description:** Results of OOPS experiment on *S. pyogenes* in Log and Stat growth phases.

---
**File:** `201016_20C028_additional_results.xlsx`

**Description:** Results of RBS-ID experiment on *S. pyogenes* in Log growth phase. The samples were prepared by two different protocols: original and OOPS followed by RBS-ID. The second variation was considered more successful and is considered as replicate 1 of RBS-ID experiment.

---
**File:** `21C021_2A_210823_interpro.txt`

**Description:** Replicates 2 and 3 of RBS-ID experiment.

---
**File:** `241120_RAB_spy_ignatov_other_studies.tsv`

**Description:** Conservation of *S. pyogenes* RBPs in other bacterial species. Data were obtained from the following publications:
- Shchepachev, V., *et al.* (2019). "Defining the RNA interactome by total RNA-associated protein purification." Mol Syst Biol 15(4): e8689. - 335 RBPs (*E. coli*)
- Queiroz, R. M. L., *et al.* (2019). "Comprehensive identification of RNA-protein interactions in any organism using orthogonal organic phase separation (OOPS)." Nat Biotechnol 37(2): 169-178. - 364 RBPs  (*E. coli*)
- Monti, M., *et al.* (2024). "Interrogation of RNA-protein interaction dynamics in bacterial growth." Mol Syst Biol. - 271 RBPs  (*E. coli*)
- Stenum, T. S., *et al.* (2023). "RNA interactome capture in Escherichia coli globally identifies RNA-binding proteins." Nucleic Acids Res 51(9): 4572-4587. - 170 RBPs  (*E. coli*)
- Urdaneta, E. C., *et al.* (2019). "Purification of cross-linked RNA-protein complexes by phenol-toluol extraction." Nat Commun 10(1): 990. - 172 RBPs (*Salmonella Typhimurium*)
- Chu, L. C., *et al.* (2022). "The RNA-bound proteome of MRSA reveals post-transcriptional roles for helix-turn-helix DNA-binding and Rossmann-fold proteins." Nat Commun 13(1): 2883. - 384 RBPs (*Staphylococcus aureus*)

## 2. Data preparation

**Prepared file:** `Proteome annotation with RBP.tsv`

**Description:** To clean `mainAnnot.streptococcus_pyogenes_serotype_m1.txt` in Excel do the following:
1. Deleted empty columns
2. Modified ENSG column so that Spy_XXXX locus tag comes first in all cells. Added Spy_XXXX to two cells, where it was missing.
3. Deleted P50470 (Immunoglobulin G-binding protein H) because could it not find the proper locus tag.
4. Saved the result as `Proteome annotation.tsv`
5. Used script `01-prepare-annotation.R` to add RBP annotation to `Proteome annotation.tsv` and saved the results as `Proteome annotation with RBP.tsv`.

---
**Prepared file:** `RBS-ID results.tsv`

**Description:**
1. Transferred the raw data of three RBS-ID replicates to `RBS-ID results.xlsx`:
    - `201016_20C028_additional_results.xlsx` (Rep_1)
    - `21C021_2A_210823_interpro.txt` (Rep_2)
    - `21C021_3A_210823_interpro.txt` (Rep_3).

2. In Rep_1 removed all rows from `Lys` sample.
3.	Merged all three replicates and created a column describing replicates. Saved as `RBS-ID results.tsv`.

## 3. Data analysis

**Script:** `02-rbsid-aas-sites-proteins.R`

**Description:** Uses `RBS-ID results.tsv` and `interpro.STRP1.tsv` as an input files:
1. Identifies cross-linked (CL) peptides and selects peptides with one CL.
2. Identifies CL proteins in each replicate and extracts their sequences from UniProt.
3. For each CL peptide identifies the CL amino acid and its coordinate within a protein.
4. Identifies the CL sites and sums up their spectral counts within replicates.
5. For each protein collects all detected CL sites.
6. Adds the domain information for CL sites
7. Plots Venn diagrams of CL sites and proteins in RBS-ID replicates: `Venn - CL sites/proteins in RBS-ID replicates.tiff`.
8. Calculates the CL frequencies of different amino acids and plots them along with CL numbers: `Crosslinked amino acids.pdf`, `Crosslinked amino acids in replicate 1/2/3.pdf`
9. Saves the files as `RBS-ID - CL sites.tsv` and `RBS-ID - CL proteins.tsv`.

---
**Script:** `03-merge-oops-rbsid-and-graphics.R`

**Description:** Uses `Proteome annotation with RBP.tsv`, `200519_20C011_results.tsv` and `RBS-ID - CL proteins.tsv` as input files:
1. Fixes the table with OOPS results.
2. Identified significantly enriched proteins in OOPS Log and Stat: fdr < 0.05 and >2-fold enrichment.
3. Merges Annotation, OOPS and RBS-ID data (RBP dataframe).
4. Identifies the candidate RBP: (i) not annotated as RBP; (ii) enriched in OOPS; (iii) CL amino acids in at least one RBS-ID replicate.
5. Saves the RBP dataframe as `Candidate RBPs.tsv`.
6. Calculates and plots the distribution of CL sites in PFAM domains: `CL sites in domains.tsv`, `CL sites in top domains.pdf` and `CL sites in domain groups.pdf`.
7. Generates Venn diagrams for OOPS enriched proteins in Log and Stat and for OOPS, RBS-ID and RBP annotation: `Venn - OOPS Log and Stat.tiff` and `Venn - RBP, OOPS, RBS-ID.tiff`.
8. Generates beeswarm plots of OOPS enrichment for total proteins, RBP types and RBS-ID CL proteins: `Beeswarm - OOPS enrichment of RBPs and CL.pdf`.
9. Generates density plots of OOPS scores for all proteins, annotated RBPs and proteins with identified CL sites: `Density plots of OOPS scores.pdf`.
10. Generates Volcano plots for OOPS enrichment of all RBPs (`Volcano - OOPS enrichment of RBPs.pdf`), different RBP types (`Volcano - OOPS enrichment of RBP types.pdf`) and proteins selected for PNK assay (`Volcano - OOPS enrichment for PNK assay.pdf`).
11. Plots a pie chart for functional categories of the candidate RBPs: `Pie chart - func cat.pdf`.
12. Calculates and plots ROC curves for RBP identification among all proteins and only OOPS enriched proteins: `ROC curve - all proteins.pdf` and `ROC curve - OOPS enriched proteins.pdf`.

---
**Script:** `04-plot-cl-domains.R`

**Description:** Uses `Candidate RBPs.tsv`, `interpro.STRP1.tsv` and `Domains.csv` as input files to plot domains and cross-linking sites for every protein with identified CL sites.
1. Using a function, domain annotation can be chosen.
2. The figures are generated for every candidate RBP and saved in the folder `CL sites for novel RBPs`.
3. The figures are generated for RBPs selected for PNK assay and saved in the folder `CL sites for PNK assay`. For the domain annotation uses `Domains.csv`.

---
**Script:** `05-rbp-conservation.R`

**Description:** Uses `241120_RAB_spy_ignatov_other_studies.tsv` and creates the final table with OOPS and RBS-ID results `OOPS_RBSID_orthologs.tsv`. Generates the figures `Number of RBPs in different studies.pdf` and `Annotated RBPs in different studies.pdf`.