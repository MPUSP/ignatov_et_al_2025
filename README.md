# Novel RNA-binding protein YebC enhances translation of proline-rich amino acid stretches in bacteria

This repository contains the R scripts and source data to reproduce the main results presented in the manuscript [Ignatov *et al.*, 2024](https://doi.org/10.1101/2024.08.26.607280).

## Description of files and structure

We structured the analysis and R scripts based on different high-throughput approaches applied. Each analysis folder contains a `README.md` file with explanatory text about the input files and how to perform the analyses.

## Installation

Follow these steps to reproduce the analysis:

### Download Repository

1. Clone or download this repository.
2. Navigate to the main project folder; it should contain the subfolders `analysis`

```bash
git clone https://github.com/MPUSP/ignatov_et_al_2025.git
cd ignatov_et_al_2025
```

### Download Alignment Data
1. Download alignment data `alignment_data.zip` from [DRYAD](https://doi.org/10.5061/dryad.j0zpc86rg).
2. Unzip the downloaded data and move the 'alignment_data' folder into the project's main folder.
3. Confirm that your project's folder now contains the subfolders: `analysis`, and `alignment_data`.
4. Inside the `alignment_data` folder, you should find two subfolders, such as `riboseq_data`, and `iclip_data`.

### Setup R Environment

In order to execute the R scripts, you will need to following R libraries:

- R (4.3.1)
- stringr ()
- tidyverse ()
- VennDiagram ()
- beeswarm ()
- ggpubr ()
- ggrepel ()
- plotROC ()
- UniProt.ws ()
- janitor ()
- drawProteins ()
- ggvenn ()
- GenomicAlignments ()
- edgeR ()
- heatmaply ()
- biostrings ()
- ggseqlogo ()

## Authors

- Dr. Dmitriy Ignatov
  - Affiliation: [Max-Planck-Unit for the Science of Pathogens](https://www.mpusp.mpg.de/) (MPUSP), Berlin, Germany
  - ORCID profile: https://orcid.org/0000-0002-2237-974X
- Dr. Rina Ahmed-Begrich
  - Affiliation: [Max-Planck-Unit for the Science of Pathogens](https://www.mpusp.mpg.de/) (MPUSP), Berlin, Germany
  - ORCID profile: https://orcid.org/0000-0002-0656-1795

Visit the MPUSP github page at https://github.com/MPUSP for more information on other projects.

## Notes

If you use the data or source code provided in this repository, please cite our manuscript: [Ignatov *et al.*, 2024](https://doi.org/10.1101/2024.08.26.607280).