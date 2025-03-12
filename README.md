# Novel RNA-binding protein YebC enhances translation of proline-rich amino acid stretches in bacteria

This repository contains the R scripts and source data necessary to reproduce the main results presented in the manuscript by [Ignatov *et al.* (2024)](https://doi.org/10.1101/2024.08.26.607280).

## Description of files and structure

We structured the analysis, R scripts, and source data based on the different high-throughput approaches applied in our study. Each analysis folder contains a `README.md` file with explanatory text about the input files and R scripts that were used to perform the analyses.

## Installation

Follow these steps to reproduce the analysis:

### Download repository

1. Clone or download this repository.
2. Navigate to the main project folder; it should contain the subfolder `analysis`.

```bash
git clone https://github.com/MPUSP/ignatov_et_al_2025.git
cd ignatov_et_al_2025
```

### Download genome and alignment data
1. Download data `data.zip` from [DRYAD](https://doi.org/10.5061/dryad.j0zpc86rg).
2. Unzip the downloaded data and move the `data` folder into the project's main folder.
3. Confirm that your project's folder now contains the subfolders `analysis` and `data`.
4. Inside the `data` folder, you should find three subfolders: `genome`, `iclip`, and `riboseq`.

Your project folder should look like this now:

```bash
ignatov_et_al_2025/
├── .gitignore
├── LICENSE
├── README.md
├── analysis/
│   └── 01_analysis_of_OOPS_and_RBS-ID
│   └── 02_analysis_of_RNAseq_for_yebC_mutant
│   └── 03_analysis_of_iCLIP_for_yebC
│   └── 04_analysis_of_riboseq_for_yebC_mutant
│   └── 05_analysis_of_proteomics_for_yebC_mutant
└── data/
    └── genome/
    └── iclip/
    └── riboseq/
```

### Setup R environment

In order to execute the R scripts, you will need the following R libraries:

- R (4.3.1)
- stringr (1.5.0)
- tidyverse (2.0.0)
- VennDiagram (1.7.3)
- beeswarm (0.4.0)
- ggpubr (0.6.0)
- ggrepel (0.9.4)
- plotROC (2.3.1)
- UniProt.ws (2.40.1)
- janitor (2.2.1)
- drawProteins (1.20.0)
- ggvenn (0.1.10)
- GenomicAlignments (1.36.0)
- edgeR (3.42.4)
- heatmaply (1.5.0)
- Biostrings (2.68.1)
- ggseqlogo (0.2.0)

## Authors

- Dr. Dmitriy Ignatov
  - Affiliation: [Max-Planck-Unit for the Science of Pathogens](https://www.mpusp.mpg.de/) (MPUSP), Berlin, Germany
  - ORCID profile: https://orcid.org/0000-0002-2237-974X
- Dr. Rina Ahmed-Begrich
  - Affiliation: [Max-Planck-Unit for the Science of Pathogens](https://www.mpusp.mpg.de/) (MPUSP), Berlin, Germany
  - ORCID profile: https://orcid.org/0000-0002-0656-1795

Visit the MPUSP github page at https://github.com/MPUSP for more information on other projects.

## Notes

If you use the source code or data provided in this repository, please cite our manuscript: [Ignatov *et al.* (2024)](https://doi.org/10.1101/2024.08.26.607280).