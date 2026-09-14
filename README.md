# PxCm-multiomics

Analysis code associated with the dataset and manuscript:

> **Temporal multi-omics dataset of a cucurbit-powdery mildew interaction**  
> Jiménez-Sánchez A., Fernández-Ortuño D., Pastor V., Polonio Á., Pérez-García A.  
> *Scientific Data* (in preparation, 2026).

## Overview

This repository contains analysis code associated with a temporal dual
RNA-seq and untargeted LC-MS/MS metabolomics dataset of the *Cucumis
melo*–*Podosphaera xanthii* (cucurbit powdery mildew) interaction.

The experiment comprises 12 time points from 12 to 144 hours post-inoculation
(hpi), providing temporal transcriptomic profiles of both host and pathogen
together with metabolomic profiles of the infected host tissue. The scripts
here cover differential expression, functional enrichment, temporal
clustering, and metabolomics QC/differential abundance analysis.

Raw and processed datasets are deposited in public repositories as detailed
below.

## Data availability

- **Raw RNA-seq data (FASTQ):** NCBI Sequence Read Archive (SRA), BioProject [PRJNA1454108](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1454108).
- **Raw LC-MS/MS data (vendor raw and CDF files):** MetaboLights, accession [MTBLS15447](https://www.ebi.ac.uk/metabolights/MTBLS15447).
- **Processed transcriptomics and metabolomics data:** Zenodo, [DOI: 10.5281/zenodo.19555639](https://doi.org/10.5281/zenodo.19555639).

The Zenodo repository includes count and TPM matrices, differential
expression results, functional enrichment results, processed metabolomics
tables, and sample metadata.

## Repository structure

```text
PxCm-multiomics/
├── README.md
├── LICENSE
├── .gitignore
├── session_info.txt
│
├── scripts/
│   ├── transcriptomics/
│   │   ├── 01_Cm_transcriptomics.R   # C. melo: TPM, QC, DESeq2, clustering, GO/KEGG
│   │   ├── 02_Px_transcriptomics.R   # P. xanthii: TPM, QC, DESeq2, clustering, KEGG
│   │   └── 03_detected_genes.R       # Detected-genes (TPM>=1) violin plot
│   │
│   └── metabolomics/
│       ├── 01_metabolomics_workflow.R  # QC, PCA/t-SNE/UMAP, DAM analysis
│       ├── 02_pathway_enrichment.R     # KEGG pathway coverage plot
│       └── 03_dynamic_metabolites.R    # Top 20 dynamic-features heatmap
│
├── data/
│   └── README.md                     # Where to download inputs, expected layout
│
└── docs/
    └── data_dictionary.md            # Per-file column/content description
```

Each script writes its outputs under `results/transcriptomics/` or
`results/metabolomics/` (git-ignored; created automatically when a script
runs).

## How to obtain the input data

The processed input files are not stored in this repository. Download them
from the Zenodo record above and place them under `data/` as described in
[`data/README.md`](data/README.md). [`docs/data_dictionary.md`](docs/data_dictionary.md)
documents the contents of each file and which script reads it.

## How to run the analyses

Run every script from the project root, e.g.:

```bash
Rscript scripts/transcriptomics/01_Cm_transcriptomics.R
Rscript scripts/transcriptomics/02_Px_transcriptomics.R
Rscript scripts/transcriptomics/03_detected_genes.R
Rscript scripts/metabolomics/01_metabolomics_workflow.R
Rscript scripts/metabolomics/02_pathway_enrichment.R
Rscript scripts/metabolomics/03_dynamic_metabolites.R
```

Each script's header documents its exact inputs, outputs, and usage. Within
each `scripts/transcriptomics/` and `scripts/metabolomics/` folder, scripts
are independent of each other (none reads another script's output).

## Dependencies

Analyses were performed in R.

Main R packages include:

`DESeq2`, `ggplot2`, `dplyr`, `tidyr`, `stringr`, `readr`, `readxl`,
`pheatmap`, `viridis`, `ggrepel`, `Rtsne`, `uwot`, `dendextend`, `scales`,
`svglite`, and `limma`.

Optional packages used for functional enrichment include:

`clusterProfiler`, `enrichplot`, `KEGGREST`, `httr`, `xml2`, and `rvest`.

Exact package versions used to run the analyses are provided in
[`session_info.txt`](session_info.txt) (R 4.4.2, matching the version cited
in the article's Methods).

## Repository status

This repository is currently under construction. The final analysis scripts associated with the article are being prepared and will be uploaded to this repository.

## License

This repository is distributed under the MIT License. See [LICENSE](LICENSE) for details.

## Cite us

If you use the code or any of the associated datasets, please cite the corresponding article:

> **Temporal multi-omics dataset of a cucurbit-powdery mildew interaction**  
> Jiménez-Sánchez A., Fernández-Ortuño D., Pastor V., Polonio Á., Pérez-García A.  
> *Scientific Data* (in preparation, 2026).

The full citation and DOI will be added upon publication.
