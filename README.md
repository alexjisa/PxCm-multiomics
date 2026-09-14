# PxCm-multiomics

Analysis code accompanying the manuscript:

> **Temporal multi-omics dataset of a cucurbit-powdery mildew interaction**
> Jimenez-Sanchez A., Fernandez-Ortuno D., Pastor V., Polonio A., Perez-Garcia A.
> *Scientific Data* (in preparation, 2026).

## Overview

Temporal dual RNA-seq (host-pathogen) and untargeted LC-MS/MS metabolomics of
the *Cucumis melo* (melon) - *Podosphaera xanthii* (cucurbit powdery mildew)
interaction. Inoculated and mock-treated melon leaves were sampled at 12 time
points from 12 to 144 hours post-inoculation (hpi), and profiled by dual
RNA-seq and untargeted metabolomics from the same biological material. This
repository contains the R scripts used to process and explore that dataset:
differential expression, functional enrichment, temporal clustering, and
metabolomics QC/differential abundance analysis.

## Data availability

- **Raw RNA-seq (FASTQ)**: NCBI SRA, BioProject
  [PRJNA1454108](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1454108).
- **Raw LC-MS/MS (metabolomics)**: MetaboLights, accession
  [MTBLS15447](https://www.ebi.ac.uk/metabolights/MTBLS15447).
- **Processed data** (count/TPM matrices, differential expression results,
  functional enrichment, metabolomics tables and complete sample metadata):
  Zenodo, [DOI 10.5281/zenodo.19555639](https://doi.org/10.5281/zenodo.19555639).

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

## Dependencies (R)

`DESeq2`, `ggplot2`, `dplyr`, `tidyr`, `stringr`, `readr`, `readxl`,
`pheatmap`, `viridis`, `ggrepel`, `Rtsne`, `uwot`, `dendextend`, `scales`,
`svglite`, `limma`. Optional (functional enrichment): `clusterProfiler`,
`enrichplot`, `KEGGREST`, `httr`, `xml2`, `rvest`.

Exact package versions verified to run these scripts are listed in
[`session_info.txt`](session_info.txt) (R 4.4.2, matching the version cited
in the article's Methods).

## Repository status

This repository is currently under construction. The final analysis scripts
associated with the article are being prepared and will be uploaded to this
repository.

## License

MIT, see [LICENSE](LICENSE).

## Citation

> **Temporal multi-omics dataset of a cucurbit-powdery mildew interaction**
> Jimenez-Sanchez A., Fernandez-Ortuno D., Pastor V., Polonio A., Perez-Garcia A.
> *Scientific Data* (in preparation, 2026).

The full citation and DOI will be added upon publication.
