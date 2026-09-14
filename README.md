# PxCm-multiomics

Analysis code associated with the dataset and manuscript:

> **Temporal multi-omics dataset of a cucurbit-powdery mildew interaction**  
> Jiménez-Sánchez A., Fernández-Ortuño D., Pastor V., Polonio Á., Pérez-García A.  
> *Scientific Data* (in preparation, 2026).

This repository contains analysis code associated with a temporal dual RNA-seq and untargeted LC-MS/MS metabolomics dataset of the *Cucumis melo*–*Podosphaera xanthii* (cucurbit powdery mildew) interaction.

The experiment comprises 12 time points from 12 to 144 hours post-inoculation (hpi), providing temporal transcriptomic profiles of both host and pathogen together with metabolomic profiles of the infected host tissue.

Raw and processed datasets are deposited in public repositories as detailed below.

## Data availability

- **Raw RNA-seq data (FASTQ):** NCBI Sequence Read Archive (SRA), BioProject [PRJNA1454108](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1454108).
- **Raw LC-MS/MS data (vendor raw and CDF files):** MetaboLights, accession [MTBLS15447](https://www.ebi.ac.uk/metabolights/MTBLS15447).
- **Processed transcriptomics and metabolomics data:** Zenodo, [DOI: 10.5281/zenodo.19555639](https://doi.org/10.5281/zenodo.19555639).

The Zenodo repository includes count and TPM matrices, differential expression results, functional enrichment results, processed metabolomics tables, and sample metadata.

## Repository structure

```text
scripts/
  transcriptomics/
    Cm_transcriptomics.R
    Px_transcriptomics.R
    detected_genes_violin.R

  metabolomics/
    metabolomics_workflow.R
    pathway_enrichment.R
    top_dynamic_metabolites_heatmap.R
```

## Dependencies

Analyses were performed in R.

Main R packages include:

`DESeq2`, `ggplot2`, `dplyr`, `tidyr`, `stringr`, `readr`, `readxl`,
`pheatmap`, `viridis`, `ggrepel`, `Rtsne`, `uwot`, `dendextend`, `scales`,
`svglite`, and `limma`.

Optional packages used for functional enrichment include:

`clusterProfiler`, `enrichplot`, `KEGGREST`, `httr`, `xml2`, and `rvest`.

Exact package versions used to run the analyses are provided in
[`session_info.txt`](session_info.txt).

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
