# PxCm-multiomics

Analysis code associated with the dataset and article:

> **Temporal multi-omics dataset of a cucurbit-powdery mildew interaction**
> Jimenez-Sanchez A., Fernandez-Ortuno D., Pastor V., Polonio A., Perez-Garcia A.
> *Scientific Data* (in preparation / 2026).

Temporal dual RNA-seq (host-pathogen) and untargeted LC-MS/MS metabolomics
resources for the *Cucumis melo* - *Podosphaera xanthii* (cucurbit powdery
mildew) interaction, sampled across 12 time points between 12 and 144 hours
post-inoculation.

This repository contains code only. Raw and processed data are deposited in
public repositories, see below.

## Data availability

- **Raw RNA-seq (FASTQ)**: NCBI SRA, BioProject
  [PRJNA1454108](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1454108).
- **Raw LC-MS/MS (metabolomics)**: MetaboLights, accession
  [MTBLS15447](https://www.ebi.ac.uk/metabolights/MTBLS15447).
- **Processed data** (count/TPM matrices, differential expression results,
  functional enrichment, metabolomics tables and complete sample metadata):
  Zenodo, [DOI 10.5281/zenodo.19555639](https://doi.org/10.5281/zenodo.19555639).

## Repository structure

```
scripts/
  transcriptomics/
    Cm_transcriptomics.R          # C. melo: TPM, PCA/t-SNE/UMAP, DESeq2, temporal clustering, GO/KEGG
    Px_transcriptomics.R          # P. xanthii: TPM, PCA/t-SNE/UMAP, DESeq2, temporal clustering, KEGG
    detected_genes_violin.R       # Detected genes (TPM>=1) per group: Cm Mock/Inoculated and Px
  metabolomics/
    metabolomics_workflow.R       # QC/exploratory: PCA, t-SNE/UMAP, sample-distance heatmaps, DAM analysis (limma)
    pathway_enrichment.R          # KEGG pathway coverage (from MarVis-Pathway output)
    top_dynamic_metabolites_heatmap.R  # Heatmap of the 20 most dynamic features
```

Each script assumes it is run with the working directory set to the folder
containing the corresponding processed data files (downloaded from Zenodo,
see above). Input files expected by each script:

| Script | Input files |
|---|---|
| `Cm_transcriptomics.R` | `01_Transcriptomics_Cm_counts.txt`, `02_Transcriptomics_Cm_sample_info.txt` |
| `Px_transcriptomics.R` | `01_Transcriptomics_Px_counts.txt`, `02_Transcriptomics_Px_sample_info.txt`, `Px_functional_annotation.xlsx` (optional, eggNOG-mapper annotation for KEGG enrichment) |
| `detected_genes_violin.R` | `03_Transcriptomics_Cm_TPM.txt`, `03_Transcriptomics_Px_TPM.txt` |
| `metabolomics_workflow.R` | `data/processed/002_METABO_NEG_AVG.txt`, `data/processed/002_METABO_POS_AVG.txt`, `data/processed/003_METABO_COMB_AVG.txt` |
| `pathway_enrichment.R` | `14_Metabolomics_Sets.xlsx` |
| `top_dynamic_metabolites_heatmap.R` | `15_Metabolomics_Significant.xlsx` |

`Cm_transcriptomics.R` looks for the input files in `data/processed/` if that
folder exists under the working directory, and otherwise in the working
directory itself. `metabolomics_workflow.R` is called as
`Rscript scripts/metabolomics/metabolomics_workflow.R` from the project root
and always reads its inputs from `data/processed/` and writes outputs under
`results/`. The other scripts look for their inputs in the working directory.

## Methods (summary)

- **Transcriptomics**: FastQC (v0.11.9) + MultiQC (v1.27) for quality
  control; Trimmomatic (v0.39) for adapter/quality trimming; HISAT2 (v2.2.1)
  for alignment against the *C. melo* genome assembly v4.0 (mock samples)
  and, for inoculated samples, sequential alignment against the *C. melo*
  genome followed by the *P. xanthii* isolate 2086 genome; SAMtools (v1.23.1)
  for alignment processing; featureCounts (v2.0.6) for read counting; DESeq2
  (v1.46.0) for differential expression (padj < 0.05, |log2FC| > 1);
  functional enrichment with DAVID (v6.8, KEGG, *C. melo*) and
  clusterProfiler (v4.14.6, GO via eggNOG-mapper v2, *P. xanthii*).
- **Metabolomics**: LC-MS/MS acquisition (ESI+/ESI-) with MassLynx (v4.2);
  processing with the `xcms` package (centWave peak detection, retention
  time correction, peak grouping, missing-value imputation); filtering and
  tentative metabolite identification with MarVis-Suite 2.0 (MarVis-Filter,
  MarVis-Pathway); differential features by one-way ANOVA with FDR
  correction (Benjamini-Hochberg, p.adj < 0.05).
- Analysis performed in R (v4.4.2).

Full parameter details are in the Methods section of the article.

## Dependencies (R)

`DESeq2`, `ggplot2`, `dplyr`, `tidyr`, `stringr`, `readr`, `readxl`,
`pheatmap`, `viridis`, `ggrepel`, `Rtsne`, `uwot`, `dendextend`, `scales`,
`svglite`, `limma`. Optional (functional enrichment): `clusterProfiler`,
`enrichplot`, `KEGGREST`, `httr`, `xml2`, `rvest`.

Exact package versions verified to run these scripts are listed in
[`session_info.txt`](session_info.txt).

## Repository status

This repository is under construction. It currently includes the
DESeq2/clustering/enrichment transcriptomics processing, a metabolomics
QC/exploratory workflow (PCA, sample-distance heatmaps, t-SNE/UMAP, and
differential accumulated metabolite [DAM] analysis), and two metabolomics
visualization scripts (pathway coverage and temporal-dynamics heatmap).

Note: `metabolomics_workflow.R` identifies DAMs with `limma` (moderated
t-statistics on median-normalized, log10-transformed, Pareto-scaled
intensities), which differs from the one-way ANOVA with Benjamini-Hochberg
FDR correction described in the article's Methods section. This has not yet
been reconciled with the manuscript text.

## License

MIT, see [LICENSE](LICENSE).

## Citation

If you use this code or the associated dataset, please cite the article (see
the header of this README) and the corresponding data repositories
(BioProject PRJNA1454108, MetaboLights MTBLS15447, Zenodo DOI
10.5281/zenodo.19555639).
