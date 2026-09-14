# Data

This folder is empty in the Git repository. The processed input files used
by the scripts in `scripts/` are hosted on Zenodo:

**https://doi.org/10.5281/zenodo.19555639**

Raw sequencing/spectrometry data are hosted separately (see the main
[README](../README.md#data-availability)):
NCBI SRA BioProject [PRJNA1454108](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1454108)
(RNA-seq) and MetaboLights [MTBLS15447](https://www.ebi.ac.uk/metabolights/MTBLS15447)
(LC-MS/MS).

## Where to place downloaded files

Download the processed files from the Zenodo record above and place them
under `data/` in this layout:

```text
data/
├── transcriptomics/
│   ├── 01_Transcriptomics_Cm_counts.txt
│   ├── 02_Transcriptomics_Cm_sample_info.txt
│   ├── 03_Transcriptomics_Cm_TPM.txt
│   ├── 01_Transcriptomics_Px_counts.txt
│   ├── 02_Transcriptomics_Px_sample_info.txt
│   ├── 03_Transcriptomics_Px_TPM.txt
│   ├── 04_Px_Annotations_eggNOG-mapper.xlsx      (optional, see note below)
│   └── Cm_functional_annotation.tsv              (optional, see note below)
│
└── metabolomics/
    ├── 002_METABO_NEG_AVG.txt
    ├── 002_METABO_POS_AVG.txt
    ├── 003_METABO_COMB_AVG.txt
    ├── 14_Metabolomics_Sets.xlsx
    └── 15_Metabolomics_Significant.xlsx
```

`data/transcriptomics/` and `data/metabolomics/` are git-ignored (see
`.gitignore`) - only this README is tracked. Create the two subfolders
yourself when you download the data.

## Which script reads which file

See [`docs/data_dictionary.md`](../docs/data_dictionary.md) for a full
description of each file's contents. Summary of what each script expects:

| File | Used by |
|---|---|
| `01_Transcriptomics_Cm_counts.txt` | `scripts/transcriptomics/01_Cm_transcriptomics.R` |
| `02_Transcriptomics_Cm_sample_info.txt` | `scripts/transcriptomics/01_Cm_transcriptomics.R` |
| `03_Transcriptomics_Cm_TPM.txt` | `scripts/transcriptomics/03_detected_genes.R` |
| `01_Transcriptomics_Px_counts.txt` | `scripts/transcriptomics/02_Px_transcriptomics.R` |
| `02_Transcriptomics_Px_sample_info.txt` | `scripts/transcriptomics/02_Px_transcriptomics.R` |
| `03_Transcriptomics_Px_TPM.txt` | `scripts/transcriptomics/03_detected_genes.R` |
| `04_Px_Annotations_eggNOG-mapper.xlsx` | `scripts/transcriptomics/02_Px_transcriptomics.R` (optional, KEGG enrichment) |
| `Cm_functional_annotation.tsv` | `scripts/transcriptomics/01_Cm_transcriptomics.R` (optional, GO/KEGG enrichment - see note) |
| `002_METABO_NEG_AVG.txt` | `scripts/metabolomics/01_metabolomics_workflow.R` |
| `002_METABO_POS_AVG.txt` | `scripts/metabolomics/01_metabolomics_workflow.R` |
| `003_METABO_COMB_AVG.txt` | `scripts/metabolomics/01_metabolomics_workflow.R` |
| `14_Metabolomics_Sets.xlsx` | `scripts/metabolomics/02_pathway_enrichment.R` |
| `15_Metabolomics_Significant.xlsx` | `scripts/metabolomics/03_dynamic_metabolites.R` |

## Known issue: `Cm_functional_annotation.tsv`

`01_Cm_transcriptomics.R` expects a **tab-separated** file with columns
`Geneid`, `GO`, `KEGG`. The eggNOG-mapper annotation actually produced for
*C. melo* (`04_Cm_Annotations_eggNOG-mapper.xlsx`) is an **Excel** file whose
columns are named `Name` (not `Geneid`), `GOs` (not `GO`), and `KEGG_ko`
(not `KEGG`). As shipped, the script's GO/KEGG enrichment step will not run
against that file. This is flagged, not fixed, in this reorganization - see
the main [README](../README.md#repository-status).
