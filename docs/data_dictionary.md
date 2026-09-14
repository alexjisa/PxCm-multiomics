# Data dictionary

Describes the processed input files consumed by the scripts in this
repository (see [`data/README.md`](../data/README.md) for where to obtain
and place them). Row/column counts and column names below were verified
directly against the actual files and cross-checked against the figures
reported in the article's Technical Validation section where noted; they are
not inferred.

## Transcriptomics (`data/transcriptomics/`)

### `01_Transcriptomics_Cm_counts.txt`
- **Omics layer / organism**: RNA-seq, *Cucumis melo* (host).
- **Rows**: 28,298 genes (melon genome assembly v4.0 gene models, e.g. `MELO3C017082.2`).
- **Columns**: `Geneid`, `Chr`, `Start`, `End`, `Strand`, `Length` (featureCounts annotation columns, `;`-delimited for multi-exon genes), followed by 72 sample columns named `<time>h_<condition>_R<replicate>` (12 time points x Mock/Inoculated x 3 replicates).
- **Values**: raw integer read counts (featureCounts).
- **Used by**: `scripts/transcriptomics/01_Cm_transcriptomics.R` (DESeq2 input; TPM is computed from this file, not read from `03_Transcriptomics_Cm_TPM.txt`).

### `02_Transcriptomics_Cm_sample_info.txt`
- **Rows**: 72 samples (one row per RNA-seq library).
- **Columns**: `sample`, `sample_name`, `condition` (`Mock`/`Inoculated`), `time_point` (`12h`-`144h`), `replicate` (`1`-`3`).
- **Used by**: `scripts/transcriptomics/01_Cm_transcriptomics.R`.

### `03_Transcriptomics_Cm_TPM.txt`
- Same row/column structure as `01_Transcriptomics_Cm_counts.txt` (28,298 genes x same 72 samples).
- **Values**: TPM (transcripts per million).
- **Used by**: `scripts/transcriptomics/03_detected_genes.R` (counts genes with TPM >= 1 per sample).

### `01_Transcriptomics_Px_counts.txt`
- **Omics layer / organism**: RNA-seq, *Podosphaera xanthii* (pathogen, isolate 2086).
- **Rows**: 16,355 genes (contig-based gene IDs, e.g. `00g000010`).
- **Columns**: `Geneid`, `Chr`, `Start`, `End`, `Strand`, `Length`, followed by 36 sample columns `<time>h_I_R<replicate>` (12 time points x 3 replicates; Inoculated only - the fungus is an obligate biotroph and is absent from mock-treated leaves).
- **Values**: raw integer read counts.
- **Used by**: `scripts/transcriptomics/02_Px_transcriptomics.R`.

### `02_Transcriptomics_Px_sample_info.txt`
- **Rows**: 36 samples.
- **Columns**: `sample`, `sample_name`, `condition` (always `Inoculated`), `time_point`, `replicate`, `infection_phase` (`Early`/`Intermediate`/`Late`, a categorical grouping of time points used for the phase-based DESeq2 contrasts).
- **Used by**: `scripts/transcriptomics/02_Px_transcriptomics.R`.

### `03_Transcriptomics_Px_TPM.txt`
- Same row/column structure as `01_Transcriptomics_Px_counts.txt` (16,355 genes x same 36 samples).
- **Values**: TPM.
- **Used by**: `scripts/transcriptomics/03_detected_genes.R`.

### `04_Px_Annotations_eggNOG-mapper.xlsx`
- **Omics layer / organism**: functional annotation, *P. xanthii* (eggNOG-mapper v2 output).
- **Rows**: 10,605 genes with an eggNOG hit (subset of the 16,355 total *P. xanthii* genes).
- **Columns**: `Geneid`, `seed_ortholog`, `evalue`, `score`, `eggNOG_OGs`, `max_annot_lvl`, `COG_category`, `Description`, `Preferred_name`, `GOs`, `EC`, `KEGG_ko`, `KEGG_Pathway`, `KEGG_Module`, `KEGG_Reaction`, `KEGG_rclass`, `BRITE`, `KEGG_TC`, `CAZy`, `BiGG_Reaction`, `PFAMs`.
- **Used by**: `scripts/transcriptomics/02_Px_transcriptomics.R` (optional; only the `Geneid` and `KEGG_ko` columns are used, for KEGG term enrichment of each temporal cluster).

### `Cm_functional_annotation.tsv` (optional, not currently usable as coded - see below)
- `scripts/transcriptomics/01_Cm_transcriptomics.R` expects a **tab-separated** file with columns `Geneid`, `GO`, `KEGG`, for optional GO/KEGG enrichment of temporal gene clusters.
- The actual eggNOG-mapper annotation available for *C. melo* (`04_Cm_Annotations_eggNOG-mapper.xlsx`, 21,649 rows) is an **Excel** file with columns `Name` (not `Geneid`), `GOs` (not `GO`), `KEGG_ko` (not `KEGG`), plus the same remaining columns as the *P. xanthii* annotation above.
- **As shipped, this enrichment step will not run against the real annotation file** - format and column names both differ. See "Repository status" in the main README; this has been left unmodified pending your review (it is not clear whether this code path is still meant to run at all, since the article's Methods states *C. melo* functional enrichment was performed with DAVID, not clusterProfiler).

## Metabolomics (`data/metabolomics/`)

### `002_METABO_NEG_AVG.txt`
- **Omics layer**: untargeted LC-MS/MS metabolomics, negative ionization mode (ESI-).
- **Format**: comma-separated despite the `.txt` extension; first column holds feature identifiers (`mz_<value>`), first data row is a `Label` sentinel row (schema check only, not used further), remaining rows are features.
- **Rows**: 3,901 features (matches the 3,901 ESI- features reported in the article's Technical Validation).
- **Columns**: 72 samples, `<time>h_<condition>_R<replicate>`.
- **Values**: peak intensities, averaged across the two technical replicates per biological sample (per the "AVG" filename and the article's Methods).
- **Used by**: `scripts/metabolomics/01_metabolomics_workflow.R` (dataset tag `NEG`).

### `002_METABO_POS_AVG.txt`
- Same format as the NEG file, positive ionization mode (ESI+).
- **Rows**: 8,215 features (matches the article's 8,215 ESI+ features).
- **Used by**: `scripts/metabolomics/01_metabolomics_workflow.R` (dataset tag `POS`).

### `003_METABO_COMB_AVG.txt`
- Same format, combining NEG and POS features.
- **Rows**: 12,116 features (= 3,901 + 8,215, i.e. the union of both ionization modes).
- **Used by**: `scripts/metabolomics/01_metabolomics_workflow.R` (dataset tag `COMB`).

### `14_Metabolomics_Sets.xlsx`
- **Omics layer**: MarVis-Pathway pathway/set enrichment output.
- **Rows**: 67 pathways/sets.
- **Columns**: `Index`, `Set ID`, `Set name`, `Score`, `Size`, `Marker hits`, `Entry hits`, `Cpd hits`, `Db/cpd/intern hits`.
- **Used by**: `scripts/metabolomics/02_pathway_enrichment.R` (KEGG-prefixed sets only, top 16 by score/hits/coverage).

### `15_Metabolomics_Significant.xlsx`
- **Omics layer**: significant differentially accumulated metabolite features (MarVis-Filter ANOVA output) and their putative metabolite annotations.
- **Sheets**:
  - `Significant`: 352 rows (matches the article's 352 significant metabolic features). Columns: `ID`, `rt`, `m/z (corrected)`, 72 sample intensity columns, `dataset`, `ANOVApValue_FDR`, `ARule`, `nC13`, `CosSum`, `FormerY`, `nC`, `MarVisPathwayEntries`.
  - `Putative`: 89 rows (matches the article's 89 features with putative metabolite annotations); not read by any script in this repository.
- **Values**: processed/normalized intensities per sample, plus per-feature ANOVA FDR-adjusted p-value and annotation metadata.
- **Used by**: `scripts/metabolomics/03_dynamic_metabolites.R` (reads the `Significant` sheet only).
