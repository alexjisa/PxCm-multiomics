# PxCm-multiomics

Codigo de analisis asociado al dataset y articulo:

> **Temporal multi-omics dataset of a cucurbit-powdery mildew interaction**
> Jimenez-Sanchez A., Fernandez-Ortuno D., Pastor V., Polonio A., Perez-Garcia A.
> *Scientific Data* (en preparacion / 2026).

Recursos temporales de RNA-seq dual (planta-hongo) y metabolomica LC-MS/MS
untargeted de la interaccion *Cucumis melo* - *Podosphaera xanthii* (oidio del
melon), muestreados en 12 tiempos entre 12 y 144 horas post-inoculacion.

Este repositorio contiene unicamente el codigo. Los datos (crudos y
procesados) estan depositados en repositorios publicos, ver mas abajo.

## Disponibilidad de datos

- **RNA-seq crudo (FASTQ)**: NCBI SRA, BioProject
  [PRJNA1454108](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1454108).
- **LC-MS/MS crudo (metabolomica)**: MetaboLights, accession
  [MTBLS15447](https://www.ebi.ac.uk/metabolights/MTBLS15447).
- **Datos procesados** (matrices de conteo/TPM, resultados de expresion
  diferencial, enriquecimiento funcional, tablas de metabolomica y metadatos
  completos): Zenodo,
  [DOI 10.5281/zenodo.19555639](https://doi.org/10.5281/zenodo.19555639).

## Estructura del repositorio

```
scripts/
  transcriptomics/
    Cm_transcriptomics.R          # C. melo: TPM, PCA/t-SNE/UMAP, DESeq2, clustering temporal, GO/KEGG
    Px_transcriptomics.R          # P. xanthii: TPM, PCA/t-SNE/UMAP, DESeq2, clustering temporal, KEGG
    detected_genes_violin.R       # Genes detectados (TPM>=1) por grupo: Cm Mock/Inoculated y Px
  metabolomics/
    pathway_enrichment.R          # Cobertura de rutas KEGG (salida de MarVis-Pathway)
    top_dynamic_metabolites_heatmap.R  # Heatmap de las 20 caracteristicas mas dinamicas
```

Cada script asume que se ejecuta con el directorio de trabajo en la carpeta
donde estan los archivos de datos procesados correspondientes (descargados
desde Zenodo, ver arriba). Archivos de entrada esperados por script:

| Script | Archivos de entrada |
|---|---|
| `Cm_transcriptomics.R` | `01_Transcriptomics_Cm_counts.txt`, `02_Transcriptomics_Cm_sample_info.txt` |
| `Px_transcriptomics.R` | `01_Transcriptomics_Px_counts.txt`, `02_Transcriptomics_Px_sample_info.txt`, `Px_functional_annotation.xlsx` (opcional, anotacion eggNOG-mapper para enriquecimiento KEGG) |
| `detected_genes_violin.R` | `03_Transcriptomics_Cm_TPM.txt`, `03_Transcriptomics_Px_TPM.txt` |
| `pathway_enrichment.R` | `14_Metabolomics_Sets.xlsx` |
| `top_dynamic_metabolites_heatmap.R` | `15_Metabolomics_Significant.xlsx` |

`Cm_transcriptomics.R` busca los archivos en `data/processed/` si existe esa
carpeta bajo el directorio de trabajo, y si no, en el directorio de trabajo
directamente. Los demas scripts los buscan en el directorio de trabajo.

## Metodologia (resumen)

- **Transcriptomica**: FastQC (v0.11.9) + MultiQC (v1.27) para control de
  calidad; Trimmomatic (v0.39) para recorte de adaptadores/calidad; HISAT2
  (v2.2.1) para alineamiento contra el genoma de *C. melo* v4.0 (muestras
  mock) y, para muestras inoculadas, alineamiento secuencial contra el genoma
  de *C. melo* y despues el de *P. xanthii* aislado 2086; SAMtools (v1.23.1)
  para el procesado de alineamientos; featureCounts (v2.0.6) para el conteo
  de lecturas; DESeq2 (v1.46.0) para expresion diferencial (padj < 0.05,
  |log2FC| > 1); enriquecimiento funcional con DAVID (v6.8, KEGG, *C. melo*)
  y clusterProfiler (v4.14.6, GO via eggNOG-mapper v2, *P. xanthii*).
- **Metabolomica**: adquisicion LC-MS/MS (ESI+/ESI-) con MassLynx (v4.2);
  procesado con el paquete `xcms` (deteccion de picos centWave, correccion de
  tiempo de retencion, agrupacion, imputacion); filtrado e identificacion
  tentativa de metabolitos con MarVis-Suite 2.0 (MarVis-Filter,
  MarVis-Pathway); caracteristicas diferenciales por ANOVA de un factor con
  correccion FDR (Benjamini-Hochberg, p.adj < 0.05).
- Analisis realizado en R (v4.4.2).

Detalles completos de parametros en la seccion Methods del articulo.

## Dependencias (R)

`DESeq2`, `ggplot2`, `dplyr`, `tidyr`, `stringr`, `readr`, `readxl`,
`pheatmap`, `viridis`, `ggrepel`, `Rtsne`, `uwot`, `dendextend`, `scales`,
`svglite`. Opcionales (enriquecimiento funcional): `clusterProfiler`,
`enrichplot`, `KEGGREST`, `httr`, `xml2`, `rvest`.

## Estado del repositorio

Este repositorio esta en construccion. Actualmente incluye el procesado
DESeq2/clustering/enriquecimiento de transcriptomica y dos scripts de
visualizacion de metabolomica (cobertura de rutas y heatmap de dinamica
temporal). Pendiente de anadir el script de control de calidad y analisis
exploratorio de metabolomica (PCA de muestras, boxplot/densidad de
intensidades, correlacion y CV entre replicas, valores ausentes, volcano
plots por punto temporal y el resumen de resultados temporales), que
actualmente solo existe como figuras/tablas de salida sin el codigo fuente
incorporado aqui todavia.

## Citacion

Si usas este codigo o el dataset asociado, por favor cita el articulo (ver
cabecera de este README) y los repositorios de datos correspondientes
(BioProject PRJNA1454108, MetaboLights MTBLS15447, Zenodo DOI
10.5281/zenodo.19555639).
