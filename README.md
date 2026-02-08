# scATAC-seq Analysis Pipeline

A comprehensive R-based workflow for single-cell ATAC-seq (scATAC-seq) data analysis using Signac and Seurat frameworks.

## Overview

This pipeline provides a complete workflow for analyzing single-cell chromatin accessibility data, from quality control to differential accessibility analysis and biological interpretation. The workflow is designed for 10X Genomics scATAC-seq data and includes integration with scRNA-seq for cell type annotation.

## Features

- ✅ **Quality Control**: Comprehensive QC metrics including TSS enrichment, nucleosome signal, and blacklist ratio
- ✅ **Data Preprocessing**: TF-IDF normalization and LSI dimensionality reduction
- ✅ **Cell Clustering**: Graph-based clustering and UMAP visualization
- ✅ **Gene Activity Quantification**: Chromatin accessibility-based gene activity scoring
- ✅ **Multi-omics Integration**: Integration with scRNA-seq data for cell type annotation
- ✅ **Differential Accessibility**: Statistical testing for differentially accessible peaks
- ✅ **Functional Annotation**: GO enrichment analysis of regulatory regions
- ✅ **Genomic Visualization**: Coverage plots and region-specific accessibility profiles

## Requirements

### R Version
- R >= 4.0.0

### Required Packages

```r
# Bioconductor packages
BiocManager::install(c('GenomeInfoDb', 'SeuratObject', 'BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86', 'BSgenome.Mmusculus.UCSC.mm10', 'EnsDb.Mmusculus.v79'))

# CRAN packages
install.packages(c('Seurat', 'Signac', 'hdf5r', 'ggplot2', 'patchwork', 'dplyr'))

# GitHub packages
remotes::install_github('stuart-lab/signac', ref = 'develop')
remotes::install_github('immunogenomics/presto')

# Additional packages for enrichment analysis
BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db', 'enrichplot'))
```

## Data Requirements

This pipeline expects 10X Genomics scATAC-seq output files:

- `*_filtered_peak_bc_matrix.h5` - Peak-barcode matrix in HDF5 format
- `*_singlecell.csv` - Cell-level metadata
- `*_fragments.tsv.gz` - Fragment file with index (`.tbi`)
- (Optional) scRNA-seq reference data for cell type annotation (`.rds` format)

### Example Data

The pipeline uses 10k PBMC scATAC-seq data from 10X Genomics as demonstration:
- Dataset: 10k Human PBMCs, ATAC v2 chemistry
- Reference genome: GRCh38 (hg38)
- Annotation: Ensembl v98

### Run Step-by-Step

The script is organized into functional sections:

#### Load Data and Create Seurat Object
```r
# Load peak matrix, metadata, and fragments
counts <- Read10X_h5(filename = "filtered_peak_bc_matrix.h5")
metadata <- read.csv("singlecell.csv", header = TRUE, row.names = 1)
chrom_assay <- CreateChromatinAssay(counts = counts, ...)
pbmc <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", ...)
```

#### Quality Control
```r
# Compute QC metrics
pbmc <- NucleosomeSignal(object = pbmc)
pbmc <- TSSEnrichment(object = pbmc)
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

# Filter cells
pbmc <- subset(pbmc, 
  subset = nCount_peaks > 9000 & 
           nCount_peaks < 100000 &
           pct_reads_in_peaks > 40 &
           blacklist_ratio < 0.01 &
           nucleosome_signal < 4 &
           TSS.enrichment > 4
)
```

#### Normalization and Dimensionality Reduction
```r
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
```

#### Clustering
```r
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
```

#### Differential Accessibility Analysis
```r
# Wilcox test (fast)
da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "CD4 Naive",
  ident.2 = "CD14+ Monocytes",
  test.use = 'wilcox',
  min.pct = 0.1
)

# Likelihood ratio test (recommended, more rigorous)
lr_da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "CD4 Naive",
  ident.2 = "CD14+ Monocytes",
  test.use = 'LR',
  latent.vars = 'nCount_peaks',
  min.pct = 0.1
)
```

## Workflow Overview

```
Raw Data (10X)
    ↓
Quality Control
    ↓
Normalization (TF-IDF)
    ↓
Dimensionality Reduction (LSI, UMAP)
    ↓
Clustering & Visualization
    ↓
Gene Activity Quantification
    ↓
scRNA-seq Integration (Cell Type Annotation)
    ↓
Differential Accessibility Analysis
    ↓
GO Enrichment Analysis
    ↓
Genomic Region Visualization
```

## Key Analysis Sections

### 1. Quality Control Metrics

- **TSS Enrichment**: Signal enrichment at transcription start sites (cutoff: > 4)
- **Nucleosome Signal**: Nucleosome banding pattern strength (cutoff: < 4)
- **Peak Region Fragments**: Percentage of reads in peaks (cutoff: > 40%)
- **Blacklist Ratio**: Fraction of reads in genomic blacklist regions (cutoff: < 0.01)

### 2. Cell Type Annotation

Integration with scRNA-seq reference data using canonical correlation analysis (CCA):
- Transfer cell type labels from RNA to ATAC
- Validate predictions using marker gene activity
- Filter low-confidence predictions

### 3. Differential Accessibility Testing

Two statistical methods provided:
- **Wilcoxon test**: Fast, suitable for exploratory analysis
- **Likelihood Ratio (LR) test**: Recommended for publication, accounts for sequencing depth

### 4. Visualization Options

- Violin plots for QC metrics
- UMAP projections for cell clusters
- Coverage plots for genomic regions
- Feature plots for gene activity/peak accessibility
- Heatmaps and dot plots for marker peaks

## Output Files

The script creates a `R.results/` directory containing:
- Processed Seurat object
- QC plots
- UMAP visualizations
- Differential accessibility results
- GO enrichment analysis results
- Coverage plots

## Example Visualizations

### Quality Control
- Violin plots showing distribution of QC metrics
- TSS enrichment profile plots
- Fragment length histograms
- Density scatter plots (peaks vs TSS enrichment)

### Cell Type Analysis
- UMAP plots colored by cluster or predicted cell type
- Gene activity feature plots for marker genes
- Prediction score distributions

### Differential Accessibility
- Volcano plots or MA plots (can be added)
- Coverage plots highlighting differential peaks
- Dot plots showing top marker peaks per cell type

### Genomic Regions
- Coverage browser for interactive exploration
- Tile plots showing individual cell accessibility
- Multi-track plots with gene annotations
