setRepositories(ind=1:3) 
install.packages("Signac")
install.packages("hdf5r")
BiocManager::install(version = "3.22")
BiocManager::install("GenomeInfoDb", force = T)
BiocManager::install("SeuratObject", force = T)
install.packages("remotes") 
remotes::install_github("stuart-lab/signac", ref = "develop")
remotes::install_github('immunogenomics/presto') 

BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86'))
BiocManager::install(c('BSgenome.Mmusculus.UCSC.mm10', 'EnsDb.Mmusculus.v79'))

# Set work place
setwd("/Users/lixia/Data/data/scATAC-seq/PBMCs.scATACseq")
# Creat reult directory
dir.create("R.results", showWarnings = FALSE)

# Load library
library(Seurat)
library(Signac)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(GenomeInfoDb)

############################################### load data and samples #######################################################
# Creat a Seurat object using the peak/cell matrix and cell metadata, and store the path to the fragment file
# Load peaks
counts <- Read10X_h5(filename = "10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5")
# View counts
dim(counts)
counts[1:5, 1:5]
length(counts@x) / (nrow(counts) * ncol(counts)) * 100 #the proportion of non-zero values

# Load metadata
metadata <- read.csv(
  file = "10k_pbmc_ATACv2_nextgem_Chromium_Controller_singlecell.csv",
  header = TRUE,
  row.names = 1
)

# View metadata
head(metadata)
colnames(metadata)

# Load fragments
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

# View chrom_assay
chrom_assay
head(chrom_assay@counts[1:5, 1:5])
granges(chrom_assay)

# Create Seurat bject
pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

# View 
Assays(pbmc)
DefaultAssay(pbmc)
pbmc

###############################################  gene annotations #######################################################

# Extract the ChromatinAssay object named peaks from the Seurat object named pbmc.
pbmc[['peaks']] # chr1:9772-10660 format

# Convert the peak coordinates to the standard GRanges format in R.
granges(pbmc)

table(seqnames(granges(pbmc))) # Count the number of peaks on each chromosome
mean(width(granges(pbmc))) # Calculate the average length of the peak

#Remove peaks on scaffolds
peaks.keep <- seqnames(granges(pbmc)) %in% standardChromosomes(granges(pbmc))
pbmc <- pbmc[as.vector(peaks.keep), ]
pbmc

#---------------------------------------------------------------------------------#
# Download genome annotations
library(AnnotationHub)
ah <- AnnotationHub(cache = "/Users/lixia/Data/database/AnnotationHub") # Store in this path.

# Locate the annotation version corresponding to the reference genome，
# and store the annotation information in ensdb_v98 (EnsDb object).
# the reference package 10x Genomics used to perform the mapping was “GRCh38-2020-A”,
# which corresponds to the Ensembl v98 patch release.
query(ah, "EnsDb.Hsapiens.v98") 
ensdb_v98 <- ah[["AH75011"]] #EnsDb object
ensdb_v98
head(genes(ensdb_v98))

# Convert EnsDb object into GRanges object
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)

# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))#Change the chromosome name from 1 to chr1
genome(annotations) <- "hg38"

seqlevels(annotations)
genome(annotations)

#---------------------------------------------------------------------------------#
# add the gene information to the object
Annotation(pbmc) <- annotations

pbmc[['peaks']]
ClosestFeature(pbmc, regions = "chr1-100-500") # the name of the gene closest to this coordinate

###############################################  Computing QC Metrics #######################################################

# Computing
pbmc <- NucleosomeSignal(object = pbmc) # compute nucleosome signal score per cell
head(pbmc$nucleosome_signal)

pbmc <- TSSEnrichment(object = pbmc) # compute TSS enrichment score per cell
head(pbmc$TSS.enrichment)

pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100 # add fraction of reads in peaks

pbmc$blacklist_ratio <- FractionCountsInRegion(  # add blacklist ratio
  object = pbmc, 
  assay = 'peaks',
  regions = blacklist_hg38_unified
)
#---------------------------------------------------------------------------------#
# Visualization : Visualiz the all QCs
VlnPlot(
  object = pbmc,
  features = c('nCount_peaks','nucleosome_signal', 'TSS.enrichment', 'pct_reads_in_peaks', 'blacklist_ratio'),
  pt.size = 0.1,
  ncol = 5
)

# Visualiz the nucleosome_signal
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')

# Visualiz the TSS.enrichment
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE) #fast = FALSE，save matrix
pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
p <- TSSPlot(pbmc, group.by = 'high.tss') + 
  ggtitle("TSS Enrichment Score Profile") +
  theme(plot.title = element_text(hjust = 0.5))
print(p)

# Visualiz two QCs
DensityScatter(pbmc, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
DensityScatter(pbmc, x = 'TSS.enrichment', y = 'nucleosome_signal', log_x = TRUE)

#---------------------------------------------------------------------------------#
# Subsetting：remove cells that are outliers for these QC metrics. 
pbmc <- subset(
  x = pbmc,
  subset = nCount_peaks > 9000 &
    nCount_peaks < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 4 &
    TSS.enrichment > 4
)

pbmc

###############################################  Normalization and linear dimensional reduction #######################################################

# (TF-IDF) normalization:sequencing depth,giving higher values to more rare peaks.
pbmc <- RunTFIDF(pbmc)

# high frequency peaks: high signal/noise ratio, cell identity
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')

# latent semantic indexing (LSI)：SVD on the TD-IDF matrix, using the feature peaks
pbmc <- RunSVD(pbmc)
# pbmc <- RunSVD(object = pbmc, n = 50)
# 默认情况下，RunSVD() 会计算前 50 个维度（Component）

# LSI_1 is indeced by sequence depth，and  remove LSI_1 for clustering
# Visualize the correlation between each LSI component and sequencing depth
DepthCor(pbmc)

# 修改 ndims 参数来显示前 50 个维度
ElbowPlot(object = pbmc, reduction = "lsi", ndims = 50)

###############################################  Non-linear dimension reduction and clustering #######################################################

#graph-based clustering and non-linear dimension reduction for visualization， which are same with scRNA-seq
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)

# pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3, resolution = 0.5)

DimPlot(object = pbmc, label = TRUE) + NoLegend()

############################################### Create a gene activity matrix #######################################################

# Quantify the activity of each gene in the genome by assessing the chromatin accessibility
gene.activities <- GeneActivity(pbmc)

# Add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)

pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)

DefaultAssay(pbmc) <- 'RNA'
#---------------------------------------------------------------------------------#
# Visualize the activities of canonical marker genes 
FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

############################################### Integrating with scRNA-seq data #######################################################

# Load the pre-processed scRNA-seq data for PBMCs
pbmc_rna <- readRDS("pbmc_10k_v3.rds")
pbmc_rna <- UpdateSeuratObject(pbmc_rna)

# CCA reduction and find anchors
transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc,
  reduction = 'cca'
)

# Predict atac cell types
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
)

#Write the cell type to pbmc
pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)
head(pbmc@meta.data)
#---------------------------------------------------------------------------------#

# Visualize the celltype annotation of RNA and ATAC
plot1 <- DimPlot(
  object = pbmc_rna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = pbmc,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2

# Plotting the prediction score for the cells assigned to each label 
VlnPlot(pbmc, 'prediction.score.max', group.by = 'predicted.id')

#---------------------------------------------------------------------------------#
# Removw cell annotations with <=20 cells total
predicted_id_counts <- table(pbmc$predicted.id)
major_predicted_ids <- names(predicted_id_counts[predicted_id_counts > 20])
pbmc <- pbmc[, pbmc$predicted.id %in% major_predicted_ids]

# change cell identities to the per-cell predicted labels
Idents(pbmc) <- pbmc$predicted.id
table(Idents(pbmc))

#---------------------------------------------------------------------------------#
# name cluster according to the most common predicted label for that cluste
# replace each cluster label with its most likely predicted label
for(i in levels(pbmc)) {
  cells_to_reid <- WhichCells(pbmc, idents = i)
  newid <- names(which.max(table(pbmc$predicted.id[cells_to_reid])))
  Idents(pbmc, cells = cells_to_reid) <- newid
}

DimPlot(object = pbmc, label = TRUE) + NoLegend()

############################################### Find differentially accessible peaks between cell types #####################################
#Comparing any groups of cells can be compared using these methods

# change back to working with peaks instead of gene activities
DefaultAssay(pbmc) <- 'peaks'

#---------------------------------------------------------------------------------#
# wilcox is the default option for test.use
library(presto)

da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "CD4 Naive",
  ident.2 = "CD14+ Monocytes",
  test.use = 'wilcox',           # Presto will automatically accelerate.
  min.pct = 0.1                  # same with scRNA
)

head(da_peaks)


plot1 <- VlnPlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("CD4 Naive","CD14+ Monocytes")
)

plot2 <- FeaturePlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)

plot1 | plot2

#---------------------------------------------------------------------------------#
# LR (Likelihood Ratio Test) test is recommended for scATAC， the binary matrix. 
# latent.vars = 'nCount_peaks', mitigates the effect of differential sequencing depth.
# More computation time is needed。
# The most rigorous and recommended way to write for difference analysis

lr_da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "CD4 Naive",
  ident.2 = "CD14+ Monocytes",
  test.use = 'LR',
  latent.vars = 'nCount_peaks',
  min.pct = 0.1
)

head(lr_da_peaks)

plot3 <- VlnPlot(
  object = pbmc,
  features = rownames(lr_da_peaks)[1],
  pt.size = 0.1,
  idents = c("CD4 Naive","CD14+ Monocytes")
)

plot4 <- FeaturePlot(
  object = pbmc,
  features = rownames(lr_da_peaks)[1],
  pt.size = 0.1
)

plot3 | plot4

#---------------------------------------------------------------------------------#
# Find all marker

all_da_peaks <- FindAllMarkers(
  object = pbmc,
  min.pct = 0.1,             
  test.use = 'wilcox',         
  latent.vars = 'nCount_peaks' 
)

# Extract the two most prominent peaks for each cell type.
library(dplyr)
top2_peaks <- all_da_peaks %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# Visualize
DotPlot(pbmc, features = top2_peaks$gene) + RotatedAxis()

#---------------------------------------------------------------------------------#
# Find the closest gene to peaks
open_cd4naive <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
open_cd14mono <- rownames(da_peaks[da_peaks$avg_log2FC < -3, ])

closest_genes_cd4naive <- ClosestFeature(pbmc, regions = open_cd4naive)
closest_genes_cd14mono <- ClosestFeature(pbmc, regions = open_cd14mono)

head(closest_genes_cd4naive)
head(closest_genes_cd14mono)

# Visualize peak and its surrounding genes
region_to_plot <- "CD4" # gene ID or coordinate string

p <- CoveragePlot(
  object = pbmc,
  region = region_to_plot,
  features = region_to_plot,    # 在下方显示该基因的表达量或活性（可选）
  annotation = TRUE,            # 显示基因结构轨道
  peaks = TRUE,                 # 显示已有的 Peak 召唤（Calling）轨道
  idents = c("CD4 Naive", "CD14+ Monocytes"), # 只对比这两类细胞
  extend.upstream = 1000,       # 向左延伸 1kb 视角
  extend.downstream = 1000      # 向右延伸 1kb 视角
)

p

#---------------------------------------------------------------------------------#

#GO enrichment analysis with clusterProfiler
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

cd4naive_ego <- enrichGO(gene = closest_genes_cd4naive$gene_id,
                         keyType = "ENSEMBL",
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05,
                         readable = TRUE)

barplot(cd4naive_ego,showCategory = 20)
dotplot(cd4naive_ego, showCategory = 20)
cnetplot(cd4naive_ego)

cd14mono_ego <- enrichGO(gene = closest_genes_cd14mono$gene_id,
                         keyType = "ENSEMBL",
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05,
                         readable = TRUE)

barplot(cd14mono_ego,showCategory = 20)
dotplot(cd14mono_ego, showCategory = 20)
cnetplot(cd14mono_ego)

############################################### Plotting genomic regions ############################################### 
# Visualize pseudo-bulk peak of each cell types in a specific region
# Rank by cell type name
pbmc <- SortIdents(pbmc)

# find DA peaks overlapping gene of interest
regions_highlight <- subsetByOverlaps(StringToGRanges(open_cd4naive), LookupGeneCoords(pbmc, "CD4"))

CoveragePlot(
  object = pbmc,
  region = "CD4",
  region.highlight = regions_highlight,
  extend.upstream = 1000,
  extend.downstream = 1000
)

regions_highlight <- subsetByOverlaps(StringToGRanges(open_cd14mono), LookupGeneCoords(pbmc, "LYZ"))

CoveragePlot(
  object = pbmc,
  region = "LYZ",
  region.highlight = regions_highlight,
  extend.upstream = 1000,
  extend.downstream = 5000
)

#---------------------------------------------------------------------------------#
# Interactive version of CoveragePlot
CoverageBrowser(pbmc, region = "CD4")

#---------------------------------------------------------------------------------#
# Visualize peak of each cell in a specific region
region_to_plot <- "CD4" # gene ID or coordinate string

p <- TilePlot(
  object = pbmc,
  region = region_to_plot,
  # idents = c("CD4 Naive", "CD14+ Monocytes"), # 对比两组
  tile.size = 100,                             # 每个方块代表 100bp
  extend.upstream = 1000,
  extend.downstream = 1000
)

p

