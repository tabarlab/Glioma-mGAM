library(data.table)
library(dplyr)
library(Signac)
library(Seurat)
library(viridis)
library(BuenColors)
library(EnsDb.Hsapiens.v86)

counts <- Read10X_h5(filename = "../data/patient3/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "../data/patient3/singlecell.csv",
  header = TRUE,
  row.names = 1
)

# make mito mutation matrix
af1 <- readRDS("../output/d3_dn_afs.rds"); colnames(af1) <- substr(start = 1, stop = 18, gsub("-1", "-1", colnames(af1)))
af2 <- readRDS("../output/d3_dp_afs.rds"); colnames(af2) <- substr(start = 1, stop = 18, gsub("-1", "-2", colnames(af2)))
af3 <- readRDS("../output/d3_mono_afs.rds"); colnames(af3) <- substr(start = 1, stop = 18, gsub("-1", "-3", colnames(af3)))
mito_mat <- data.matrix(cbind(cbind(af1, af2), af3))

cells <- intersect(colnames(counts), colnames(mito_mat))
metadata2 <- metadata[cells,]

chrom_assay <- CreateChromatinAssay(
  counts = counts[grepl("^chr", rownames(counts)),cells],
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = '../data/patient3/fragments.tsv.gz',
  min.cells = 0,
  min.features = 0
)
dim(chrom_assay)
so <- CreateSeuratObject(
  counts = chrom_assay,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata2
)

so$pct_reads_in_peaks <- so$peak_region_fragments / so$passed_filters * 100

so_filt <- subset(
  x = so,
  peak_region_fragments < 100000 &
    pct_reads_in_peaks > 25 &
    passed_filters > 10^3
)
dim(so_filt)
dim(so)

# Annotate samples
idx <- substr(colnames(so_filt), 18, 18)
so_filt$channel <- case_when(
  idx == "1" ~ "DN", 
  idx == "2" ~ "DP", 
  idx == "3" ~ "Mono"
)

# Do dimension reduction
so_filt <- RunTFIDF(so_filt)
so_filt <- FindTopFeatures(so_filt, min.cutoff = 'q50')
so_filt <- RunSVD(
  object = so_filt,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

so_filt <- RunUMAP(object = so_filt, reduction = 'lsi', dims = 2:30)
so_filt <- FindNeighbors(object = so_filt, reduction = 'lsi', dims = 2:30)
so_filt <- FindClusters(object = so_filt, verbose = FALSE, algorithm = 3, resolution = 0.5)

DimPlot(so_filt, group.by = c("channel", "seurat_clusters"), label = TRUE)

# Compute gene scores
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(so_filt) <- annotations
gene.activities <- GeneActivity(so_filt)

# Scale gene score data
so_filt[['RNA']] <- CreateAssayObject(counts = gene.activities)
so_filt <- NormalizeData(
  object = so_filt,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(so_filt$nCount_RNA)
)

DefaultAssay(so_filt) <- 'RNA'

# Add gene sets for module scoring
df <- fread("../data/longboi_genesets.tsv", header = FALSE) %>%
  data.frame() 
gene_sets <- unique(df$V2)
lapply(gene_sets, function(gs){
  df %>% dplyr::filter(V2 == gs) %>% pull(V1)
}) -> list_genes
names(list_genes) <- gene_sets
so_filt <- AddModuleScore(so_filt, list_genes, name = gene_sets)
FeaturePlot(so_filt, c("FOSL2.small1", "FOSL2.large2", "pos.discrim3", "neg.discrim4"), min.cutoff = "q05", max.cutoff = "q95") &
  scale_color_viridis()

VlnPlot(
  so_filt,
  c("FOSL2.small1", "FOSL2.large2", "pos.discrim3", "neg.discrim4"),
  group.by = "channel", pt.size = 0)

fm <- FindMarkers(so_filt, group.by = "channel", ident.1 = "DP", ident.2 = "DN")
head(fm, 20)

FeaturePlot(so_filt, c("SHOX", "TERT", "C8A"), min.cutoff = "q05", max.cutoff = "q95", sort = TRUE) &
  scale_color_viridis()

so_filt[['mito']] <- CreateAssayObject(counts = mito_mat[,colnames(so_filt)])
DefaultAssay(so_filt) <- 'mito'

FeaturePlot(so_filt, sort(rowSums(mito_mat > .5)) %>% tail(3) %>% names(), sort.cell = TRUE, split.by = "channel") &
  scale_color_viridis()

