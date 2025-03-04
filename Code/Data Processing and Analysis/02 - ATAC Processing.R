library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(ggrepel)
set.seed(1234)


#use futures for parallel processing to decrease time

library(future)
plan()
plan("multicore", workers = 10)
plan()
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM


Filtering_TSS_Nuc <- function(object, annotations){
    plot_list<-list()  

    # add the gene information to the object
    Annotation(object) <- annotations

    # compute nucleosome signal score per cell
    object <- NucleosomeSignal(object = object)
    object <- TSSEnrichment(object = object, fast = FALSE)
    # add blacklist ratio and fraction of reads in peaks
    object$pct_reads_in_peaks <- object$peak_region_fragments / object$passed_filters * 100
    object$blacklist_ratio <- object$blacklist_region_fragments / object$peak_region_fragments
    object$high.tss <- ifelse(object$TSS.enrichment > 2, 'High', 'Low')
    plot_list[["TSS"]]<-TSSPlot(object, group.by = 'high.tss') + NoLegend()

    object$nucleosome_group <- ifelse(object$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
    plot_list[["Fragments"]]<-FragmentHistogram(object = object, group.by = 'nucleosome_group')
    plot_list[["Summary_plot"]]<-VlnPlot(
    object = object,
    features = c('pct_reads_in_peaks', 'peak_region_fragments',
                'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
    pt.size = 0.1,
    ncol = 5
    )
    plot_list[["object"]]<-object

    return(plot_list)
}

# Combining the MSK101 original and re-do datasets template
counts <- Read10X_h5(filename = "./ATAC/MSK101_ATAC/MSK101_ATAC_filtered_peak_bc_matrix.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = './ATAC/MSK101_ATAC/MSK101_ATAC_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)
ATAC_MSK101 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)

# Defining annotations variable 
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

meta_ATAC_MSK101<-read.csv("./ATAC/MSK101_ATAC_singlecell.csv",header = T,row.names = 1)


ATAC_MSK101<-AddMetaData(ATAC_MSK101,meta_ATAC_MSK101)
filter_list_MSK101<-Filtering_TSS_Nuc(ATAC_MSK101,annotations=annotations) 
ATAC_MSK101<-filter_list_MSK101$object

ATAC_MSK101 <- subset(
  x = ATAC_MSK101,
  subset = peak_region_fragments > 100&
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

counts <- Read10X_h5(filename = "./ATAC/MSK101_ATAC/second_round/MSK101_Redo_ATAC_filtered_peak_bc_matrix.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = './ATAC/MSK101_ATAC/second_round/MSK101_Redo_ATAC_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)
ATAC_MSK101_redo <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)
meta_ATAC_MSK101_redo<-read.csv("./ATAC/MSK101_ATAC/second_round/MSK101_Redo_ATAC_singlecell.csv",header = T,row.names = 1)

ATAC_MSK101_redo<-AddMetaData(ATAC_MSK101_redo,meta_ATAC_MSK101_redo)
filter_list_MSK101_redo<-Filtering_TSS_Nuc(object=ATAC_MSK101_redo,annotations = annotations) 
ATAC_MSK101_redo<-filter_list_MSK101_redo$object


ATAC_MSK101_redo <- subset(
  x = ATAC_MSK101_redo,
  subset = peak_region_fragments > 100 &
    peak_region_fragments < 60000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)


rangesMSK101<-reduce(c(StringToGRanges(row.names(ATAC_MSK101)),StringToGRanges(row.names(ATAC_MSK101_redo))))

# Filter out bad peaks based on length
peakwidths <- width(rangesMSK101)
combined.peaks <- rangesMSK101[peakwidths  < 10000 & peakwidths > 20]
#combined.peaks

feature_MSK101 <- FeatureMatrix(
  fragments = Fragments(ATAC_MSK101),
  features = combined.peaks,
  cells = colnames(ATAC_MSK101)
)

feature_MSK101_redo<- FeatureMatrix(
  fragments = Fragments(ATAC_MSK101_redo),
  features = combined.peaks,
  cells = colnames(ATAC_MSK101_redo)
)

ATAC_MSK101[["merged"]]<-CreateChromatinAssay(feature_MSK101)

ATAC_MSK101_redo[["merged"]]<-CreateChromatinAssay(feature_MSK101_redo)


# add information to identify dataset of origin
ATAC_MSK101$dataset <- 'Original'
ATAC_MSK101_redo$dataset <- 'Redo'

DefaultAssay(ATAC_MSK101)<-"merged"

DefaultAssay(ATAC_MSK101_redo)<-"merged"
# merge all datasets, adding a cell ID to make sure cell names are unique
combined_MSK101 <- merge(
  x = ATAC_MSK101,
  y = list(ATAC_MSK101_redo),
  add.cell.ids = c("orig", "redo")
)

combined_MSK101 <- RunTFIDF(combined_MSK101)
combined_MSK101 <- FindTopFeatures(combined_MSK101, min.cutoff = 'q0')
combined_MSK101 <- RunSVD(combined_MSK101)

combined_MSK101 <- RunUMAP(object = combined_MSK101, reduction = 'lsi', dims = 2:30)
combined_MSK101 <- FindNeighbors(object = combined_MSK101, reduction = 'lsi', dims = 2:30)
combined_MSK101 <- FindClusters(object = combined_MSK101, verbose = FALSE, algorithm = 3)
DimPlot(object = combined_MSK101, label = TRUE,group.by = "dataset")

#add back fragments

DefaultAssay(combined_MSK101)<-"peaks"
fragments<-Fragments(combined_MSK101)
DefaultAssay(combined_MSK101)<-"merged"
Fragments(combined_MSK101)<-fragments

# Combining the MSK102 original and re-do datasets

counts <- Read10X_h5(filename = "./ATAC/MSK102_ATAC/filtered_peak_bc_matrix.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = './ATAC/MSK102_ATAC/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)
ATAC_MSK102 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)

meta_ATAC_MSK102<-read.csv("./ATAC/MSK102_ATAC_singlecell.csv",header = T,row.names = 1)


ATAC_MSK102<-AddMetaData(ATAC_MSK102,meta_ATAC_MSK102)
filter_list_MSK102<-Filtering_TSS_Nuc(ATAC_MSK102,annotations=annotations) 
ATAC_MSK102<-filter_list_MSK102$object

# adjust these values based on graphs obtained in step above
ATAC_MSK102 <- subset(
  x = ATAC_MSK102,
  subset = peak_region_fragments > 100 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

counts <- Read10X_h5(filename = "./ATAC/MSK102_ATAC/second_round/MSK102_Redo_ATAC_filtered_peak_bc_matrix.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = './ATAC/MSK102_ATAC/second_round/MSK102_Redo_ATAC_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)
ATAC_MSK102_redo <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)

meta_ATAC_MSK102_redo<-read.csv("./ATAC/MSK102_ATAC/second_round/MSK102_Redo_ATAC_singlecell.csv",header = T,row.names = 1)

ATAC_MSK102_redo<-AddMetaData(ATAC_MSK102_redo,meta_ATAC_MSK102_redo)
filter_list_MSK102_redo<-Filtering_TSS_Nuc(object=ATAC_MSK102_redo,annotations = annotations) 
ATAC_MSK102_redo<-filter_list_MSK102_redo$object


ATAC_MSK102_redo <- subset(
  x = ATAC_MSK102_redo,
  subset = peak_region_fragments > 100 &
    peak_region_fragments < 60000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

rangesMSK102<-reduce(c(StringToGRanges(row.names(ATAC_MSK102)),StringToGRanges(row.names(ATAC_MSK102_redo))))

# Filter out bad peaks based on length
peakwidths <- width(rangesMSK102)
combined.peaks <- rangesMSK102[peakwidths  < 10000 & peakwidths > 20]

feature_MSK102 <- FeatureMatrix(
  fragments = Fragments(ATAC_MSK102),
  features = combined.peaks,
  cells = colnames(ATAC_MSK102)
)

feature_MSK102_redo<- FeatureMatrix(
  fragments = Fragments(ATAC_MSK102_redo),
  features = combined.peaks,
  cells = colnames(ATAC_MSK102_redo)
)

ATAC_MSK102[["merged"]]<-CreateChromatinAssay(feature_MSK102)

ATAC_MSK102_redo[["merged"]]<-CreateChromatinAssay(feature_MSK102_redo)


# add information to identify dataset of origin
ATAC_MSK102$dataset <- 'Original'
ATAC_MSK102_redo$dataset <- 'Redo'

DefaultAssay(ATAC_MSK102)<-"merged"

DefaultAssay(ATAC_MSK102_redo)<-"merged"
# merge all datasets, adding a cell ID to make sure cell names are unique
combined_MSK102 <- merge(
  x = ATAC_MSK102,
  y = list(ATAC_MSK102_redo),
  add.cell.ids = c("orig", "redo")
)

combined_MSK102 <- RunTFIDF(combined_MSK102)
combined_MSK102 <- FindTopFeatures(combined_MSK102, min.cutoff = 'q0')
combined_MSK102 <- RunSVD(combined_MSK102)

combined_MSK102 <- RunUMAP(object = combined_MSK102, reduction = 'lsi', dims = 2:30)
combined_MSK102 <- FindNeighbors(object = combined_MSK102, reduction = 'lsi', dims = 2:30)
combined_MSK102 <- FindClusters(object = combined_MSK102, verbose = FALSE, algorithm = 3)
DimPlot(object = combined_MSK102, label = TRUE,group.by = "dataset")

DefaultAssay(combined_MSK102)<-"peaks"
fragments<-Fragments(combined_MSK102)
DefaultAssay(combined_MSK102)<-"merged"
Fragments(combined_MSK102)<-fragments



# Combining the MSK108 original and re-do datasets

counts <- Read10X_h5(filename = "./ATAC/MSK108_ATAC/filtered_peak_bc_matrix.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = './ATAC/MSK108_ATAC/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)
ATAC_MSK108 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)

meta_ATAC_MSK108<-read.csv("./ATAC/MSK108_ATAC_singlecell.csv",header = T,row.names = 1)


ATAC_MSK108<-AddMetaData(ATAC_MSK108,meta_ATAC_MSK108)
filter_list_MSK108<-Filtering_TSS_Nuc(ATAC_MSK108,annotations=annotations) 
ATAC_MSK108<-filter_list_MSK108$object

# adjust these values based on graphs obtained in step above
ATAC_MSK108 <- subset(
  x = ATAC_MSK108,
  subset = peak_region_fragments > 100 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

counts <- Read10X_h5(filename = "./ATAC/MSK108_ATAC/second_round/MSK108_Redo_ATAC_filtered_peak_bc_matrix.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = './ATAC/MSK108_ATAC/second_round/MSK108_Redo_ATAC_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)
ATAC_MSK108_redo <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)

meta_ATAC_MSK108_redo<-read.csv("./ATAC/MSK108_ATAC/second_round/MSK108_Redo_ATAC_singlecell.csv",header = T,row.names = 1)

ATAC_MSK108_redo<-AddMetaData(ATAC_MSK108_redo,meta_ATAC_MSK108_redo)
filter_list_MSK108_redo<-Filtering_TSS_Nuc(object=ATAC_MSK108_redo,annotations = annotations) 
ATAC_MSK108_redo<-filter_list_MSK108_redo$object


ATAC_MSK108_redo <- subset(
  x = ATAC_MSK108_redo,
  subset = peak_region_fragments > 100 &
    peak_region_fragments < 60000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

rangesMSK108<-reduce(c(StringToGRanges(row.names(ATAC_MSK108)),StringToGRanges(row.names(ATAC_MSK108_redo))))

# Filter out bad peaks based on length
peakwidths <- width(rangesMSK108)
combined.peaks <- rangesMSK108[peakwidths  < 10000 & peakwidths > 20]

feature_MSK108 <- FeatureMatrix(
  fragments = Fragments(ATAC_MSK108),
  features = combined.peaks,
  cells = colnames(ATAC_MSK108)
)

feature_MSK108_redo<- FeatureMatrix(
  fragments = Fragments(ATAC_MSK108_redo),
  features = combined.peaks,
  cells = colnames(ATAC_MSK108_redo)
)

ATAC_MSK108[["merged"]]<-CreateChromatinAssay(feature_MSK108)

ATAC_MSK108_redo[["merged"]]<-CreateChromatinAssay(feature_MSK108_redo)


# add information to identify dataset of origin
ATAC_MSK108$dataset <- 'Original'
ATAC_MSK108_redo$dataset <- 'Redo'

DefaultAssay(ATAC_MSK108)<-"merged"

DefaultAssay(ATAC_MSK108_redo)<-"merged"
# merge all datasets, adding a cell ID to make sure cell names are unique
combined_MSK108 <- merge(
  x = ATAC_MSK108,
  y = list(ATAC_MSK108_redo),
  add.cell.ids = c("orig", "redo")
)

combined_MSK108 <- RunTFIDF(combined_MSK108)
combined_MSK108 <- FindTopFeatures(combined_MSK108, min.cutoff = 'q0')
combined_MSK108 <- RunSVD(combined_MSK108)

combined_MSK108 <- RunUMAP(object = combined_MSK108, reduction = 'lsi', dims = 2:30)
combined_MSK108 <- FindNeighbors(object = combined_MSK108, reduction = 'lsi', dims = 2:30)
combined_MSK108 <- FindClusters(object = combined_MSK108, verbose = FALSE, algorithm = 3)
DimPlot(object = combined_MSK108, label = TRUE,group.by = "dataset")

DefaultAssay(combined_MSK108)<-"peaks"
fragments<-Fragments(combined_MSK108)
DefaultAssay(combined_MSK108)<-"merged"
Fragments(combined_MSK108)<-fragments



# Combining the MSK109 original and re-do datasets

counts <- Read10X_h5(filename = "./ATAC/MSK109-PM_ATAC/filtered_peak_bc_matrix.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = './ATAC/MSK109-PM_ATAC/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)
ATAC_MSK109 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)

# Check the file path if correct
meta_ATAC_MSK109<-read.csv("./ATAC/MSK109_PM_ATAC_singlecell.csv",header = T,row.names = 1)


ATAC_MSK109<-AddMetaData(ATAC_MSK109,meta_ATAC_MSK109)
filter_list_MSK109<-Filtering_TSS_Nuc(ATAC_MSK109,annotations=annotations) 
ATAC_MSK109<-filter_list_MSK109$object

# adjust these values based on graphs obtained in step above
ATAC_MSK109 <- subset(
  x = ATAC_MSK109,
  subset = peak_region_fragments > 100 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

counts <- Read10X_h5(filename = "./ATAC/MSK109-PM_Redo_ATAC_filtered_peak_bc_matrix.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = './ATAC/MSK109-PM_Redo_ATAC_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)
ATAC_MSK109_redo <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)

meta_ATAC_MSK109_redo<-read.csv("./ATAC/MSK109-PM_Redo_ATAC_singlecell.csv",header = T,row.names = 1)

ATAC_MSK109_redo<-AddMetaData(ATAC_MSK109_redo,meta_ATAC_MSK109_redo)
filter_list_MSK109_redo<-Filtering_TSS_Nuc(object=ATAC_MSK109_redo,annotations = annotations) 
ATAC_MSK109_redo<-filter_list_MSK109_redo$object


ATAC_MSK109_redo <- subset(
  x = ATAC_MSK109_redo,
  subset = peak_region_fragments > 100 &
    peak_region_fragments < 60000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

rangesMSK109<-reduce(c(StringToGRanges(row.names(ATAC_MSK109)),StringToGRanges(row.names(ATAC_MSK109_redo))))

# Filter out bad peaks based on length
peakwidths <- width(rangesMSK109)
combined.peaks <- rangesMSK109[peakwidths  < 10000 & peakwidths > 20]

feature_MSK109 <- FeatureMatrix(
  fragments = Fragments(ATAC_MSK109),
  features = combined.peaks,
  cells = colnames(ATAC_MSK109)
)

feature_MSK109_redo<- FeatureMatrix(
  fragments = Fragments(ATAC_MSK109_redo),
  features = combined.peaks,
  cells = colnames(ATAC_MSK109_redo)
)

ATAC_MSK109[["merged"]]<-CreateChromatinAssay(feature_MSK109)

ATAC_MSK109_redo[["merged"]]<-CreateChromatinAssay(feature_MSK109_redo)


# add information to identify dataset of origin
ATAC_MSK109$dataset <- 'Original'
ATAC_MSK109_redo$dataset <- 'Redo'

DefaultAssay(ATAC_MSK109)<-"merged"

DefaultAssay(ATAC_MSK109_redo)<-"merged"
# merge all datasets, adding a cell ID to make sure cell names are unique
combined_MSK109 <- merge(
  x = ATAC_MSK109,
  y = list(ATAC_MSK109_redo),
  add.cell.ids = c("orig", "redo")
)

combined_MSK109 <- RunTFIDF(combined_MSK109)
combined_MSK109 <- FindTopFeatures(combined_MSK109, min.cutoff = 'q0')
combined_MSK109 <- RunSVD(combined_MSK109)

combined_MSK109 <- RunUMAP(object = combined_MSK109, reduction = 'lsi', dims = 2:30)
combined_MSK109 <- FindNeighbors(object = combined_MSK109, reduction = 'lsi', dims = 2:30)
combined_MSK109 <- FindClusters(object = combined_MSK109, verbose = FALSE, algorithm = 3)
DimPlot(object = combined_MSK109, label = TRUE,group.by = "dataset")

DefaultAssay(combined_MSK109)<-"peaks"
fragments<-Fragments(combined_MSK109)
DefaultAssay(combined_MSK109)<-"merged"
Fragments(combined_MSK109)<-fragments


# add metadata to combined batch 1 datasets
combined_MSK101$batch = "1"
combined_MSK101$type = "WT"
combined_MSK101$patient = "MSK101"
combined_MSK101$histo = "gbm"
combined_MSK101$grade = "4"

combined_MSK102$batch = "1"
combined_MSK102$type = "WT"
combined_MSK102$patient = "MSK102"
combined_MSK102$histo = "gbm"
combined_MSK102$grade = "4"

combined_MSK108$batch = "1"
combined_MSK108$type = "MUT"
combined_MSK108$patient = "MSK108"
combined_MSK108$histo = "astro"
combined_MSK108$grade = "3"

combined_MSK109$batch = "1"
combined_MSK109$type = "MUT"
combined_MSK109$patient = "117pm"
combined_MSK109$histo = "astro"
combined_MSK109$grade = "2"



# Batch 2
# MSK115 ATAC script

counts <- Read10X_h5(filename = "./ATAC/MSK115-AM_ATAC/filtered_peak_bc_matrix.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = './ATAC/MSK115-AM_ATAC/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)



ATAC_MSK115 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)
ATAC_MSK115$patient<-"MSK115"
ATAC_MSK115$type<- "MUT"
ATAC_MSK115$batch<- "2"
ATAC_MSK115$histo<- "oligo"
ATAC_MSK115$grade<- "3"


meta_ATAC_MSK115<-read.csv("./ATAC/metadata/singlecell_MSK115.csv",header = T,row.names = 1)

ATAC_MSK115<-AddMetaData(ATAC_MSK115,meta_ATAC_MSK115)

filter_list_MSK115 <- Filtering_TSS_Nuc(ATAC_MSK115, annotations=annotations) 
ATAC_MSK115 <- filter_list_MSK115$object

# Adjust cut-off values based on graphs generated above
ATAC_MSK115 <- subset(
  x = ATAC_MSK115,
  subset = peak_region_fragments > 100&
    peak_region_fragments < 40000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)


# MSK114 ATAC script

counts <- Read10X_h5(filename = "./ATAC/MSK114_ATAC/filtered_peak_bc_matrix.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = './ATAC/MSK114_ATAC/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

ATAC_MSK114 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)
ATAC_MSK114$patient<-"MSK114"
ATAC_MSK114$type<- "MUT"
ATAC_MSK114$batch<- "2"
ATAC_MSK114$histo<- "astro"
ATAC_MSK114$grade<- "4"


meta_ATAC_119am<-read.csv("./ATAC/metadata/singlecell-19.csv",header = T,row.names = 1)

ATAC_MSK114<-AddMetaData(ATAC_MSK114,meta_ATAC_119am)

filter_list_MSK114 <- Filtering_TSS_Nuc(ATAC_MSK114, annotations=annotations) 
ATAC_MSK114 <- filter_list_MSK114$object

# Adjust cut-off values based on graphs generated above
ATAC_MSK114 <- subset(
  x = ATAC_MSK114,
  subset = peak_region_fragments > 100&
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

# MSK110 ATAC script

counts <- Read10X_h5(filename = "./ATAC/MSK110_ATAC/filtered_peak_bc_matrix.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = './ATAC/MSK110_ATAC/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

ATAC_MSK110 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)
ATAC_MSK110$patient<-"MSK110"
ATAC_MSK110$type<- "MUT"
ATAC_MSK110$batch<-"2"
ATAC_MSK110$histo <- "oligo"
ATAC_MSK110$grade <- "3"


meta_ATAC_MSK110<-read.csv("./ATAC/metadata/singlecell_MSK110.csv",header = T,row.names = 1)

ATAC_MSK110<-AddMetaData(ATAC_MSK110,meta_ATAC_MSK110)


filter_list_MSK110 <- Filtering_TSS_Nuc(ATAC_MSK110, annotations=annotations) 
ATAC_MSK110 <- filter_list_MSK110$object

# Adjust cut-off values based on graphs generated above
ATAC_MSK110 <- subset(
  x = ATAC_MSK110,
  subset = peak_region_fragments > 100 &
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

# MSK111 ATAC script

counts <- Read10X_h5(filename = "./ATAC/MSK111_ATAC/filtered_peak_bc_matrix.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = './ATAC/MSK111_ATAC/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

ATAC_MSK111 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)
ATAC_MSK111$patient<-"MSK111"
ATAC_MSK111$type<- "MUT"
ATAC_MSK111$batch<-"2"
ATAC_MSK111$histo<-"astro"
ATAC_MSK111$grade<-"3"

meta_ATAC_MSK111<-read.csv("./ATAC/metadata/singlecell_MSK111.csv",header = T,row.names = 1)

ATAC_MSK111<-AddMetaData(ATAC_MSK111,meta_ATAC_MSK111)

filter_list_MSK111 <- Filtering_TSS_Nuc(ATAC_MSK111, annotations=annotations) 
ATAC_MSK111 <- filter_list_MSK111$object

# Adjust cut-off values based on graphs generated above
ATAC_MSK111 <- subset(
  x = ATAC_MSK111,
  subset = peak_region_fragments > 100 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

# MSK112 ATAC script

counts <- Read10X_h5(filename = "./ATAC/MSK112_ATAC/filtered_peak_bc_matrix.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = './ATAC/MSK112_ATAC/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

ATAC_MSK112 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)
ATAC_MSK112$patient<-"MSK112"
ATAC_MSK112$type<- "MUT"
ATAC_MSK112$batch<-"2"
ATAC_MSK112$histo<-"astro"
ATAC_MSK112$grade<-"3"

meta_ATAC_MSK112<-read.csv("./ATAC/metadata/singlecell_MSK112.csv",header = T,row.names = 1)

ATAC_MSK112<-AddMetaData(ATAC_MSK112,meta_ATAC_MSK112)

filter_list_MSK112 <- Filtering_TSS_Nuc(ATAC_MSK112, annotations=annotations) 
ATAC_MSK112 <- filter_list_MSK112$object

# Adjust cut-off values based on graphs generated above
ATAC_MSK112 <- subset(
  x = ATAC_MSK112,
  subset = peak_region_fragments > 100&
    peak_region_fragments < 30000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

# MSK103 ATAC script

counts <- Read10X_h5(filename = "./ATAC/MSK103_ATAC/filtered_peak_bc_matrix.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = './ATAC/MSK103_ATAC/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

ATAC_MSK103 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)
ATAC_MSK103$patient<-"MSK103"
ATAC_MSK103$type<- "WT"
ATAC_MSK103$batch<-"2"
ATAC_MSK103$histo<-"gbm"
ATAC_MSK103$grade<-"4"


meta_ATAC_MSK103<-read.csv("./ATAC/metadata/singlecell_MSK103.csv",header = T,row.names = 1)

ATAC_MSK103<-AddMetaData(ATAC_MSK103,meta_ATAC_MSK103)

filter_list_MSK103 <- Filtering_TSS_Nuc(ATAC_MSK103, annotations=annotations) 
ATAC_MSK103 <- filter_list_MSK103$object

# Adjust cut-off values based on graphs generated above
ATAC_MSK103 <- subset(
  x = ATAC_MSK103,
  subset = peak_region_fragments > 100 &
    peak_region_fragments < 15000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

# MSK113 ATAC script

counts <- Read10X_h5(filename = "./ATAC/MSK113_ATAC/filtered_peak_bc_matrix.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = './ATAC/MSK113_ATAC/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)


ATAC_MSK113 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)
ATAC_MSK113$patient<-"MSK113"
ATAC_MSK113$type<- "MUT"
ATAC_MSK113$batch<- "2"
ATAC_MSK113$histo<- "astro"
ATAC_MSK113$grade<- "4"



meta_ATAC_MSK113<-read.csv("./ATAC/metadata/singlecell_MSK113.csv",header = T,row.names = 1)

ATAC_MSK113<-AddMetaData(ATAC_MSK113,meta_ATAC_MSK113)

filter_list_MSK113 <- Filtering_TSS_Nuc(ATAC_MSK113, annotations=annotations) 
ATAC_MSK113 <- filter_list_MSK113$object

# Adjust cut-off values based on graphs generated above
ATAC_MSK113 <- subset(
  x = ATAC_MSK113,
  subset = peak_region_fragments > 100&
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

# Batch 3

# MSK107 ATAC script

counts <- Read10X_h5(filename = "./ATAC/MSK107_ATAC/filtered_peak_bc_matrix_MSK107.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = './ATAC/MSK107_ATAC/fragments_MSK107.tsv.gz',
  min.cells = 10,
  min.features = 200
)



ATAC_MSK107 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)

ATAC_MSK107$patient<-"MSK107"
ATAC_MSK107$type<- "WT"
ATAC_MSK107$batch<- "3"
ATAC_MSK107$histo<- "gbm"
ATAC_MSK107$grade<- "4"


meta_ATAC_MSK107<-read.csv("./ATAC/MSK107_ATAC/singlecell_MSK107.csv",header = T,row.names = 1)

ATAC_MSK107<-AddMetaData(ATAC_MSK107,meta_ATAC_MSK107)

filter_list_MSK107 <- Filtering_TSS_Nuc(ATAC_MSK107, annotations=annotations) 
ATAC_MSK107 <- filter_list_MSK107$object

# Adjust cut-off values based on graphs generated above
ATAC_MSK107 <- subset(
  x = ATAC_MSK107,
  subset = peak_region_fragments > 100&
    peak_region_fragments < 40000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

# MSK105 ATAC script

counts <- Read10X_h5(filename = "./ATAC/MSK105_ATAC/filtered_peak_bc_matrix_MSK105.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = './ATAC/MSK105_ATAC/fragments_MSK105.tsv.gz',
  min.cells = 10,
  min.features = 200
)



ATAC_MSK105 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)


ATAC_MSK105$patient<-"MSK105"
ATAC_MSK105$type<- "WT"
ATAC_MSK105$batch<- "3"
ATAC_MSK105$histo<- "gbm"
ATAC_MSK105$grade<- "4"



meta_ATAC_MSK105<-read.csv("./ATAC/MSK105_ATAC/singlecell_MSK105.csv",header = T,row.names = 1)

ATAC_MSK105<-AddMetaData(ATAC_MSK105,meta_ATAC_MSK105)

filter_list_MSK105 <- Filtering_TSS_Nuc(ATAC_MSK105, annotations=annotations) 
ATAC_MSK105 <- filter_list_MSK105$object

# Adjust cut-off values based on graphs generated above
ATAC_MSK105 <- subset(
  x = ATAC_MSK105,
  subset = peak_region_fragments > 100&
    peak_region_fragments < 40000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

# MSK106 ATAC script

counts <- Read10X_h5(filename = "./ATAC/MSK106_ATAC/filtered_peak_bc_matrix_MSK106.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = './ATAC/MSK106_ATAC/fragments_MSK106.tsv.gz',
  min.cells = 10,
  min.features = 200
)



ATAC_MSK106 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)

ATAC_MSK106$patient<-"MSK106"
ATAC_MSK106$type<- "WT"
ATAC_MSK106$batch<- "3"
ATAC_MSK106$histo<- "gbm"
ATAC_MSK106$grade<- "4"


meta_ATAC_MSK106<-read.csv("./ATAC/MSK106_ATAC/singlecell_MSK106.csv",header = T,row.names = 1)

ATAC_MSK106<-AddMetaData(ATAC_MSK106,meta_ATAC_MSK106)

filter_list_MSK106 <- Filtering_TSS_Nuc(ATAC_MSK106, annotations=annotations) 
ATAC_MSK106 <- filter_list_MSK106$object

# Adjust cut-off values based on graphs generated above
ATAC_MSK106 <- subset(
  x = ATAC_MSK106,
  subset = peak_region_fragments > 100&
    peak_region_fragments < 40000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

# MSK104 ATAC script

counts <- Read10X_h5(filename = "./ATAC/MSK104_ATAC/filtered_peak_bc_matrix_MSK104.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = './ATAC/MSK104_ATAC/fragments_MSK104.tsv.gz',
  min.cells = 10,
  min.features = 200
)


ATAC_MSK104 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)

ATAC_MSK104$patient<-"MSK104"
ATAC_MSK104$type<- "WT"
ATAC_MSK104$batch<- "3"
ATAC_MSK104$histo<- "gbm"
ATAC_MSK104$grade<- "4"

meta_ATAC_MSK104<-read.csv("./ATAC/MSK104_ATAC/singlecell_MSK104.csv",header = T,row.names = 1)

ATAC_MSK104<-AddMetaData(ATAC_MSK104,meta_ATAC_MSK104)

filter_list_MSK104 <- Filtering_TSS_Nuc(ATAC_MSK104, annotations=annotations) 
ATAC_MSK104 <- filter_list_MSK104$object

# Adjust cut-off values based on graphs generated above
ATAC_MSK104 <- subset(
  x = ATAC_MSK104,
  subset = peak_region_fragments > 100&
    peak_region_fragments < 40000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

################### Merging samples ###################

# convert to genomic ranges
gr.MSK109 <- StringToGRanges(row.names(combined_MSK109))
gr.MSK108 <- StringToGRanges(row.names(combined_MSK108))
gr.MSK101 <- StringToGRanges(row.names(combined_MSK101))
gr.MSK102 <- StringToGRanges(row.names(combined_MSK102))
gr.MSK115 <- StringToGRanges(row.names(ATAC_MSK115))
gr.MSK114 <- StringToGRanges(row.names(ATAC_MSK114))
gr.MSK110 <- StringToGRanges(row.names(ATAC_MSK110))
gr.MSK111 <- StringToGRanges(row.names(ATAC_MSK111))
gr.MSK112 <- StringToGRanges(row.names(ATAC_MSK112))
gr.MSK103 <- StringToGRanges(row.names(ATAC_MSK103))
gr.MSK113 <- StringToGRanges(row.names(ATAC_MSK113))
gr.MSK104 <- StringToGRanges(row.names(ATAC_MSK104))
gr.MSK106 <- StringToGRanges(row.names(ATAC_MSK106))
gr.MSK107 <- StringToGRanges(row.names(ATAC_MSK107))
gr.MSK105 <- StringToGRanges(row.names(ATAC_MSK105))
gr.nb1 <- StringToGRanges(row.names(ATAC_NB1))
gr.nb2 <- StringToGRanges(row.names(ATAC_NB2))

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.MSK109, gr.MSK108, gr.MSK101, gr.MSK102, gr.MSK115, gr.MSK114, gr.MSK110, gr.MSK111, gr.MSK112, gr.MSK103, gr.MSK113,gr.MSK104,gr.MSK106,gr.MSK107,gr.MSK105, gr.nb1, gr.nb2))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# adjust core usage for computationally heavy steps
plan("multiprocess", workers = 8)
options(future.globals.maxSize = 50 * 1024 ^ 3)

feature_MSK101 <-FeatureMatrix(
  fragments = Fragments(combined_MSK101[['peaks']]),
  features = combined.peaks,
  cells = colnames(combined_MSK101[['peaks']])
)

combined_MSK101[["merge.peaks"]]<-CreateChromatinAssay(feature_MSK101, fragments = Fragments(combined_MSK101[['merged']]))
DefaultAssay(combined_MSK101)<-"merge.peaks"

feature_MSK102 <-FeatureMatrix(
  fragments = Fragments(combined_MSK102[['peaks']]),
  features = combined.peaks,
  cells = colnames(combined_MSK102[['peaks']])
)

combined_MSK102[["merge.peaks"]]<-CreateChromatinAssay(feature_MSK102, fragments = Fragments(combined_MSK102[['merged']]))
DefaultAssay(combined_MSK102)<-"merge.peaks"

feature_MSK108 <-FeatureMatrix(
  fragments = Fragments(combined_MSK108[['peaks']]),
  features = combined.peaks,
  cells = colnames(combined_MSK108[['peaks']])
)

combined_MSK108[["merge.peaks"]]<-CreateChromatinAssay(feature_MSK108, fragments = Fragments(combined_MSK108[['peaks']]))
DefaultAssay(combined_MSK108)<-"merge.peaks"

feature_MSK109 <-FeatureMatrix(
  fragments = Fragments(combined_MSK109[['peaks']]),
  features = combined.peaks,
  cells = colnames(combined_MSK109[['peaks']])
)

combined_MSK109[["merge.peaks"]]<-CreateChromatinAssay(feature_MSK109, fragments = Fragments(combined_MSK109[['merged']]))
DefaultAssay(combined_MSK109)<-"merge.peaks"

feature_MSK115 <-FeatureMatrix(
  fragments = Fragments(ATAC_MSK115),
  features = combined.peaks,
  cells = colnames(ATAC_MSK115)
)

ATAC_MSK115[["merge.peaks"]]<-CreateChromatinAssay(feature_MSK115, fragments = Fragments(ATAC_MSK115))
DefaultAssay(ATAC_MSK115)<-"merge.peaks"

feature_MSK114 <-FeatureMatrix(
  fragments = Fragments(ATAC_MSK114),
  features = combined.peaks,
  cells = colnames(ATAC_MSK114)
)

ATAC_MSK114[["merge.peaks"]]<-CreateChromatinAssay(feature_MSK114, fragments = Fragments(ATAC_MSK114))
DefaultAssay(ATAC_MSK114)<-"merge.peaks"

feature_MSK110 <-FeatureMatrix(
  fragments = Fragments(ATAC_MSK110),
  features = combined.peaks,
  cells = colnames(ATAC_MSK110)
)

ATAC_MSK110[["merge.peaks"]]<-CreateChromatinAssay(feature_MSK110, fragments = Fragments(ATAC_MSK110))
DefaultAssay(ATAC_MSK110)<-"merge.peaks"

feature_MSK111 <-FeatureMatrix(
  fragments = Fragments(ATAC_MSK111),
  features = combined.peaks,
  cells = colnames(ATAC_MSK111)
)

ATAC_MSK111[["merge.peaks"]]<-CreateChromatinAssay(feature_MSK111, fragments = Fragments(ATAC_MSK111))
DefaultAssay(ATAC_MSK111)<-"merge.peaks"

feature_MSK112 <-FeatureMatrix(
  fragments = Fragments(ATAC_MSK112),
  features = combined.peaks,
  cells = colnames(ATAC_MSK112)
)

ATAC_MSK112[["merge.peaks"]]<-CreateChromatinAssay(feature_MSK112, fragments = Fragments(ATAC_MSK112))
DefaultAssay(ATAC_MSK112)<-"merge.peaks"

feature_MSK103 <-FeatureMatrix(
  fragments = Fragments(ATAC_MSK103),
  features = combined.peaks,
  cells = colnames(ATAC_MSK103)
)

ATAC_MSK103[["merge.peaks"]]<-CreateChromatinAssay(feature_MSK103, fragments = Fragments(ATAC_MSK103))
DefaultAssay(ATAC_MSK103)<-"merge.peaks"

feature_MSK113 <-FeatureMatrix(
  fragments = Fragments(ATAC_MSK113),
  features = combined.peaks,
  cells = colnames(ATAC_MSK113)
)

ATAC_MSK113[["merge.peaks"]]<-CreateChromatinAssay(feature_MSK113, fragments = Fragments(ATAC_MSK113))
DefaultAssay(ATAC_MSK113)<-"merge.peaks"

feature_MSK107 <-FeatureMatrix(
  fragments = Fragments(ATAC_MSK107),
  features = combined.peaks,
  cells = colnames(ATAC_MSK107)
)

ATAC_MSK107[["merge.peaks"]]<-CreateChromatinAssay(feature_MSK107, fragments = Fragments(ATAC_MSK107))
DefaultAssay(ATAC_MSK107)<-"merge.peaks"

feature_MSK104 <-FeatureMatrix(
  fragments = Fragments(ATAC_MSK104),
  features = combined.peaks,
  cells = colnames(ATAC_MSK104)
)

ATAC_MSK104[["merge.peaks"]]<-CreateChromatinAssay(feature_MSK104, fragments = Fragments(ATAC_MSK104))
DefaultAssay(ATAC_MSK104)<-"merge.peaks"

feature_MSK106 <-FeatureMatrix(
  fragments = Fragments(ATAC_MSK106),
  features = combined.peaks,
  cells = colnames(ATAC_MSK106)
)

ATAC_MSK106[["merge.peaks"]]<-CreateChromatinAssay(feature_MSK106, fragments = Fragments(ATAC_MSK106))
DefaultAssay(ATAC_MSK106)<-"merge.peaks"

feature_MSK105 <-FeatureMatrix(
  fragments = Fragments(ATAC_MSK105),
  features = combined.peaks,
  cells = colnames(ATAC_MSK105)
)

ATAC_MSK105[["merge.peaks"]]<-CreateChromatinAssay(feature_MSK105, fragments = Fragments(ATAC_MSK105))
DefaultAssay(ATAC_MSK105)<-"merge.peaks"


unintegrated_b1_b2[["merged"]]<-NULL
unintegrated_b1_b2[["peaks"]]<-NULL


unintegrated_b1_b2 <-merge(combined_MSK101, c(combined_MSK102,combined_MSK108, combined_MSK109, ATAC_MSK115, ATAC_MSK114,ATAC_MSK110, ATAC_MSK111, ATAC_MSK112, ATAC_MSK103, ATAC_MSK113,ATAC_MSK104,ATAC_MSK106,ATAC_MSK107,ATAC_MSK105), add.cell.ids = c("combined_MSK101", "combined_MSK102","combined_MSK108", "combined_MSK109","ATAC_MSK115", "ATAC_MSK114","ATAC_MSK110", "ATAC_MSK111", "ATAC_MSK112", "ATAC_MSK103", "ATAC_MSK113","ATAC_MSK104","ATAC_MSK106","ATAC_MSK107","ATAC_MSK105"))
unintegrated_b1_b2 <-merge(unintegrated_b1_b2, c(ATAC_MSK104), add.cell.ids = c("","ATAC_MSK104"))
unintegrated_b1_b2 <-merge(unintegrated_b1_b2, c(ATAC_MSK106), add.cell.ids = c("","ATAC_MSK106"))
unintegrated_b1_b2 <-merge(unintegrated_b1_b2, c(ATAC_MSK107), add.cell.ids = c("","ATAC_MSK107"))
unintegrated_b1_b2 <-merge(unintegrated_b1_b2, c(ATAC_MSK105), add.cell.ids = c("","ATAC_MSK105"))
unintegrated_b1_b2 <-merge(unintegrated_b1_b2, c(ATAC_NB1), add.cell.ids = c("","ATAC_NB1"))
unintegrated_b1_b2 <-merge(unintegrated_b1_b2, c(ATAC_NB2), add.cell.ids = c("","ATAC_NB2"))


DefaultAssay(unintegrated_b1_b2) <- "merge.peaks"

#filtering step added in based on nucleosome signal

# compute nucleosome signal score per cell
unintegrated_b1_b2 <- NucleosomeSignal(object = unintegrated_b1_b2)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(unintegrated_b1_b2) <- annotations

unintegrated_b1_b2 <- RunTFIDF(unintegrated_b1_b2)
unintegrated_b1_b2 <- FindTopFeatures(unintegrated_b1_b2, min.cutoff = 'q50')
unintegrated_b1_b2 <- RunSVD(unintegrated_b1_b2,n=100)

unintegrated_b1_b2 <- RunUMAP(object = unintegrated_b1_b2, reduction = 'lsi', dims = 2:100)
unintegrated_b1_b2 <- FindNeighbors(object = unintegrated_b1_b2, reduction = 'lsi', dims = 2:100)

unintegrated_b1_b2 <- FindClusters(object = unintegrated_b1_b2, verbose = FALSE, algorithm = 3)
pa1<-DimPlot(object = unintegrated_b1_b2, label = TRUE)
pa2<-DimPlot(object = unintegrated_b1_b2, label = TRUE, group.by = "type")
pa3<-DimPlot(object = unintegrated_b1_b2, label = TRUE, group.by = "batch")
pa4<-DimPlot(object = unintegrated_b1_b2, label = TRUE, group.by = "histo")
pa5<-DimPlot(object = unintegrated_b1_b2, label = TRUE, group.by = "grade")

pa6<-DimPlot(object = unintegrated_b1_b2, label = TRUE, group.by = "patient")

unint_list2<-list(pa1,pa2,pa3,pa4,pa5,pa6)
wrap_plots(unint_list2, nrow = 2,ncol = 3)


saveRDS(unintegrated_b1_b2,file="./Unintegrated_all_ATAC.rds")

#we ID'd the ES cells previously so remove those
old <- readRDS("Coembed_RNA_ATAC_all_unint.rds")
unintegrated_b1_b2_ATAC<-readRDS(file="./ATAC_all_unint.rds")
Remove_these_cells<-Idents(unintegrated_b1_b2_ATAC)[Idents(unintegrated_b1_b2_ATAC)=="Brain/ES"]
remove(unintegrated_b1_b2_ATAC)


rename_df<-data.frame(orig_names=row.names(unintegrated_b1_b2[[]]))

rename_df$changed_name<-gsub(x = rename_df$orig_names,pattern = "^_+",replacement = "",perl = T)

rename_df$isES<-rename_df$changed_name %in% names(Remove_these_cells)

row.names(rename_df)<-rename_df$orig_names

sum(rename_df$isES)

unintegrated_b1_b2<-AddMetaData(unintegrated_b1_b2,rename_df)

unintegrated_b1_b2_ATAC<-subset(unintegrated_b1_b2,subset= isES == FALSE)

unintegrated_b1_b2_ATAC<-RenameCells(unintegrated_b1_b2_ATAC, new.names = unintegrated_b1_b2_ATAC$changed_name)

saveRDS(unintegrated_b1_b2_ATAC,file="./Unintegrated_all_ATAC.rds")


# Run Harmony on the unintegrated.test object

library(harmony)
library(stringr)
hm.integrated <- RunHarmony(
  object = unintegrated_b1_b2,
  group.by.vars = 'batch',
  reduction = 'lsi',
  assay.use = 'merge.peaks',
  project.dim = FALSE
)

hm.integrated <- RunUMAP(hm.integrated, dims = 2:30, reduction = 'harmony')
pa1<-DimPlot(object = hm.integrated, label = TRUE)
pa2<-DimPlot(object = hm.integrated, label = TRUE, group.by = "type")
pa3<-DimPlot(object = hm.integrated, label = TRUE, group.by = "batch")
pa4<-DimPlot(object = hm.integrated, label = TRUE, group.by = "histo")
pa5<-DimPlot(object = hm.integrated, label = TRUE, group.by = "grade")

pa6<-DimPlot(object = hm.integrated, label = TRUE, group.by = "patient")

int_list2<-list(pa1,pa2,pa3,pa4,pa5,pa6)
wrap_plots(int_list2, nrow = 2,ncol = 3)

saveRDS(hm.integrated,file="./all_ATAC.rds")

################### ChromVar Analysis ###################

library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)

int_b1_b2_ATAC<-readRDS(file="./all_ATAC.rds")

ATAC_GAM<-subset(int_b1_b2_ATAC,subset= celltype_K  %in% c("GAMS", "Microglia_Mono"))
DefaultAssay(ATAC_GAM) <- 'merge.peaks'

ATAC_GAM <- FindTopFeatures(ATAC_GAM, min.cutoff = "q90" )
ATAC_GAM <- RunTFIDF(ATAC_GAM)
ATAC_GAM <- RunSVD(ATAC_GAM)
ATAC_GAM <- RunUMAP(ATAC_GAM, reduction = 'lsi', dims = 2:30, reduction.name = 'umap.GAM')
DimPlot(ATAC_GAM,group.by = "type")

#analysis of DA markers
Idents(ATAC_GAM)<-ATAC_GAM$type
GAM_WT_vs_MUT_markers<-FindMarkers(ATAC_GAM,only.pos = F,min.pct = 0,logfc.threshold = 0.05,ident.1 = "WT",ident.2 = "MUT")
GAM_WT_vs_MUT_markers2 <- presto:::wilcoxauc.Seurat(X = ATAC_GAM, group_by = 'type',assay = 'data', seurat_assay = 'merge.peaks',min.pct = 0)
GAM_WT_vs_MUT_markers<-GAM_WT_vs_MUT_markers2[GAM_WT_vs_MUT_markers2$group=="MUT",]
CoveragePlot(ATAC_GAM,region = StringToGRanges("chrX-73850594-73852041"),group.by = "patient")
CoveragePlot(ATAC_GAM,region = StringToGRanges("chrX-131740926-131742410"),group.by = "patient")
CoveragePlot(ATAC_GAM,region = StringToGRanges("chrX-131740926-131742410"),group.by = "patient")
CoveragePlot(ATAC_GAM,region = StringToGRanges("chr7-121519672-121520827"),group.by = "type")

ClosestFeature(object = ATAC_GAM,regions = StringToGRanges("chr7-121519672-121520827"))
row.names(GAM_WT_vs_MUT_markers)<-GAM_WT_vs_MUT_markers$feature
GAM_WT_vs_MUT_markers$chr<-str_split_fixed(row.names(GAM_WT_vs_MUT_markers),pattern = "-",n = 2)[,1]
MUT_SIG<-GAM_WT_vs_MUT_markers[GAM_WT_vs_MUT_markers$padj<0.05 & GAM_WT_vs_MUT_markers$logFC < -0.2,]
WT_SIG<-GAM_WT_vs_MUT_markers[GAM_WT_vs_MUT_markers$padj<0.05 & GAM_WT_vs_MUT_markers$logFC>0.2,]

GAM_closest_features<-ClosestFeature(regions = StringToGRanges(row.names(GAM_WT_vs_MUT_markers)),object = ATAC_GAM)
GAM_WT_vs_MUT_markers$query_region<-row.names(GAM_WT_vs_MUT_markers)
GAM_WT_vs_MUT_markers<-merge(GAM_WT_vs_MUT_markers,GAM_closest_features,by="query_region",all.x = T)
GAM_WT_vs_MUT_markers<-na.omit(GAM_WT_vs_MUT_markers)
GAM_WT_vs_MUT_markers$gene_dist<-paste0(GAM_WT_vs_MUT_markers$gene_name,"_",GAM_WT_vs_MUT_markers$distance)

GAM_WT_vs_MUT_markers$delabel<-GAM_WT_vs_MUT_markers$padj<0.05 & GAM_WT_vs_MUT_markers$logFC < -0.1
GAM_WT_vs_MUT_markers$delabel[GAM_WT_vs_MUT_markers$delabel==TRUE]<-GAM_WT_vs_MUT_markers$gene_dist[GAM_WT_vs_MUT_markers$delabel==TRUE]
GAM_WT_vs_MUT_markers$delabel[GAM_WT_vs_MUT_markers$delabel==FALSE]<-NA
GAM_WT_vs_MUT_markers$delabel[GAM_WT_vs_MUT_markers$padj<0.05 & GAM_WT_vs_MUT_markers$logFC > 0.1] <- GAM_WT_vs_MUT_markers$gene_dist[GAM_WT_vs_MUT_markers$padj<0.05 & GAM_WT_vs_MUT_markers$logFC > 0.15]

GAM_WT_vs_MUT_markers$diffexpressed<-GAM_WT_vs_MUT_markers$padj<0.05 & GAM_WT_vs_MUT_markers$logFC < -0.1
GAM_WT_vs_MUT_markers$diffexpressed[GAM_WT_vs_MUT_markers$diffexpressed==TRUE]<-"WT"
GAM_WT_vs_MUT_markers$diffexpressed[GAM_WT_vs_MUT_markers$diffexpressed==FALSE]<-"not_DA"
GAM_WT_vs_MUT_markers$diffexpressed[GAM_WT_vs_MUT_markers$padj<0.05 & GAM_WT_vs_MUT_markers$logFC > 0.1] <-"MUT"

library(ggrepel)
# plot adding up all layers we have seen so far
p<-ggplot(data=GAM_WT_vs_MUT_markers, aes(x=logFC, y=-log10(padj), col=diffexpressed, label=delabel)) +
        geom_point() + 
        theme_minimal() +
        geom_text_repel() +
        scale_color_manual(values=c("blue", "black", "red")) +
        geom_vline(xintercept=c(-0.1, 0.1), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red")

#doGO
library(rGREAT)
library(BSgenome.Hsapiens.UCSC.hg38)
GO_DA <- function (input_peaks,background_peaks=NULL,genome="hg38")
{
if (is.null(background_peaks)){
job = submitGreatJob(input_peaks, species = genome) 
} else
{
job = submitGreatJob(input_peaks, species = genome,bg=background_peaks)
}
tb = getEnrichmentTables(job)
return(tb)
}


GAM_WT_vs_MUT_markers_DA<-GAM_WT_vs_MUT_markers[GAM_WT_vs_MUT_markers$diffexpressed!="not_DA",]
write.table(GAM_WT_vs_MUT_markers_DA,file = "./DA_GAMS_WT_vs_MUT.txt",quote = F,row.names = F,col.names = T,sep = "\t")


#remove x chrom diff
rowRanges = StringToGRanges(row.names(ATAC_GAM))
rowRanges<-rowRanges[as.character(seqnames(rowRanges)) %in% c(paste0("chr",1:22),"chrX","chrY")]
bg_peaks<-GRangesToString(rowRanges)
WT<-GAM_WT_vs_MUT_markers_DA[GAM_WT_vs_MUT_markers_DA$diffexpressed=="WT",]
MUT<-GAM_WT_vs_MUT_markers_DA[GAM_WT_vs_MUT_markers_DA$diffexpressed=="MUT",]


WT_GREAT<-GO_DA(input_peaks = StringToGRanges(WT$feature),background_peaks = rowRanges,genome = "hg38")

WT_plot<- ggplot(data = head(WT_GREAT$`GO Molecular Function`,n=20),aes(x=name,fill = Hyper_Fold_Enrichment,y=-log10(Hyper_Adjp_BH))) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))

MUT_GREAT<-GO_DA(input_peaks = StringToGRanges(MUT$feature),background_peaks = rowRanges,genome = "hg38")

MUT_plot<- ggplot(data = head(MUT_GREAT$`GO Molecular Function`,n=20),aes(x=name,fill = Hyper_Fold_Enrichment,y=-log10(Hyper_Adjp_BH))) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))


counts=GetAssayData(ATAC_GAM,slot = "counts",assay = "merge.peaks")[GRangesToString(rowRanges),]
ATAC_GAM[["peaks_corr"]]<-CreateChromatinAssay(
  counts = counts,
  sep = c("-", "-"),ranges = rowRanges,
  genome = 'hg38',
)

DefaultAssay(ATAC_GAM)<-"peaks_corr"

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
ATAC_GAM <- AddMotifs(
  object = ATAC_GAM,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

ATAC_GAM <- RunChromVAR(
  object = ATAC_GAM,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(ATAC_GAM) <- 'chromvar'

differential.activity <- FindMarkers(
  object = ATAC_GAM,group.by = "type",
  ident.1 = "WT",
  ident.2 = 'MUT',
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

MotifPlot(
  object = ATAC_GAM,
  motifs = head(rownames(differential.activity),n=10),
  assay = 'peaks_corr'
)

MotifPlot(
  object = ATAC_GAM,
  motifs = tail(rownames(differential.activity),n=10),
  assay = 'peaks_corr'
)


saveRDS(ATAC_GAM,file="./GAMS_all_ATAC_chromvar.rds")

