library(data.table)
library(dplyr)
library(Signac)
library(Seurat)
library(viridis)
library(BuenColors)
library(EnsDb.Hsapiens.v86)
library(chromVAR)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg38)

counts <- Read10X_h5(filename = "../data/patient1/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "../data/patient1/singlecell.csv",
  header = TRUE,
  row.names = 1
)

# make mito mutation matrix
af1 <- readRDS("../output/dn_afs.rds"); colnames(af1) <- substr(start = 1, stop = 18, gsub("-1", "-1", colnames(af1)))
af2 <- readRDS("../output/dp_afs.rds"); colnames(af2) <- substr(start = 1, stop = 18, gsub("-1", "-2", colnames(af2)))
af3 <- readRDS("../output/mono_afs.rds"); colnames(af3) <- substr(start = 1, stop = 18, gsub("-1", "-3", colnames(af3)))
mito_mat <- data.matrix(cbind(cbind(af1, af2), af3))

cells <- intersect(colnames(counts), colnames(mito_mat))
metadata2 <- metadata[cells,]

chrom_assay <- CreateChromatinAssay(
  counts = counts[grepl("^chr", rownames(counts)),cells],
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = '../data/patient1/fragments.tsv.gz',
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

pA <- DimPlot(so_filt, group.by = c("channel"), label = FALSE) +
  scale_color_manual(values = c("#6FC276","#B19CD9","#A0D5F6")) +
  theme_void() + ggtitle("")
cowplot::ggsave2(pA, file = "../figures/3clusters.pdf", width = 2.2, height = 1.8)

if(FALSE){
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
}

# Look at mitochondrial mutations
so_filt[['mito']] <- CreateAssayObject(counts = mito_mat[,colnames(so_filt)])
DefaultAssay(so_filt) <- 'mito'

pX <- FeaturePlot(so_filt, c("9253G>A", "8780T>C"),
            sort.cell = TRUE, split.by = "channel", max.cutoff = 0.2) &
  scale_color_viridis() & theme_void() & theme(legend.position = "none")
cowplot::ggsave2(pX, file = "../figures/mito_vars2.png", width = 6, height = 4, dpi = 300)

# Do chromVAR scoring
gr_mtscatac <- so_filt@assays[["peaks"]]@ranges

fosl_diff_peak <- fread("../data/FOSL2_high_vs_low_da_peaks_pure_030624.csv")
fosl_diff_peak_up <- fosl_diff_peak %>% dplyr::filter(p_val_adj < 0.01 & avg_log2FC > 0)
fosl_diff_peak_down <- fosl_diff_peak %>% dplyr::filter(p_val_adj < 0.01 & avg_log2FC < 0)

ss_up <- stringr::str_split_fixed(fosl_diff_peak_up[["V1"]], "-", 3)
ss_down <- stringr::str_split_fixed(fosl_diff_peak_down[["V1"]], "-", 3)
data.frame(
  chr = ss_up[,1],
  start = as.numeric(ss_up[,2]),
  end = as.numeric(ss_up[,3])
) %>% makeGRangesFromDataFrame() -> gr_up

data.frame(
  chr = ss_down[,1],
  start = as.numeric(ss_down[,2]),
  end = as.numeric(ss_down[,3])
) %>% makeGRangesFromDataFrame() -> gr_down

ov_up <- findOverlaps(gr_up, gr_mtscatac)
ov_down <- findOverlaps(gr_down, gr_mtscatac)
matches <- data.frame(
  up = as.numeric(1:length(gr_mtscatac) %in% subjectHits(ov_up)),
  down = as.numeric(1:length(gr_mtscatac) %in% subjectHits(ov_down))
) %>% data.matrix()
SE <- SummarizedExperiment(
  assays = list(counts = so_filt@assays$peaks@counts),
  rowRanges = gr_mtscatac, colData = so_filt@meta.data 
)
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg38)
#motif_ix <- matchMotifs(human_pwms_v1, counts, genome = BSgenome.Hsapiens.UCSC.hg38)
dev <- computeDeviations(object = SE, annotations = matches)
so_filt$up_chromatin <- t(dev@assays@data$deviations)[,1]
so_filt$down_chromatin <- t(dev@assays@data$deviations)[,2]

FeaturePlot(so_filt, c("down_chromatin", "up_chromatin"), min.cutoff = "q10", max.cutoff = "q90") &
  scale_color_viridis()

wilcox.test(so_filt$up_chromatin[so_filt$channel == "DP"],
            so_filt$up_chromatin[so_filt$channel == "DN"]) %>%

pA <- ggplot(so_filt@meta.data, aes(x = channel, y = up_chromatin, color = channel)) +
  geom_violin(color = "lightgrey") + geom_boxplot(outlier.shape = NA, width = 0.3) +
  pretty_plot(fontsize = 7) + L_border() +
  labs(x = "sorted sample", y = "FOSL2 hi ATAC signature (deviation)") +
  scale_color_manual(values = c("#6FC276","#B19CD9","#A0D5F6")) +
  theme(legend.position = "none")
cowplot::ggsave2(pA, file = "../figures/chromatin_signature.pdf", width = 1.8, height = 1.8)
