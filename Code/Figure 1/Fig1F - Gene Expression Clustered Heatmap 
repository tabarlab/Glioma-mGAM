library(Signac)
library(Seurat)
library(ComplexHeatmap)
library(patchwork)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(dplyr)

set.seed(1234)

################## Color Maps ##################
small_color_map <- c('#00B2CA','#FBDF88', '#DF5368')
################################################

RNA_GAMs <- readRDS("./RNA_GAMs.rds")

Idents(RNA_GAMs) <- RNA_GAMs$type

# Get GAM markers from each condition (WT, MUT, NB)
markers <- FindAllMarkers(RNA_GAMs, min.pct = 0.25, logfc.threshold = 0.25)

# Filter markers to top 50
markers %>%
    group_by(cluster) %>%
    slice_max(n = 50, order_by = avg_log2FC) -> markers.top.50

options(repr.plot.width=15, repr.plot.height=18)
DoHeatmap(subset(RNA_GAMs, downsample = 300), features = unique(markers.top.50$gene),
          size = 3, group.colors = small_color_map) -> p1
