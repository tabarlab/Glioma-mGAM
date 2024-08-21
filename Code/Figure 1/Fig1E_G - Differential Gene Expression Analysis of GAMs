library(Signac)
library(Seurat)
library(ComplexHeatmap)
library(patchwork)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(dplyr)

set.seed(1234)

################## Markers of GAMs from different tumor types (Related to Fig1E) ##################
RNA_GAMs <- readRDS("./RNA_GAMs.rds")

GAM_markers_WT <- FindMarkers(RNA_GAMs, `ident.1` = 'WT', min.pct = 0.25, logfc.threshold = 0.25)
GAM_markers_MUT <- FindMarkers(RNA_GAMs, `ident.1` = 'MUT', min.pct = 0.25, logfc.threshold = 0.25)
GAM_markers_NB <- FindMarkers(RNA_GAMs, `ident.1` = 'NB', min.pct = 0.25, logfc.threshold = 0.25)


################## Gene Expression Correlation for GAMs from different tumor types (Fig1G) ##################
GAM_markers_1 <- FindMarkers(RNA_GAMs, `ident.1` = 'WT', `ident.2` = 'NB', min.pct = 0.25, logfc.threshold = 0.25)
GAM_markers_2 <- FindMarkers(RNA_GAMs, `ident.1` = 'MUT', `ident.2` = 'NB', min.pct = 0.25, logfc.threshold = 0.25)
GAM_markers_3 <- FindMarkers(RNA_GAMs, `ident.1` = 'WT', `ident.2` = 'MUT', min.pct = 0.25, logfc.threshold = 0.25)

GAM_markers_1 <- GAM_markers_1[GAM_markers_1$p_val_adj <= 0.05,]
GAM_markers_2 <- GAM_markers_2[GAM_markers_2$p_val_adj <= 0.05,]
GAM_markers_3 <- GAM_markers_3[GAM_markers_3$p_val_adj <= 0.05,]

options(repr.plot.width=5, repr.plot.height=15)
ggplot2::theme_set(theme_cowplot())
Idents(RNA_GAMs) <- "type"
avg.RNA_GAMs <- as.data.frame(log10(AverageExpression(RNA_GAMs, verbose = FALSE)$RNA))
avg.RNA_GAMs$gene <- rownames(avg.RNA_GAMs)

# Filter MALAT1 and Mitocondrial for better visualization
avg.RNA_GAMs <- avg.RNA_GAMs[!grepl("MALAT1", avg.RNA_GAMs$gene),]
avg.RNA_GAMs <- avg.RNA_GAMs[!grepl("^MT-", avg.RNA_GAMs$gene),]

genes.to.label = c('NR4A3','TGFBI','ANXA1','CD163','FCGR1A','PLK2','C1QC','P2RY12','BHLHE41')
p1 <- ggplot2::ggplot(avg.RNA_GAMs[!is.infinite(rowSums(avg.RNA_GAMs[c('NB', 'WT')])),], ggplot2::aes(NB, WT)) + ggplot2::geom_point(col='grey', alpha = 0.5) + ggplot2::ggtitle("NB vs WT")
p1 <- p1 + ggplot2::geom_abline(slope = 1)
p1 <- p1 + ggpubr::stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)
p1 <- p1 + ggplot2::geom_point(data = avg.RNA_GAMs[avg.RNA_GAMs$gene %in% rownames(GAM_markers_1[GAM_markers_1$avg_log2FC <= -1,]),], ggplot2::aes(NB, WT), col='#FBDF88')
p1 <- p1 + ggplot2::geom_point(data = avg.RNA_GAMs[avg.RNA_GAMs$gene %in% rownames(GAM_markers_1[GAM_markers_1$avg_log2FC >= 1,]),], ggplot2::aes(NB, WT), col='#DF5368')
p1 <- LabelPoints(plot = p1, points = genes.to.label[1:4], repel = TRUE, color = 'black',  xnudge = -0.8, ynudge = 0.5)
p1 <- LabelPoints(plot = p1, points = genes.to.label[5:9], repel = TRUE, color = 'black',  xnudge = 0.5, ynudge = -0.5)
p1 <- p1 + ggplot2::xlab('Log10(Average Gene Expression NB)') + ggplot2::ylab('Log10(Average Gene Expression WT)')

genes.to.label = c('PDK4','CD83','NR4A3','CCL4','CCL3','MIF','AIF1','BHLHE41')
p2 <- ggplot2::ggplot(avg.RNA_GAMs[!is.infinite(rowSums(avg.RNA_GAMs[c('NB', 'MUT')])),], ggplot2::aes(NB, MUT)) + ggplot2::geom_point(col='grey', alpha = 0.5) + ggplot2::ggtitle("NB vs MUT")
p2 <- p2 + ggplot2::geom_abline(slope = 1)
p2 <- p2 + ggpubr::stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)
p2 <- p2 + ggplot2::geom_point(data = avg.RNA_GAMs[avg.RNA_GAMs$gene %in% rownames(GAM_markers_2[GAM_markers_2$avg_log2FC <= -1,]),], ggplot2::aes(NB, MUT), col='#FBDF88')
p2 <- p2 + ggplot2::geom_point(data = avg.RNA_GAMs[avg.RNA_GAMs$gene %in% rownames(GAM_markers_2[GAM_markers_2$avg_log2FC >= 1,]),], ggplot2::aes(NB, MUT), col='#00B2CA')
p2 <- LabelPoints(plot = p2, points = genes.to.label[1:5], repel = TRUE, color = 'black',  xnudge = -0.4, ynudge = 0.5)
p2 <- LabelPoints(plot = p2, points = genes.to.label[6:8], repel = TRUE, color = 'black',  xnudge = 0.6, ynudge = -0.5)
p2 <- p2 + ggplot2::xlab('Log10(Average Gene Expression NB)') + ggplot2::ylab('Log10(Average Gene Expression MUT)')

genes.to.label = c('PDK4','BHLHE41','VIM','CD44','VCAN','FCGR2B','LGALS1','SERPINE1','FN1','TGFBI')
p3 <- ggplot2::ggplot(avg.RNA_GAMs[!is.infinite(rowSums(avg.RNA_GAMs[c('WT', 'MUT')])),], ggplot2::aes(WT, MUT)) + ggplot2::geom_point(col='grey', alpha = 0.5) + ggplot2::ggtitle("WT vs MUT")
p3 <- p3 + ggplot2::geom_abline(slope = 1)
p3 <- p3 + ggpubr::stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)
p3 <- p3 + ggplot2::geom_point(data = avg.RNA_GAMs[avg.RNA_GAMs$gene %in% rownames(GAM_markers_3[GAM_markers_3$avg_log2FC <= -1,]),], ggplot2::aes(WT, MUT), col='#00B2CA')
p3 <- p3 + ggplot2::geom_point(data = avg.RNA_GAMs[avg.RNA_GAMs$gene %in% rownames(GAM_markers_3[GAM_markers_3$avg_log2FC >= 1,]),], ggplot2::aes(WT, MUT), col='#DF5368')
p3 <- LabelPoints(plot = p3, points = genes.to.label[1:2], repel = TRUE, color = 'black',  xnudge = -0.8, ynudge = 0.5)
p3 <- LabelPoints(plot = p3, points = genes.to.label[3:10], repel = TRUE, color = 'black',  xnudge = 0, ynudge = -1)
p3 <- p3 + ggplot2::xlab('Log10(Average Gene Expression WT)') + ggplot2::ylab('Log10(Average Gene Expression MUT)')

p_all <- p1 / p2 / p3
p_all