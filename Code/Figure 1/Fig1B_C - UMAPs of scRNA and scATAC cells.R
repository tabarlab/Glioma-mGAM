library(Signac)
library(Seurat)
library(ComplexHeatmap)
library(patchwork)
library(ggplot2)

set.seed(1234)

################## UMAPs ##################
################## Color Maps ##################
small_color_map <- c('#00B2CA','#FBDF88', '#DF5368')
large_color_map <- c("#C91D42", "#2E45B8", "#1DC9A4", "#F97A1F", "#F9C31F", 
                     "#e67889", "#9CA9E7", "#93F0DB", "#FCB583", "#fcde88")

################## RNA UMAPs ##################
RNA_all <- readRDS("./RNA_all.rds")

p1 <- DimPlot(RNA_all, cols = small_color_map, group.by = "type") + 
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(), 
              axis.text.y=element_blank(),
              axis.ticks=element_blank(), 
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.key.size = unit(10, 'cm'), #change legend key size
              legend.key.height = unit(1.2, 'cm'), #change legend key height
              legend.key.width = unit(0.5, 'cm'), #change legend key width
              legend.title = element_text(size=14), #change legend title font size
              legend.text = element_text(size=20, face = 'bold'), #change legend text font size 
              plot.title = element_blank())

p2 <- DimPlot(RNA_all, cols = large_color_map, group.by = "celltype", pt.size = 1) + 
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(), 
              axis.text.y=element_blank(),
              axis.ticks=element_blank(), 
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.key.size = unit(10, 'cm'), #change legend key size
              legend.key.height = unit(1.2, 'cm'), #change legend key height
              legend.key.width = unit(0.5, 'cm'), #change legend key width
              legend.title = element_text(size=14), #change legend title font size
              legend.text = element_text(size=20, face = 'bold'), #change legend text font size 
              plot.title = element_blank())

################## ATAC UMAPs ##################
ATAC_all <- readRDS("./ATAC_all.rds")

p3 <- DimPlot(ATAC_all, cols = small_color_map, group.by = "type") + 
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(), 
              axis.text.y=element_blank(),
              axis.ticks=element_blank(), 
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.key.size = unit(10, 'cm'), #change legend key size
              legend.key.height = unit(1.2, 'cm'), #change legend key height
              legend.key.width = unit(0.5, 'cm'), #change legend key width
              legend.title = element_text(size=14), #change legend title font size
              legend.text = element_text(size=20, face = 'bold'), #change legend text font size 
              plot.title = element_blank())

p4 <- DimPlot(ATAC_all, cols = large_color_map, group.by = "celltype", pt.size = 1) + 
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(), 
              axis.text.y=element_blank(),
              axis.ticks=element_blank(), 
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.key.size = unit(10, 'cm'), #change legend key size
              legend.key.height = unit(1.2, 'cm'), #change legend key height
              legend.key.width = unit(0.5, 'cm'), #change legend key width
              legend.title = element_text(size=14), #change legend title font size
              legend.text = element_text(size=20, face = 'bold'), #change legend text font size 
              plot.title = element_blank())

################## Dot Plot ##################
################## Gene List ##################
markers <- c("PDGFRA","SOX2","NES","OLIG1","GFAP","S100B","EGFR","CD4","GZMB",
             "CD8A","CD3D","CD3E","PECAM1","VWF","CD74","AIF1","CD14","CSF1R",
             "TREM2","CD163","PTPRC","CD68", "CX3CR1","LYZ","ITGAX","SEPP1",
             "MERTK","SPI1","TGFBI","TMEM119","P2RY12","ITGAM","CD33")

################## RNA Dotplot ##################
p5 <- DotPlot(RNA_all, features = markers, group.by="celltype") + RotatedAxis()

