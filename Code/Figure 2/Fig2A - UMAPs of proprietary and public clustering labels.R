library(Signac)
library(Seurat)
library(ComplexHeatmap)
library(patchwork)
library(ggplot2)

set.seed(1234)

################## UMAPs ##################
################## Color Maps ##################
cluster_map <- c("#ff595e","#ffca3a","#ff924c","#c5ca30","#36949d","#8ac926","#1982c4","#565aa0","#4267ac","#6a4c93")
gbmpap_map <- c('#7FB2F0','#f78b88','#EFC164', '#F3835D', '#36175E','#4192D9','#F35955','#43A047','#00434C','#FD7400', '#00ABD8', '#F2B705')

################## Metadata labels ##################
# Done through Azimuth, refer to Data processing folder
RNA_GAMs <- readRDS('./Azimuth/RNA_GAMs_Azimuth_labels.rds')

################## RNA UMAPs ##################
p1 <- DimPlot(RNA_GAMs, pt.size = 1, group.by = 'RNA_snn_res.0.3', cols = cluster_map) + 
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

p2 <- DimPlot(RNA_GAMs, pt.size = 1, group.by ='ABD_cluster') + 
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

p3 <- DimPlot(RNA_GAMs, pt.size = 1, group.by ='predicted.annotation_level_4', cols = rev(gbmpap_map)) + 
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