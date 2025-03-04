library(Signac)
library(Seurat)
library(patchwork)
library(ggplot2)

set.seed(1234)

################## UMAPs ##################
################## Color Maps ##################
cluster_map <- c("#ff595e","#ffca3a","#ff924c","#c5ca30","#36949d","#8ac926","#1982c4","#565aa0","#4267ac","#6a4c93")
gbmpap_map <- c('#7FB2F0','#f78b88','#EFC164', '#F3835D', '#36175E','#4192D9','#F35955','#43A047','#00434C','#FD7400', '#00ABD8', '#F2B705')

################## Merging Datasets ##################
# Refer to Data processing folder
RNA_all <- readRDS('./RNA_all.rds')

################## RNA UMAPs ##################
p1 <- DimPlot(RNA_all, cols = c("#247ba0","#ffd97d","#60d394", "#ee6055","#ddb5e4"),pt.size = .6, shuffle = TRUE) + 
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
