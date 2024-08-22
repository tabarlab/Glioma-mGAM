library(harmony)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(Matrix)
library(patchwork)
library(dplyr)
library(ComplexHeatmap)
library(reshape2)

set.seed(1234)

################## Import and Processing of Seurat Objects ##################

################## MSK 101 ##################
genes<-read.csv("./MSK101_RNA/2502_MSK101_IGO_11874_20_sparse_counts_genes.csv",header = F)
barc<-read.csv("./MSK101_RNA/2502_MSK101_IGO_11874_20_sparse_counts_barcodes.csv",header = F)
matrix<-read.table(file = "./MSK101_RNA/2502_MSK101_IGO_11874_20_sparse_read_counts.mtx",skip = 3,header = F)

barc$cell<-paste("Cell",barc$V1,sep=".")

IN<-sparseMatrix(i=matrix[,2],j=matrix[,1],x=matrix[,3])
colnames(IN)<-barc$cell
row.names(IN)<-genes$V2


RNA_MSK101<-CreateSeuratObject(counts = IN,assay = "RNA",min.cells = 3, min.features = 200)
RNA_MSK101[["percent.mt"]] <- PercentageFeatureSet(RNA_MSK101, pattern = "^MT-")
RNA_MSK101$percent.rp <- PercentageFeatureSet(RNA_MSK101,pattern="^(RPS|RPL)")

RNA_MSK101$batch = "1"
RNA_MSK101$type = "WT"
RNA_MSK101$patient = "MSK101"
RNA_MSK101$histo = "gbm"
RNA_MSK101$grade = "4"

VlnPlot(RNA_MSK101, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
FeatureScatter(RNA_MSK101, feature1 = "nCount_RNA", feature2 = "percent.mt")
RNA_MSK101 <- subset(RNA_MSK101, subset = nFeature_RNA > 200 & nFeature_RNA < 7500  & percent.mt < 25)

saveRDS(RNA_MSK101,file="RNA_MSK101.rds")

################## MSK 102 ##################
genes<-read.csv("./MSK102_RNA/2503_MSK102_IGO_11874_21_sparse_counts_genes.csv",header = F)
barc<-read.csv("./MSK102_RNA/2503_MSK102_IGO_11874_21_sparse_counts_barcodes.csv",header = F)
matrix<-read.table(file = "./MSK102_RNA/2503_MSK102_IGO_11874_21_sparse_read_counts.mtx",skip = 3,header = F)

barc$cell<-paste("Cell",barc$V1,sep=".")

IN<-sparseMatrix(i=matrix[,2],j=matrix[,1],x=matrix[,3])
colnames(IN)<-barc$cell
row.names(IN)<-genes$V2

RNA_MSK102<-CreateSeuratObject(counts = IN,assay = "RNA",min.cells = 3, min.features = 200)
RNA_MSK102[["percent.mt"]] <- PercentageFeatureSet(RNA_MSK102, pattern = "^MT-")
RNA_MSK102$percent.rp <- PercentageFeatureSet(RNA_MSK102,pattern="^(RPS|RPL)")

RNA_MSK102$batch = "1"
RNA_MSK102$type = "WT"
RNA_MSK102$patient = "MSK102"
RNA_MSK102$histo = "gbm"
RNA_MSK102$grade = "4"


VlnPlot(RNA_MSK102, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)

VlnPlot(RNA_MSK102, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(RNA_MSK102, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
RNA_MSK102 <- subset(RNA_MSK102, subset = nFeature_RNA > 200 & nFeature_RNA < 7500  & percent.mt <25)


saveRDS(RNA_MSK102,file="RNA_MSK102.rds")

################## MSK 108 ##################
genes<-read.csv("./MSK108_RNA/2504_MSK108_IGO_11874_22_sparse_counts_genes.csv",header = F)
barc<-read.csv("./MSK108_RNA/2504_MSK108_IGO_11874_22_sparse_counts_barcodes.csv",header = F)
matrix<-read.table(file = "./MSK108_RNA/2504_MSK108_IGO_11874_22_sparse_read_counts.mtx",skip = 3,header = F)

barc$cell<-paste("Cell",barc$V1,sep=".")

IN<-sparseMatrix(i=matrix[,2],j=matrix[,1],x=matrix[,3])
colnames(IN)<-barc$cell
row.names(IN)<-genes$V2


RNA_MSK108<-CreateSeuratObject(counts = IN,assay = "RNA",min.cells = 3, min.features = 200)
RNA_MSK108[["percent.mt"]] <- PercentageFeatureSet(RNA_MSK108, pattern = "^MT-")
RNA_MSK108$percent.rp <- PercentageFeatureSet(RNA_MSK108,pattern="^(RPS|RPL)")

RNA_MSK108$batch = "1"
RNA_MSK108$type = "MUT"
RNA_MSK108$patient = "MSK108"
RNA_MSK108$histo = "astro"
RNA_MSK108$grade = "3"

VlnPlot(RNA_MSK108, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
FeatureScatter(RNA_MSK108, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
RNA_MSK108 <- subset(RNA_MSK108, subset = nFeature_RNA > 200 & nFeature_RNA < 7500  & percent.mt < 25)

saveRDS(RNA_MSK108,file="RNA_MSK108.rds")

################## MSK 109 ##################
genes<-read.csv("./MSK109_RNA/2505_MSK109-PM_IGO_11874_23_sparse_counts_genes.csv",header = F)
barc<-read.csv("./MSK109_RNA/2505_MSK109-PM_IGO_11874_23_sparse_counts_barcodes.csv",header = F)
matrix<-read.table(file = "./MSK109_RNA/2505_MSK109-PM_IGO_11874_23_sparse_read_counts.mtx",skip = 3,header = F)

barc$cell<-paste("Cell",barc$V1,sep=".")

IN<-sparseMatrix(i=matrix[,2],j=matrix[,1],x=matrix[,3])
colnames(IN)<-barc$cell
row.names(IN)<-genes$V2


RNA_MSK109<-CreateSeuratObject(counts = IN,assay = "RNA",min.cells = 3, min.features = 200)
RNA_MSK109[["percent.mt"]] <- PercentageFeatureSet(RNA_MSK109, pattern = "^MT-")
RNA_MSK109$percent.rp <- PercentageFeatureSet(RNA_MSK109,pattern="^(RPS|RPL)")


RNA_MSK109$batch = "1"
RNA_MSK109$type = "MUT"
RNA_MSK109$patient = "MSK109"
RNA_MSK109$histo = "astro"
RNA_MSK109$grade = "2"


VlnPlot(RNA_MSK109, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
FeatureScatter(RNA_MSK109, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
RNA_MSK109 <- subset(RNA_MSK109, subset = nFeature_RNA > 200 & nFeature_RNA < 7500  & percent.mt < 25)

saveRDS(RNA_MSK109,file="RNA_MSK109.rds")

################## MSK 110 ##################
genes<-read.csv("./RNA/2984_MSK110_IGO_12381_28_sparse_counts_genes.csv",header = F)
barc<-read.csv("./RNA/2984_MSK110_IGO_12381_28_sparse_counts_barcodes.csv",header = F)
matrix<-read.table(file = "./RNA/2984_MSK110_IGO_12381_28_sparse_read_counts.mtx",skip = 3,header = F)

barc$cell<-paste("Cell",barc$V1,sep=".")

IN<-sparseMatrix(i=matrix[,2],j=matrix[,1],x=matrix[,3])
colnames(IN)<-barc$cell
row.names(IN)<-genes$V2


RNA_MSK110<-CreateSeuratObject(counts = IN,assay = "RNA",min.cells = 3, min.features = 200)
RNA_MSK110[["percent.mt"]] <- PercentageFeatureSet(RNA_MSK110, pattern = "^MT-")
RNA_MSK110$percent.rp <- PercentageFeatureSet(RNA_MSK110,pattern="^(RPS|RPL)")

RNA_MSK110$patient<-"MSK110"
RNA_MSK110$type<- "MUT"
RNA_MSK110$batch<-"2"
RNA_MSK110$histo <- "oligo"
RNA_MSK110$grade <- "3"


VlnPlot(RNA_MSK110, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
FeatureScatter(RNA_MSK110, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
RNA_MSK110 <- subset(RNA_MSK110, subset = nFeature_RNA > 200 & nFeature_RNA < 7500  & percent.mt < 25)

saveRDS(RNA_MSK110,file="RNA_MSK110.rds")

################## MSK 114 ##################
genes<-read.csv("./RNA/2989_MSK114_IGO_12381_30_sparse_counts_genes.csv",header = F)
barc<-read.csv("./RNA/2989_MSK114_IGO_12381_30_sparse_counts_barcodes.csv",header = F)
matrix<-read.table(file = "./RNA/2989_MSK114_IGO_12381_30_sparse_read_counts.mtx",skip = 3,header = F)

barc$cell<-paste("Cell",barc$V1,sep=".")

IN<-sparseMatrix(i=matrix[,2],j=matrix[,1],x=matrix[,3])
colnames(IN)<-barc$cell
row.names(IN)<-genes$V2


RNA_MSK114<-CreateSeuratObject(counts = IN,assay = "RNA",min.cells = 3, min.features = 200)
RNA_MSK114[["percent.mt"]] <- PercentageFeatureSet(RNA_MSK114, pattern = "^MT-")
RNA_MSK114$percent.rp <- PercentageFeatureSet(RNA_MSK114,pattern="^(RPS|RPL)")

RNA_MSK114$patient<-"MSK114"
RNA_MSK114$type<- "MUT"
RNA_MSK114$batch<- "2"
RNA_MSK114$histo<- "astro"
RNA_MSK114$grade<- "4"

VlnPlot(RNA_MSK114, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
FeatureScatter(RNA_MSK114, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
RNA_MSK114 <- subset(RNA_MSK114, subset = nFeature_RNA > 200 & nFeature_RNA < 7500  & percent.mt < 25)

saveRDS(RNA_MSK114,file="RNA_MSK114.rds")

################## MSK 111 ##################
genes<-read.csv("./RNA/2986_MSK111_IGO_12381_27_sparse_counts_genes.csv",header = F)
barc<-read.csv("./RNA/2986_MSK111_IGO_12381_27_sparse_counts_barcodes.csv",header = F)
matrix<-read.table(file = "./RNA/2986_MSK111_IGO_12381_27_sparse_read_counts.mtx",skip = 3,header = F)

barc$cell<-paste("Cell",barc$V1,sep=".")

IN<-sparseMatrix(i=matrix[,2],j=matrix[,1],x=matrix[,3])
colnames(IN)<-barc$cell
row.names(IN)<-genes$V2


RNA_MSK111<-CreateSeuratObject(counts = IN,assay = "RNA",min.cells = 3, min.features = 200)
RNA_MSK111[["percent.mt"]] <- PercentageFeatureSet(RNA_MSK111, pattern = "^MT-")
RNA_MSK111$percent.rp <- PercentageFeatureSet(RNA_MSK111,pattern="^(RPS|RPL)")

RNA_MSK111$patient<-"MSK111"
RNA_MSK111$type<- "MUT"
RNA_MSK111$batch<-"2"
RNA_MSK111$histo<-"astro"
RNA_MSK111$grade<-"3"


VlnPlot(RNA_MSK111, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
FeatureScatter(RNA_MSK111, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
RNA_MSK111 <- subset(RNA_MSK111, subset = nFeature_RNA > 200 & nFeature_RNA < 7500  & percent.mt < 25)

saveRDS(RNA_MSK111,file="RNA_MSK111.rds")

################## MSK 112 ##################
genes<-read.csv("./RNA/2987_MSK112_IGO_12381_31_sparse_counts_genes.csv",header = F)
barc<-read.csv("./RNA/2987_MSK112_IGO_12381_31_sparse_counts_barcodes.csv",header = F)
matrix<-read.table(file = "./RNA/2987_MSK112_IGO_12381_31_sparse_read_counts.mtx",skip = 3,header = F)

barc$cell<-paste("Cell",barc$V1,sep=".")

IN<-sparseMatrix(i=matrix[,2],j=matrix[,1],x=matrix[,3])
colnames(IN)<-barc$cell
row.names(IN)<-genes$V2


RNA_MSK112<-CreateSeuratObject(counts = IN,assay = "RNA",min.cells = 3, min.features = 200)
RNA_MSK112[["percent.mt"]] <- PercentageFeatureSet(RNA_MSK112, pattern = "^MT-")
RNA_MSK112$percent.rp <- PercentageFeatureSet(RNA_MSK112,pattern="^(RPS|RPL)")

RNA_MSK112$patient<-"MSK112"
RNA_MSK112$type<- "MUT"
RNA_MSK112$batch<-"2"
RNA_MSK112$histo<-"astro"
RNA_MSK112$grade<-"3"

VlnPlot(RNA_MSK112, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
FeatureScatter(RNA_MSK112, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
RNA_MSK112 <- subset(RNA_MSK112, subset = nFeature_RNA > 200 & nFeature_RNA < 7500  & percent.mt < 25)

saveRDS(RNA_MSK112,file="RNA_MSK112.rds")

################## MSK 113 ##################
genes<-read.csv("./RNA/2988_MSK113_IGO_12381_32_sparse_counts_genes.csv",header = F)
barc<-read.csv("./RNA/2988_MSK113_IGO_12381_32_sparse_counts_barcodes.csv",header = F)
matrix<-read.table(file = "./RNA/2988_MSK113_IGO_12381_32_sparse_read_counts.mtx",skip = 3,header = F)

barc$cell<-paste("Cell",barc$V1,sep=".")

IN<-sparseMatrix(i=matrix[,2],j=matrix[,1],x=matrix[,3])
colnames(IN)<-barc$cell
row.names(IN)<-genes$V2


RNA_MSK113<-CreateSeuratObject(counts = IN,assay = "RNA",min.cells = 3, min.features = 200)
RNA_MSK113[["percent.mt"]] <- PercentageFeatureSet(RNA_MSK113, pattern = "^MT-")
RNA_MSK113$percent.rp <- PercentageFeatureSet(RNA_MSK113,pattern="^(RPS|RPL)")

RNA_MSK113$patient<-"MSK113"
RNA_MSK113$type<- "MUT"
RNA_MSK113$batch<- "2"
RNA_MSK113$histo<- "astro"
RNA_MSK113$grade<- "4"

VlnPlot(RNA_MSK113, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
FeatureScatter(RNA_MSK113, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
RNA_MSK113 <- subset(RNA_MSK113, subset = nFeature_RNA > 200 & nFeature_RNA < 7500  & percent.mt < 25)

saveRDS(RNA_MSK113,file="RNA_MSK113.rds")

################## MSK 103 ##################
genes<-read.csv("./RNA/2990_MSK103_IGO_12381_33_sparse_counts_genes.csv",header = F)
barc<-read.csv("./RNA/2990_MSK103_IGO_12381_33_sparse_counts_barcodes.csv",header = F)
matrix<-read.table(file = "./RNA/2990_MSK103_IGO_12381_33_sparse_read_counts.mtx",skip = 3,header = F)

barc$cell<-paste("Cell",barc$V1,sep=".")

IN<-sparseMatrix(i=matrix[,2],j=matrix[,1],x=matrix[,3])
colnames(IN)<-barc$cell
row.names(IN)<-genes$V2


RNA_MSK103<-CreateSeuratObject(counts = IN,assay = "RNA",min.cells = 3, min.features = 200)
RNA_MSK103[["percent.mt"]] <- PercentageFeatureSet(RNA_MSK103, pattern = "^MT-")
RNA_MSK103$percent.rp <- PercentageFeatureSet(RNA_MSK103,pattern="^(RPS|RPL)")

RNA_MSK103$patient<-"MSK103"
RNA_MSK103$type<- "WT"
RNA_MSK103$batch<- "2"
RNA_MSK103$histo<- "gbm"
RNA_MSK103$grade<- "4"


VlnPlot(RNA_MSK103, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
FeatureScatter(RNA_MSK103, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
RNA_MSK103 <- subset(RNA_MSK103, subset = nFeature_RNA > 200 & nFeature_RNA < 7500  & percent.mt < 25)

saveRDS(RNA_MSK103,file="RNA_MSK103.rds")

################## MSK 115 ##################
genes<-read.csv("./RNA/2985_MSK115_IGO_12381_29_sparse_counts_genes.csv",header = F)
barc<-read.csv("./RNA/2985_MSK115_IGO_12381_29_sparse_counts_barcodes.csv",header = F)
matrix<-read.table(file = "./RNA/2985_MSK115_IGO_12381_29_sparse_read_counts.mtx",skip = 3,header = F)

barc$cell<-paste("Cell",barc$V1,sep=".")

IN<-sparseMatrix(i=matrix[,2],j=matrix[,1],x=matrix[,3])
colnames(IN)<-barc$cell
row.names(IN)<-genes$V2


RNA_MSK115<-CreateSeuratObject(counts = IN,assay = "RNA",min.cells = 3, min.features = 200)
RNA_MSK115[["percent.mt"]] <- PercentageFeatureSet(RNA_MSK115, pattern = "^MT-")
RNA_MSK115$percent.rp <- PercentageFeatureSet(RNA_MSK115,pattern="^(RPS|RPL)")

RNA_MSK115$patient<-"MSK115"
RNA_MSK115$type<- "MUT"
RNA_MSK115$batch<- "2"
RNA_MSK115$histo<- "oligo"
RNA_MSK115$grade<- "3"

VlnPlot(RNA_MSK115, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
FeatureScatter(RNA_MSK115, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
RNA_MSK115 <- subset(RNA_MSK115, subset = nFeature_RNA > 200 & nFeature_RNA < 7500  & percent.mt < 25)

saveRDS(RNA_MSK115,file="RNA_MSK115.rds")

################## MSK 107 ##################
genes<-read.csv("./RNA/3324_KY-1222_MSK107_IGO_12437_J_7_sparse_counts_genes.csv",header = F)
barc<-read.csv("./RNA/3324_KY-1222_MSK107_IGO_12437_J_7_sparse_counts_barcodes.csv",header = F)
matrix<-read.table(file = "./RNA/3324_KY-1222_MSK107_IGO_12437_J_7_sparse_read_counts.mtx",skip = 3,header = F)

barc$cell<-paste("Cell",barc$V1,sep=".")

IN<-sparseMatrix(i=matrix[,2],j=matrix[,1],x=matrix[,3])
colnames(IN)<-barc$cell
row.names(IN)<-genes$V2


RNA_MSK107<-CreateSeuratObject(counts = IN,assay = "RNA",min.cells = 3, min.features = 200)
RNA_MSK107[["percent.mt"]] <- PercentageFeatureSet(RNA_MSK107, pattern = "^MT-")
RNA_MSK107$percent.rp <- PercentageFeatureSet(RNA_MSK107,pattern="^(RPS|RPL)")

RNA_MSK107$patient<-"MSK107"
RNA_MSK107$type<- "WT"
RNA_MSK107$batch<- "3"
RNA_MSK107$histo<- "gbm"
RNA_MSK107$grade<- "4"

VlnPlot(RNA_MSK107, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
FeatureScatter(RNA_MSK107, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
RNA_MSK107 <- subset(RNA_MSK107, subset = nFeature_RNA > 200 & nFeature_RNA < 7500  & percent.mt < 25)

saveRDS(RNA_MSK107,file="RNA_MSK107.rds")

################## MSK 106 ##################
genes<-read.csv("./RNA/3325_KY-1222_MSK106_IGO_12437_J_8_sparse_counts_genes.csv",header = F)
barc<-read.csv("./RNA/3325_KY-1222_MSK106_IGO_12437_J_8_sparse_counts_barcodes.csv",header = F)
matrix<-read.table(file = "./RNA/3325_KY-1222_MSK106_IGO_12437_J_8_sparse_read_counts.mtx",skip = 3,header = F)

barc$cell<-paste("Cell",barc$V1,sep=".")

IN<-sparseMatrix(i=matrix[,2],j=matrix[,1],x=matrix[,3])
colnames(IN)<-barc$cell
row.names(IN)<-genes$V2


RNA_MSK106<-CreateSeuratObject(counts = IN,assay = "RNA",min.cells = 3, min.features = 200)
RNA_MSK106[["percent.mt"]] <- PercentageFeatureSet(RNA_MSK106, pattern = "^MT-")
RNA_MSK106$percent.rp <- PercentageFeatureSet(RNA_MSK106,pattern="^(RPS|RPL)")

RNA_MSK106$patient<-"MSK106"
RNA_MSK106$type<- "WT"
RNA_MSK106$batch<- "3"
RNA_MSK106$histo<- "gbm"
RNA_MSK106$grade<- "4"

VlnPlot(RNA_MSK106, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
FeatureScatter(RNA_MSK106, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
RNA_MSK106 <- subset(RNA_MSK106, subset = nFeature_RNA > 200 & nFeature_RNA < 7500  & percent.mt < 25)

saveRDS(RNA_MSK106,file="RNA_MSK106.rds")

################## MSK 105 ##################
genes<-read.csv("./RNA/3326_KY-1222_MSK105_IGO_12437_J_9_sparse_counts_genes.csv",header = F)
barc<-read.csv("./RNA/3326_KY-1222_MSK105_IGO_12437_J_9_sparse_counts_barcodes.csv",header = F)
matrix<-read.table(file = "./RNA/3326_KY-1222_MSK105_IGO_12437_J_9_sparse_read_counts.mtx",skip = 3,header = F)

barc$cell<-paste("Cell",barc$V1,sep=".")

IN<-sparseMatrix(i=matrix[,2],j=matrix[,1],x=matrix[,3])
colnames(IN)<-barc$cell
row.names(IN)<-genes$V2


RNA_MSK105<-CreateSeuratObject(counts = IN,assay = "RNA",min.cells = 3, min.features = 200)
RNA_MSK105[["percent.mt"]] <- PercentageFeatureSet(RNA_MSK105, pattern = "^MT-")
RNA_MSK105$percent.rp <- PercentageFeatureSet(RNA_MSK105,pattern="^(RPS|RPL)")

RNA_MSK105$patient<-"MSK105"
RNA_MSK105$type<- "WT"
RNA_MSK105$batch<- "3"
RNA_MSK105$histo<- "gbm"
RNA_MSK105$grade<- "4"

VlnPlot(RNA_MSK105, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
FeatureScatter(RNA_MSK105, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
RNA_MSK105 <- subset(RNA_MSK105, subset = nFeature_RNA > 200 & nFeature_RNA < 7500  & percent.mt < 25)

saveRDS(RNA_MSK105,file="RNA_MSK105.rds")

################## MSK 104 ##################
genes<-read.csv("./RNA/3327_KY-1222_MSK104_IGO_12437_L_16_sparse_counts_genes.csv",header = F)
barc<-read.csv("./RNA/3327_KY-1222_MSK104_IGO_12437_L_16_sparse_counts_barcodes.csv",header = F)
matrix<-read.table(file = "./RNA/3327_KY-1222_MSK104_IGO_12437_L_16_sparse_read_counts.mtx",skip = 3,header = F)

barc$cell<-paste("Cell",barc$V1,sep=".")

IN<-sparseMatrix(i=matrix[,2],j=matrix[,1],x=matrix[,3])
colnames(IN)<-barc$cell
row.names(IN)<-genes$V2


RNA_MSK104<-CreateSeuratObject(counts = IN,assay = "RNA",min.cells = 3, min.features = 200)
RNA_MSK104[["percent.mt"]] <- PercentageFeatureSet(RNA_MSK104, pattern = "^MT-")
RNA_MSK104$percent.rp <- PercentageFeatureSet(RNA_MSK104,pattern="^(RPS|RPL)")

RNA_MSK104$patient<-"MSK104"
RNA_MSK104$type<- "WT"
RNA_MSK104$batch<- "3"
RNA_MSK104$histo<- "gbm"
RNA_MSK104$grade<- "4"

VlnPlot(RNA_MSK104, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
FeatureScatter(RNA_MSK104, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
RNA_MSK104 <- subset(RNA_MSK104, subset = nFeature_RNA > 200 & nFeature_RNA < 7500  & percent.mt < 25)

saveRDS(RNA_MSK104,file="RNA_MSK104.rds")

################## MSK NB 1 ##################
genes<-read.csv("./RNA/NB1/3815_KY-1666_NB1_IGO_12437_AO_1_sparse_counts_genes.csv",header = F)
barc<-read.csv("./RNA/NB1/3815_KY-1666_NB1_IGO_12437_AO_1_sparse_counts_barcodes.csv",header = F)
matrix<-read.table(file = "./RNA/NB1/3815_KY-1666_NB1_IGO_12437_AO_1_sparse_molecule_counts.mtx",skip = 3,header = F)

barc$cell<-paste("Cell",barc$V1,sep=".")

IN<-sparseMatrix(i=matrix[,2],j=matrix[,1],x=matrix[,3])
colnames(IN)<-barc$cell
row.names(IN)<-genes$V2

RNA_NB1<-CreateSeuratObject(counts = IN,assay = "RNA",min.cells = 3, min.features = 200)
RNA_NB1[["percent.mt"]] <- PercentageFeatureSet(RNA_NB1, pattern = "^MT-")
RNA_NB1$percent.rp <- PercentageFeatureSet(RNA_NB1,pattern="^(RPS|RPL)")

RNA_NB1$patient<-"NB1"
RNA_NB1$type<- "NB"
RNA_NB1$batch<-"4"
RNA_NB1$histo<-"NB"
RNA_NB1$grade<-"NB"

VlnPlot(RNA_NB1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
FeatureScatter(RNA_NB1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
RNA_NB1 <- subset(RNA_NB1, subset = nFeature_RNA > 200 & nFeature_RNA < 7500  & percent.mt < 25)

saveRDS(RNA_NB1,file="RNA_NB1.rds")

################## MSK NB 2 ##################
genes<-read.csv("./RNA/NB2/3816_KY-1666_NB2_IGO_12437_AO_2_sparse_counts_genes.csv",header = F)
barc<-read.csv("./RNA/NB2/3816_KY-1666_NB2_IGO_12437_AO_2_sparse_counts_barcodes.csv",header = F)
matrix<-read.table(file = "./RNA/NB2/3816_KY-1666_NB2_IGO_12437_AO_2_sparse_molecule_counts.mtx",skip = 3,header = F)

barc$cell<-paste("Cell",barc$V1,sep=".")

IN<-sparseMatrix(i=matrix[,2],j=matrix[,1],x=matrix[,3])
colnames(IN)<-barc$cell
row.names(IN)<-genes$V2

RNA_NB2<-CreateSeuratObject(counts = IN,assay = "RNA",min.cells = 3, min.features = 200)
RNA_NB2[["percent.mt"]] <- PercentageFeatureSet(RNA_NB2, pattern = "^MT-")
RNA_NB2$percent.rp <- PercentageFeatureSet(RNA_NB2,pattern="^(RPS|RPL)")

RNA_NB2$patient<-"NB2"
RNA_NB2$type<- "NB"
RNA_NB2$batch<-"4"
RNA_NB2$histo<-"NB"
RNA_NB2$grade<-"NB"


VlnPlot(RNA_NB2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
FeatureScatter(RNA_NB2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
RNA_NB2 <- subset(RNA_NB2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500  & percent.mt < 25)

saveRDS(RNA_NB2,file="RNA_NB2.rds")

#######################################
RNA_all<-merge(x=RNA_20201012,y=(list(RNA_20201013,RNA_20201116,RNA_20201117,RNA_0218,RNA_0308,RNA_0330,RNA_0405,RNA_0603,RNA_1119,RNA_117am,RNA_0413,RNA_0628,RNA_0416,RNA_0428)),add.cell.ids = c("RNA_12","RNA_13","RNA_16","RNA_17PM","RNA_18","RNA_08","RNA_30","RNA_05","RNA_03","RNA_19","RNA_17am","RNA_0413","RNA_0628","RNA_0416","RNA_0428"))

# Run the standard workflow for visualization and clustering
RNA_all <- NormalizeData(RNA_all)
RNA_all <- FindVariableFeatures(RNA_all, selection.method = "vst", nfeatures = 2000)
RNA_all <- ScaleData(RNA_all, verbose = FALSE)
RNA_all <- RunPCA(RNA_all, npcs = 30, verbose = FALSE)
RNA_all <- RunUMAP(RNA_all, reduction = "pca", dims = 1:30)
RNA_all <- FindNeighbors(RNA_all, reduction = "pca", dims = 1:30)
RNA_all <- FindClusters(RNA_all, resolution = 0.5)


saveRDS(RNA_all,file="./RNA_all_unint.rds")


p1<-DimPlot(object = RNA_all, label = FALSE) + theme(axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), axis.title.x=element_blank(),
          axis.title.y=element_blank(), legend.position = 'none', plot.title = element_blank()) #, raster = FALSE, group.by = 'celltype.stim') 
p2<-DimPlot(object = RNA_all, label = FALSE, group.by = "type") + theme(axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), axis.title.x=element_blank(),
          axis.title.y=element_blank(), legend.position = 'none', plot.title = element_blank())
p3<-DimPlot(object = RNA_all, label = TRUE, group.by = "batch")
p4<-DimPlot(object = RNA_all, label = TRUE, group.by = "histo")
p5<-DimPlot(object = RNA_all, label = TRUE, group.by = "grade")
p6<-DimPlot(object = RNA_all, label = TRUE, group.by = "patient")


# Subset and subclustering
RNA_GAMs <- subset(RNA_all, idents = c('combined_GAMs_WT', 'combined_GAMs_MUT'))
RNA_GAMs <- readRDS("/data/tabar/SingleCell/RNA_GAMs_2.rds")

# Run the standard workflow for visualization and clustering
RNA_GAMs <- NormalizeData(RNA_GAMs)
RNA_GAMs <- FindVariableFeatures(RNA_GAMs, selection.method = "vst", nfeatures = 2000)
RNA_GAMs <- ScaleData(RNA_GAMs, verbose = FALSE)
RNA_GAMs <- RunPCA(RNA_GAMs, npcs = 30, verbose = FALSE)
RNA_GAMs <- RunHarmony(RNA_GAMs, c('patient'))
RNA_GAMs <- RunUMAP(RNA_GAMs, reduction = "harmony", dims = 1:30)
RNA_GAMs <- FindNeighbors(RNA_GAMs, reduction = "harmony", dims = 1:30)
RNA_GAMs <- FindClusters(RNA_GAMs, resolution = 0.3)

RNA.markers <- FindAllMarkers(RNA_all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

GAM.markers.08 <- FindAllMarkers(RNA_GAMs_new, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
GAM.markers.08 %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC) -> GAM.markers.table.08

features <- unique(GAM.markers.table.02$gene)

p1 <- DotPlot(RNA_GAMs, features = gene_list_FigR, group.by = 'type') + RotatedAxis()


saveRDS(RNA_GAMs,file="./RNA_GAMs.rds")


#perform analysis similar to friedrich
micr_list <- list(RNA_20201012,RNA_20201013,RNA_20201116,RNA_20201117,RNA_0218,RNA_0308,RNA_0330,RNA_0405,RNA_0603,RNA_1119,RNA_117am,RNA_0413,RNA_0628,RNA_0416,RNA_0428)

  #sctransform
  for (i in 1:length(micr_list)) {
    micr_list[[i]] <- SCTransform(micr_list[[i]], 
                                  variable.features.n = 10000,
                                  verbose = T,
                                  return.only.var.genes = F,
                                  vars.to.regress = c("percent.mt")
    )
  }




anchor_features <- grep("^(FOS|JUN|RP|ZFP36|EGR|MALAT1|XIST|HSP|MT-|HIST)", SelectIntegrationFeatures(micr_list, nfeatures = 5000), invert = T, value = T)
  
  micr_list <- PrepSCTIntegration(object.list = micr_list, 
                                  anchor.features = anchor_features, 
                                  verbose = T)
  
  immune.anchors <- FindIntegrationAnchors(object.list = micr_list, 
                                           normalization.method = "SCT",
                                           dims = 1:20,
                                           anchor.features = anchor_features)
  
saveRDS(immune.anchors,file = "./immune_anchors_int.rds")
  immune.combined <- IntegrateData(anchorset = immune.anchors, 
                                   normalization.method = "SCT", 
                                   dims = 1:20,features = anchor_features[1:2000])

 