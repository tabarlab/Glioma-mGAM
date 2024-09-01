library(edgeR)
library(limma)
library(RColorBrewer)
library(tidyverse)
library(matrixStats)
library("DESeq2")
library(dplyr)
library(AUCell)
library("Rsubread")
library(org.Hs.eg.db)
library(annotate)

# mono mGAMs
bam1 <- '../bulk_monoGAMs/N_mono1_IGO_15954_1.bam'
bam2 <- '../bulk_monoGAMs/H_mono1_IGO_15954_2.bam'
bam3 <- '../bulk_monoGAMs/NC_mono1_IGO_15954_3.bam'
bam4 <- '../bulk_monoGAMs/HC_mono1_IGO_15954_4.bam'
bam5 <- '../bulk_monoGAMs/NC_mono2_IGO_15954_7.bam'
bam6 <- '../bulk_monoGAMs/HC_mono2_IGO_15954_8.bam'
bam7 <- '../bulk_monoGAMs/NC_mono3_IGO_15954_11.bam'
bam8 <- '../bulk_monoGAMs/HC_mono3_IGO_15954_12.bam'

bam_list <- c(bam1, bam2, bam3, bam4, bam5, bam6, bam7, bam8)

gtffile <- '../bulk_monoGAMs/Homo_sapiens.GRCh38.112.gtf.gz'

fc <- featureCounts(files=bam_list, 
                    annot.ext=gtffile, 
                    isGTFAnnotationFile=TRUE,
                    isPairedEnd=TRUE)

batch_counts <- (fc$counts)

write.csv(batch_counts, '../bulk_monoGAMs/FeatureCounts_allsamples.csv')
batch_counts <- read.csv('../bulk_monoGAMs/FeatureCounts_allsamples.csv', row.names=1)

batch_counts <- as.data.frame(batch_counts)
batch_counts$ensembl_gene_id <- rownames(batch_counts)

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(batch_counts)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)

res_merge <- merge(batch_counts,G_list,by='ensembl_gene_id')
res_merge$ensembl_gene_id <- NULL
res_merge <- aggregate(res_merge[,-9], list(Symbol=res_merge[,9]), FUN = sum)
rownames(res_merge) <- res_merge$Symbol
res_merge$Symbol <- NULL
head(res_merge)

cpmdata <- cpm(res_merge)
thresh <- cpmdata > 0.5
keep <- rowSums(thresh) >=2
counts.keep <- cpmdata[keep, ]

coldata <- data.frame(samples = colnames(all_cts),
                      condition = c('Normoxia_mono', 'Hypoxia_mono','Normoxia_cocult', 'Hypoxia_cocult','Normoxia_cocult', 'Hypoxia_cocult','Normoxia_cocult', 'Hypoxia_cocult'),
                      hypoxia = c('Normoxia', 'Hypoxia','Normoxia', 'Hypoxia','Normoxia', 'Hypoxia','Normoxia', 'Hypoxia'),
                      culture = c('Mono', 'Mono', 'Co', 'Co', 'Co', 'Co', 'Co', 'Co'),
                      patient = c('m1','m1','m1','m1','m2','m2','m3','m3'),
                      batch = c('2','2','2','2','2','2','2','2')
                     )

row.names(coldata) <- coldata$samples
coldata$condition <- factor(coldata$condition)
coldata$hypoxia <- factor(coldata$hypoxia)
coldata$culture <- factor(coldata$culture)
coldata$patient <- factor(coldata$patient)
coldata$batch <- factor(coldata$batch)

dds <- DESeqDataSetFromMatrix(countData = all_cts,
                              colData = coldata,
                              design = ~patient + condition)

dds$condition <- relevel(dds$condition, ref = "Normoxia_mono")
dds <- DESeq(dds)

res1 <- results(dds, name="condition_Hypoxia_mono_vs_Normoxia_mono")
res2 <- results(dds, name="condition_Normoxia_cocult_vs_Normoxia_mono")
res3 <- results(dds, name="condition_Hypoxia_cocult_vs_Normoxia_mono")

write.csv(res1, './mono_mGAMs/H_vs_N.csv')
write.csv(res2, './mono_mGAMs/NC_vs_N.csv')
write.csv(res3, './mono_mGAMs/HC_vs_N.csv')

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

#################### PCA Plots (Fig6D) ####################
options(repr.plot.width=10, repr.plot.height=5)
pcaData <- plotPCA(vsd, intgroup=c("culture", 'hypoxia', 'patient'), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p1 <- ggplot(pcaData, aes(PC1, PC2, color=culture, shape=hypoxia))+
        geom_point(size=4)+
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance"))+
        scale_color_manual(values=c('#d62828', '#003049'))  + 
        theme_classic() + 
        theme(aspect.ratio=1)
p2 <- ggplot(pcaData, aes(PC1, PC2, color=patient, shape=hypoxia))+
        geom_point(size=4)+
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance"))+
        scale_color_manual(values=c('green', 'blue', 'pink'))  + 
        theme_classic() + 
        theme(aspect.ratio=1)

options(repr.plot.width=10, repr.plot.height=5)
p3 <- (plotPCA(vsd, intgroup=c("culture")) + scale_fill_manual(values=c('black', 'red'))  + theme_classic() + theme(aspect.ratio=1)) | 
        (plotPCA(vsd, intgroup=c("hypoxia")) + theme_classic() + theme(aspect.ratio=1))

