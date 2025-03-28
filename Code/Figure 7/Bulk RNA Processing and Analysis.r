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

bam1 <- '../Project_16770_B/20241121_DP_1_IGO_16770_B_1.bam'
bam2 <- '../Project_16770_B/20241121_DP_2_IGO_16770_B_2.bam'
bam3 <- '../Project_16770_B/20241121_DP_3_IGO_16770_B_3.bam'
bam4 <- '../Project_16770_B/20241121_VEH_1_IGO_16770_B_4.bam'
bam5 <- '../Project_16770_B/20241121_VEH_2_IGO_16770_B_5.bam'
bam6 <- '../Project_16770_B/20241121_VEH_3_IGO_16770_B_6.bam'

bam_list <- c(bam1, bam2, bam3, bam4, bam5, bam6)

gtffile <- '../BT142_culture/Homo_sapiens.GRCh38.112.gtf.gz'

fc <- featureCounts(files=bam_list, 
                    annot.ext=gtffile, 
                    isGTFAnnotationFile=TRUE,
                    isPairedEnd=TRUE)

batch_counts <- (fc$counts)

write.csv(batch_counts, '../BT142_culture/FeatureCounts_allsamples.csv')
batch_counts <- read.csv('../BT142_culture/FeatureCounts_allsamples.csv', row.names=1)

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
                      condition = c('DP', 'DP', 'DP', 'Veh', 'Veh', 'Veh'),
                      replicate = c('DP1', 'DP2', 'DP3', 'Veh1', 'Veh2', 'Veh3')
                     )

row.names(coldata) <- coldata$samples

coldata$condition <- factor(coldata$condition)

coldata$replicate <- factor(coldata$replicate)


dds <- DESeqDataSetFromMatrix(countData = all_cts,
                              colData = coldata,
                              design = ~condition)

dds$condition <- relevel(dds$condition, ref = "Veh")
dds <- DESeq(dds)
res1 <- results(dds)

write.csv(res1, './DP_vs_Veh_DEG.csv')

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

#################### Volcano Plot (Fig7E) ####################
MES = c('HILPDA','ADM','DDIT3','NDRG1','HERPUD1','DNAJB9','TRIB3','ENO2','AKAP12','SQSTM1','MT1X','ATF3','NAMPT','NRN1','SLC2A1','BNIP3','LGALS3','INSIG2','IGFBP3','PPP1R15A','VIM','PLOD2','GBE1','SLC2A3','FTL','WARS','ERO1L','XPOT','HSPA5','GDF15','ANXA2','EPAS1','LDHA','P4HA1','SERTAD1','PFKP','PGK1','EGLN3','SLC6A6','CA9','BNIP3L','RPL21','TRAM1','UFM1','ASNS','GOLT1B','ANGPTL4','SLC39A14','CDKN1A','HSPA9','CHI3L1','ANXA2','ANXA1','CD44','VIM','MT2A','C1S','NAMPT','EFEMP1','C1R','SOD2','IFITM3','TIMP1','SPP1','A2M','S100A11','MT1X','S100A10','FN1','LGALS1','S100A16','CLIC1','MGST1','RCAN1','TAGLN2','NPC2','SERPING1','C8orf4','EMP1','APOE','CTSB','C3','LGALS3','MT1E','EMP3','SERPINA3','ACTN1','PRDX6','IGFBP7','SERPINE1','PLP2','MGP','CLIC4','GFPT2','GSN','NNMT','TUBA1C','GJA1','TNFRSF1A','WWTR1')
AC = c('CST3','S100B','SLC1A3','HEPN1','HOPX','MT3','SPARCL1','MLC1','GFAP','FABP7','BCAN','PON2','METTL7B','SPARC','GATM','RAMP1','PMP2','AQP4','DBI','EDNRB','PTPRZ1','CLU','PMP22','ATP1A2','S100A16','HEY1','PCDHGC3','TTYH1','NDRG2','PRCP','ATP1B2','AGT','PLTP','GPM6B','F3','RAB31','PPAP2B','ANXA5','TSPAN7')
OPC = c('BCAN','PLP1','GPR17','FIBIN','LHFPL3','OLIG1','PSAT1','SCRG1','OMG','APOD','SIRT2','TNR','THY1','PHYHIPL','SOX2-OT','NKAIN4','LPPR1','PTPRZ1','VCAN','DBI','PMP2','CNP','TNS3','LIMA1','CA10','PCDHGC3','CNTN1','SCD5','P2RX7','CADM2','TTYH1','FGF12','TMEM206','NEU4','FXYD6','RNF13','RTKN','GPM6B','LMF1','ALCAM','PGRMC1','HRASLS','BCAS1','RAB31','PLLP','FABP5','NLGN3','SERINC5','EPB41L2','GPR37L1')
NPC = c('DLL3','DLL1','SOX4','TUBB3','HES6','TAGLN3','NEU4','MARCKSL1','CD24','STMN1','TCF12','BEX1','OLIG1','MAP2','FXYD6','PTPRS','MLLT11','NPPA','BCAN','MEST','ASCL1','BTG2','DCX','NXPH1','HN1','PFN2','SCG3','MYT1','CHD7','GPR56','TUBA1A','PCBP4','ETV1','SHD','TNR','AMOTL2','DBN1','HIP1','ABAT','ELAVL4','LMF1','GRIK2','SERINC5','TSPAN13','ELMO1','GLCCI1','SEZ6L','LRRN1','SEZ6','SOX11','STMN2','CD24','RND3','HMP19','TUBB3','MIAT','DCX','NSG1','ELAVL4','MLLT11','DLX6-AS1','SOX11','NREP','FNBP1L','TAGLN3','STMN4','DLX5','SOX4','MAP1B','RBFOX2','IGFBPL1','STMN1','HN1','TMEM161B-AS1','DPYSL3','SEPT3','PKIA','ATP1B1','DYNC1I1','CD200','SNAP25','PAK3','NDRG4','KIF5A','UCHL1','ENO2','KIF5C','DDAH2','TUBB2A','LBH','LOC150568','TCF4','GNG3','NFIB','DPYSL5','CRABP1','DBN1','NFIX','CEP170','BLCAP')

keyvals.colors <- ifelse(
    res1$symbol %in% MES, '#ef476f',
      ifelse(res1$symbol %in% AC, '#26547c',
        ifelse(res1$symbol %in% OPC, '#ffd166',
          ifelse(res1$symbol %in% NPC, '#06d6a0',
        'grey'))))
keyvals.colors[is.na(keyvals.colors)] <- 'grey'
names(keyvals.colors)[keyvals.colors == '#06d6a0'] <- 'NPC-like'
names(keyvals.colors)[keyvals.colors == '#ffd166'] <- 'OPC-like'
names(keyvals.colors)[keyvals.colors == '#26547c'] <- 'AC-like'
names(keyvals.colors)[keyvals.colors == '#ef476f'] <- 'MES-like'

options(repr.plot.width=10, repr.plot.height=10)
EnhancedVolcano(as.data.frame(res1), lab = res1$symbol, 
                x='log2FoldChange', 
                y='pvalue', 
                pCutoff = 0.05, 
                title = 'DP vs Veh',
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                selectLab = c('CD44', 'VIM', 'CHI3L1', 'ADM', 'GFAP', 'OLIG1', 'SOX4', 'DCX', 'DLL3', 'CD24'),
                drawConnectors = TRUE,
                colCustom = keyvals.colors,
                colAlpha=0.8,
                subtitle = NULL,
                ) -> p2
p2 <- p2 + theme(aspect.ratio = 1)
p2

#################### Module Scoring (Fig7F) ####################
library("FerrenaBulkRNAseq")

mes.df <- as.data.frame(modulescore(df.final, MES, numbins = 2, numcontrolgenesperbin = 100))
ac.df <- as.data.frame(modulescore(df.final, AC, numbins = 2, numcontrolgenesperbin = 100))
npc.df <- as.data.frame(modulescore(df.final, NPC, numbins = 2, numcontrolgenesperbin = 100))
opc.df <- as.data.frame(modulescore(df.final, OPC, numbins = 2, numcontrolgenesperbin = 100))

scores.df <- mes.df
scores.df$ac <- ac.df$`modulescore(df.final, AC, numbins = 2, numcontrolgenesperbin = 100)`
scores.df$opc <- opc.df$`modulescore(df.final, OPC, numbins = 2, numcontrolgenesperbin = 100)`
scores.df$npc <- npc.df$`modulescore(df.final, NPC, numbins = 2, numcontrolgenesperbin = 100)`

colnames(scores.df) <- c('mes', 'ac', 'opc', 'npc')
rownames(scores.df) <- c('DP1', 'DP2', 'DP3', 'Veh1', 'Veh2', 'Veh3')
