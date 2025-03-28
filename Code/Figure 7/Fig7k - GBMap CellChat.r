library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(ComplexHeatmap)
library(patchwork)
library(ggrepel)
library(harmony)
set.seed(1234)

gbmap <- readRDS('../../azimuth_GBmap.rds')

# Add mGAM label to GAMs
FOSL2_outs <- c("VIM", "NPC2","LAPTM5","CTSD","HMOX1","MCL1","CD163","CXCL8","PLIN2","FOSB","CXCR4","S100A10","S100A6","ATP1B3","TFRC","IL1B","STAB1","ANXA1","SQSTM1","TGFBI","CAST","ARL4C","IQGAP1","ANXA2","RGCC","SERPINE1","SLC16A10","VEGFA","CD44","LY6E","HK2","METRNL","ATP6V1G1","CD9","MDM2","FOSL2","SOCS3","CD93","THBD","NRP1","NDRG1","ANKH","UPP1","LMNA","DHRS3","ATXN1","MYO5A","RPN1","MIR24-2","VKORC1","ASPH","FAM20C","PDE8A","SLC6A6")

gbmap <- AddModuleScore(
  object = gbmap,
  features = list(FOSL2_outs),
  ctrl = 100,
  name = 'FOSL2_module_v', 
  search = TRUE)

gbmap$FOSL2_mod_positive <- as.integer(gbmap$FOSL2_module_v1 >= 0.4)
gbmap$mGAM <- as.integer(gbmap$FOSL2_mod_positive == 1 & gbmap$annotation_level_3 %in% c('Mono', 'TAM-BDM', 'TAM-MG'))
gbmap$new_annot <- gbmap$annotation_level_4
gbmap$new_annot <- as.character(gbmap$new_annot)
gbmap$new_annot[gbmap$new_annot == '1'] <- 'mGAM'

# Run CellChat on new anotations
library(CellChat)

cellchat <- createCellChat(object = gbmap, group.by = "new_annot")

CellChatDB <- CellChatDB.human

## Add custom known interactions then reimport
write.csv(CellChatDB$interaction, './interaction_input_CellChat.csv')
write.csv(CellChatDB$complex, './complex_input_CellChat.csv')
write.csv(CellChatDB$cofactor, './cofactor_input_CellChat.csv')

options(stringsAsFactors = FALSE)  
interaction_input <- read.csv(file = 'interaction_input_CellChat.csv') 
row.names(interaction_input) <- interaction_input[,1] 
CellChatDB$interaction <- interaction_input

# Continue
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)

future::plan("multisession", workers = 4)

cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)
saveRDS(cellchat, './cellchat_GBmap_mGAM_precompute.rds')

cellchat <- computeCommunProb(cellchat)
saveRDS(cellchat, './cellchat_GBmap_mGAM_postcompute_interactive.rds')

cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

options(repr.plot.width=20, repr.plot.height=20)
mat <- cellchat@net$weight
par(mfrow = c(1,2), xpd=TRUE)
for (i in 19:20) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


