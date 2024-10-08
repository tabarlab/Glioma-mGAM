# Script housing main FigR wrapper functions

### Author: Vinay Kartha
### Contact: <vinay_kartha@g.harvard.edu>
### Affiliation: Buenrostro Lab, Department of Stem Cell and Regenerative Biology, Harvard University

library(ggplot2)
library(ggrepel)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(networkD3)
library(GGally)
library(igraph)
library(network)

library(ggrastr)
library(BuenColors)

source("./R/utils.R")



# ---------------------------------------------------------------------------------- FUNCTION DEFINITIONS

runFigR <- function(ATAC.se, # SE of scATAC peak counts. Needed for chromVAR bg peaks etc.
                    dorcK=30, # How many dorc kNNs are we using to pool peaks
                    dorcTab, # peak x DORC connections (should contain indices relative to peaks in ATAC.se)
                    n_bg=50, # No. of background peaks to use for motif Z test
                    genome, # One of mm10, hg19, hg38, with no default
                    dorcMat, # Expect smoothed
                    rnaMat, # Expect smoothed
                    dorcGenes=NULL, # If only running on a subset of genes
                    nCores=1
){
  # Must be matched data
  stopifnot(all.equal(ncol(dorcMat),ncol(rnaMat)))
  
  # Expects "Gene" / "Peak" in dorcTab
  if(!all(c("Peak","Gene") %in% colnames(dorcTab)))
    stop("Expecting fields Peak and Gene in dorcTab data.frame .. see runGenePeakcorr function in BuenRTools")
  
  if(all(grepl("chr",dorcTab$Peak,ignore.case = TRUE))) {
    usePeakNames <- TRUE
    message("Detected peak region names in Peak field")
    
    if(!(all(grepl("chr",rownames(ATAC.se),ignore.case = TRUE))))
      stop("Peak regions provided in dorcTab data.frame but not found as rownames in input SE")
      
    if(!all(dorcTab$Peak %in% rownames(ATAC.se)))
      stop("Found DORC peak region not present in input SE.. make sure DORC calling output corresponds to same input SE as the one provided here ..")
  } else{
    usePeakNames <- FALSE
    message("Assuming peak indices in Peak field")
  # If using index, make sure no indices are outside range of SE
    if(max(dorcTab$Peak) > nrow(ATAC.se))
      stop("Found DORC peak index outside range of input SE.. make sure DORC calling output corresponds to same input SE as the one provided here ..")
  }  
  
  
  if(is.null(dorcGenes)) {
    dorcGenes <- rownames(dorcMat)
  } else {
    cat("Using specified list of dorc genes ..\n")
    if (!(all(dorcGenes %in% rownames(dorcMat)))) {
      cat("One or more of the gene names supplied is not present in the DORC matrix provided: \n")
      cat(dorcGenes[!dorcGenes %in% rownames(dorcMat)], sep = ", ")
      cat("\n")
      stop()
    }
  }
  
  DORC.knn <- FNN::get.knn(data = t(scale(Matrix::t(dorcMat))),k = dorcK)$nn.index # Scaled
  rownames(DORC.knn) <- rownames(dorcMat)
  
  if (is.null(rowData(ATAC.se)$bias)) {
    if (genome %in% "hg19") 
      myGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    if (genome %in% "mm10") 
      myGenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    if (genome %in% "hg38") 
      myGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    ATAC.se <- chromVAR::addGCBias(ATAC.se, genome = myGenome)
  }
  
  if(grepl("hg",genome)){
    #pwm <- chromVARmotifs::human_pwms_v2
    pwm <- readRDS("./data/cisBP_human_pfms_2021.rds")
  } else {
    #pwm <- BuenRTools::mouse_pwms_v3
    pwm <- readRDS("./data/cisBP_mouse_pfms_2021.rds")
  }
  
  # Old motif naming convention
  if(all(grepl("_",names(pwm),fixed = TRUE)))
     names(pwm) <- extractTFNames(names(pwm))
  
  message("Removing genes with 0 expression across cells ..\n")
  rnaMat <- rnaMat[Matrix::rowSums(rnaMat)!=0,]
  myGeneNames <- gsub(x = rownames(rnaMat),pattern = "-",replacement = "") # NKX2-1 to NKX21 (e.g.)
  rownames(rnaMat) <- myGeneNames
  
  # Only non-zero expression TFs (also found in rnaMat)
  motifsToKeep <- intersect(names(pwm),myGeneNames)
    
  # This has to be done on the full SE (same peakset used as input to dorc calling)
  cat("Getting peak x motif matches ..\n")
  motif_ix <- motifmatchr::matchMotifs(subject = ATAC.se,pwms = pwm[motifsToKeep],genome=genome)
  
  # Keep TFs with some peak x motif match
  motif_ix <- motif_ix[,Matrix::colSums(assay(motif_ix))!=0] 
  
  cat("Determining background peaks ..\n")
  cat("Using ", n_bg, " iterations ..\n\n")
  if(any(Matrix::rowSums(assay(ATAC.se))==0)){
    ATAC.mat <- assay(ATAC.se)
    ATAC.mat <- cbind(ATAC.mat,1)
    ATAC.se.new <- SummarizedExperiment::SummarizedExperiment(assays = list(counts=ATAC.mat),rowRanges = granges(ATAC.se))
    set.seed(123)
    bg <- chromVAR::getBackgroundPeaks(ATAC.se.new, niterations = n_bg)
  } else {
    set.seed(123)
    bg <- chromVAR::getBackgroundPeaks(ATAC.se, niterations = n_bg)
  }
  
  # For each DORC, do motif enrichment among dorc sig Peaks, and correlation of DORC accessibility (smoothed) to TF RNA levels
  
  cat("Testing ",length(motifsToKeep)," TFs\n")
  cat("Testing ",nrow(dorcMat)," DORCs\n")
  library(doParallel)
  if(nCores > 1)
    message("Running FigR using ",nCores," cores ..\n")
  opts <- list()
  pb <- txtProgressBar(min = 0, max = length(dorcGenes), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  time_elapsed <- Sys.time()
  cl <- parallel::makeCluster(nCores)
  doSNOW::registerDoSNOW(cl)
  mZtest.list <- foreach(g=dorcGenes,
                         .options.snow = opts, 
                         .packages = c("dplyr","Matrix","chromVAR","Signac","Seurat","motifmatchr","SummarizedExperiment")) %dopar%   {
                           
                           source("./R/utils.R")
                           
                           

                           
                           
                           
                           
                           # Take peaks associated with gene and its k neighbors
                           # Pool and use union for motif enrichment
                           DORCNNpeaks <- unique(dorcTab$Peak[dorcTab$Gene %in% c(g,dorcGenes[DORC.knn[g,]])])
                           
                           if(usePeakNames)
                             DORCNNpeaks <- which(rownames(ATAC.se) %in% DORCNNpeaks) # Convert to index relative to input
                             
                           mZ <- motifPeakZtest(peakSet = DORCNNpeaks,
                                                bgPeaks = bg,
                                                tfMat = assay(motif_ix))
                           mZ <- mZ[,c("gene","z_test")]
                           colnames(mZ)[1] <- "Motif"
                           colnames(mZ)[2] <- "Enrichment.Z"
                           mZ$Enrichment.P <- 2*pnorm(abs(mZ$Enrichment.Z),lower.tail = FALSE) # One-tailed
                           mZ$Enrichment.log10P <- sign(mZ$Enrichment.Z) * -log10(mZ$Enrichment.P)
                           mZ <- cbind("DORC"=g,mZ)
                           # Correlate smoothed dorc with smoothed expression, with spearman
                           corr.r <- cor(dorcMat[g,],t(as.matrix(rnaMat[mZ$Motif,])),method = "spearman")
                           stopifnot(all.equal(colnames(corr.r),mZ$Motif))
                           mZ$Corr <- corr.r[1,] # Correlation coefficient
                           mZ$Corr.Z <- scale(mZ$Corr,center = TRUE,scale = TRUE)[,1] # Z-score among all TF correlations
                           mZ$Corr.P <- 2*pnorm(abs(mZ$Corr.Z),lower.tail = FALSE) # One-tailed
                           mZ$Corr.log10P <- sign(mZ$Corr.Z)*-log10(mZ$Corr.P)
                           return(mZ)
                         }
  cat("Finished!\n")
  cat("Merging results ..\n")
  # Merge and save table for downstream filtering and plotting (network)
  TFenrich.d <- do.call('rbind',mZtest.list)
  dim(TFenrich.d)
  rownames(TFenrich.d) <- NULL
  
  # Make combined score based on multiplication
  # Here, we only sign by corr
  # Since sometimes we lose digit precision (1 - v small number is 1, instead of 0.9999999..)
  # Use Rmpfr, increase precision limits above default (100 here)
  TFenrich.d <- TFenrich.d %>% mutate("Score"=sign(Corr)*as.numeric(-log10(1-(1-Rmpfr::mpfr(Enrichment.P,100))*(1-Rmpfr::mpfr(Corr.P,100)))))
  TFenrich.d$Score[TFenrich.d$Enrichment.Z < 0] <- 0
  TFenrich.d
}


# Bar plot of ranked activator / repressors
rankDrivers <- function(figR.d,
                        myLabels=NULL){
  
    figR.summ <- figR.d %>%  group_by(Motif) %>% 
      summarise(Score=mean(Score)) %>% 
      arrange(Score) %>% 
      mutate(Motif=factor(Motif,levels=as.character(Motif)))
    
    # Top and bottom %ile labels
    figR.summ$Label <- as.character(figR.summ$Motif)
    
    if(is.null(myLabels)){
      # Use quantiles to define what labels are drawn
      figR.summ$Label[figR.summ$Score >= quantile(figR.summ$Score,0.05) & figR.summ$Score <= quantile(figR.summ$Score,0.95)] <- ""
    } else {
      # Only highlight user-specified
      figR.summ$Label[!figR.summ$Label %in% myLabels] <- ""
    }
      
    
    gAll <- ggplot(figR.summ,aes(x=Motif,y=Score,label=Label)) + 
      geom_bar(size=0.1,stat="identity",fill="darkorange",color=NA) + 
      theme_classic() + theme(axis.text.x = element_blank()) + 
      geom_text_repel(size=3,min.segment.length = 0.1,segment.size = 0.2,max.overlaps = 20) + 
      geom_hline(yintercept = 0) + labs(x="TF Motifs",y="Regulation Score")
    
    gAll
}



# Scatter
plotDrivers <- function(figR.d,
                        marker,
                        score.cut=1,
                        label=TRUE){

  if(!marker %in% figR.d$DORC)
    stop("Marker specified is not a valid DORC symbol found in the data.frame")

  d <- figR.d %>% filter(DORC %in% marker) %>% mutate(isSig=ifelse(abs(Score) >= score.cut,"Yes","No"))
  d$Label <- d$Motif
  d$Label[d$isSig %in% "No"] <- ""
  gScatter <- d %>%  ggplot(aes(x=Corr.log10P,y=Enrichment.log10P,color=isSig,label=Label)) + 
    geom_hline(yintercept = 0,color="gray60",linetype="dashed") + 
    geom_vline(xintercept = 0,color="gray60",linetype="dashed") + 
    geom_point(size=0.8,shape=16) + theme_classic() + 
    scale_color_manual(values=c("gray66","firebrick3"))+
    scale_x_continuous(breaks=scales::pretty_breaks()) +
    scale_y_continuous(breaks=scales::pretty_breaks()) + 
    labs(y="Enrichment log10 P",x="Correlation log10 P",title=marker) + 
    ylim(-ceiling(max(abs(d$Enrichment.log10P))), ceiling(max(abs(d$Enrichment.log10P))))+
    xlim(-ceiling(max(abs(d$Corr.log10P))), ceiling(max(abs(d$Corr.log10P))))+
    theme(legend.position = "none",axis.text = element_text(color="black"),
          plot.title = element_text(hjust=0.5,face="italic"),panel.background = element_rect(fill=NA)) + 
    geom_text(hjust=1.1,fontface="italic",color="black",size=3)
  
  gScatter
  
}

plotfigRHeatmap <- function(figR.d,
                            score.cut=1,
                            DORCs=NULL,
                            TFs=NULL,
                            ... # Additional params passed to ComplexHeatmap
){
  
  
  DORCsToKeep <- figR.d %>% filter(abs(Score) >= score.cut) %>% pull(DORC) %>% unique()
  TFsToKeep <- figR.d %>% filter(abs(Score) >= score.cut) %>% pull(Motif) %>% unique()
  
  
  if(!is.null(DORCs)){
    if(!all(DORCs %in% figR.d$DORC))
      stop("One or more DORCs specified is not a valid DORC symbol found in the data.frame")
    DORCsToKeep <- intersect(DORCsToKeep,DORCs)
    TFsToKeep <- figR.d %>% filter(abs(Score) >= score.cut & DORC %in% DORCsToKeep) %>% pull(Motif) %>% unique()
  }
  
  
  if(!is.null(TFs)){
    if(!all(TFs %in% figR.d$Motif))
      stop("One or more TFs specified is not a valid TF symbol found in the data.frame")
    TFsToKeep <- intersect(TFsToKeep,TFs)
    DORCsToKeep <- figR.d %>% filter(abs(Score) >= score.cut & Motif %in% TFsToKeep) %>% pull(DORC) %>% unique()
  }
 
  
  net.d <- figR.d %>% filter(DORC %in% DORCsToKeep & Motif %in% TFsToKeep) %>% 
    reshape2::dcast(DORC ~ Motif) %>% 
    tibble::column_to_rownames("DORC") %>% as.matrix()

  message("Plotting ",nrow(net.d)," DORCs x ",ncol(net.d), "TFs\n")
  
  # Heatmap view
  
  myCols <- colorRamp2(seq(-2,2,length.out = 9),colors = jdb_palette("solar_flare"))
  myHeat <- Heatmap(net.d,
          col=myCols,
          clustering_distance_rows = "pearson",
          clustering_distance_columns = "pearson",
          name="Score",border = TRUE,
          row_names_gp = gpar(fontsize=5,fontface="italic"),...)
  
  myHeat
  
}

plotfigRNetwork <- function(figR.d,
                        score.cut=1,
                        DORCs=NULL,
                        TFs=NULL,
                        weight.edges=FALSE){
# Network view

# Filter
net.dat <-  figR.d %>% filter(abs(Score) >= score.cut)

if(!is.null(DORCs))
  net.dat <- net.dat %>% filter(DORC %in% DORCs)

if(!is.null(TFs))
  net.dat <- net.dat %>% filter(Motif %in% TFs)

net.dat$Motif <- paste0(net.dat$Motif, ".")
net.dat$DORC <- paste0(net.dat$DORC)


dorcs <- data.frame(name = unique(net.dat$DORC), group = "DORC", size = 8)
tfs <- data.frame(name = unique(net.dat$Motif), group = "TF", size = 3)
nodes <- rbind(dorcs,tfs)

edges <- as.data.frame(net.dat)

# Make edges into links (subtract 1 for 0 indexing)
links <- data.frame(source=unlist(lapply(edges$Motif, function(x) {which(nodes$name==x)-1})), 
                    target=unlist(lapply(edges$DORC, function(x) {which(nodes$name==x)-1})), 
                    corr=edges$Corr,
                    enrichment=edges$Enrichment.P)

links$Value <- scales::rescale(edges$Score)*20

# Set of colors you can choose from for TF/DORC nodes
colors <- c("Red", "Orange", "Yellow", "Green", "Blue", "Purple", "Tomato", "Forest Green", "Sky Blue","Gray","Steelblue3","Firebrick2","Brown")
nodeColorMap <- data.frame(color = colors, hex = gplots::col2hex(colors))

getColors <- function(tfColor, dorcColor = NULL) {
  temp <- c(as.character(nodeColorMap[nodeColorMap$color==tfColor,]$hex),
            as.character(nodeColorMap[nodeColorMap$color==dorcColor,]$hex))
  if (is.null(dorcColor)) {
    temp <- temp[1]
  }
  colors <- paste(temp, collapse = '", "')
  colorJS <- paste('d3.scaleOrdinal(["', colors, '"])')
  colorJS
}

forceNetwork(Links = links, 
             Nodes = nodes,
             Source = "target",
             Target = "source",
             NodeID = "name",
             #NodeID="myLabel",
             Group = "group",
             Value = "Value",
             Nodesize = "size",
             #linkWidth = 1.5, # Fixed link weight
             radiusCalculation = "Math.sqrt(d.nodesize)*2",
             arrows = FALSE,
             opacityNoHover = 0.6,
             opacity = 1,
             zoom = TRUE,
             bounded = TRUE,
             charge = -15, 
             fontSize = 13,
             legend = TRUE,
             colourScale = getColors(tfColor = "Tomato",dorcColor =  "Sky Blue"), # TF then DORC
             linkColour = ifelse(links$corr > 0, as.character(nodeColorMap[nodeColorMap$color=="Forest Green",]$hex),
                                 as.character(nodeColorMap[nodeColorMap$color=="Purple",]$hex)))

}