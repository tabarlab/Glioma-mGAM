# Script housing dependency and wrapper function to perform optmatching between scATAC and scRNA cells
# Assumes CCA has been run on both assays (together) already (see Methods of associated manuscript for details)

### Author: Vinay Kartha, Yan Hu
### Contact: <vinay_kartha@g.harvard.edu>
### Affiliation: Buenrostro Lab, Department of Stem Cell and Regenerative Biology, Harvard University



library(optmatch)
library(Matrix)
library(FNN)
library(dplyr)

# Function to plot resulting pairs
# Assumes condordant vectors for ATAC and RNA labels
# Also requires a data frame of UMAP coordinates with ATAC RNA labels as rownames
plotPairs <- function(ATAC,
                      RNA,
                      umap.df, # Assumes rownames are cell barcodes
                      max.show=300,
                      seed=123,
                      pairPoiintSize=1
                      ){

  cat("Assuming concordance in pairs ..\n")

  stopifnot(!is.null(rownames(umap.df)))
  stopifnot(all.unique(rownames(umap.df)))
  stopifnot(length(ATAC)==length(RNA))

  if(length(ATAC) > max.show){
    cat("Only highlighting",max.show,"pairs at random..\n")
    set.seed(seed)
    subsample <- sample(1:length(ATAC),max.show,replace = FALSE)
    ATAC <- ATAC[subsample]
    RNA <- RNA[subsample]
  }

  d1 <- umap.df[ATAC,1:2]
  d2 <- umap.df[RNA,1:2]

  d1$pair <- factor(1:nrow(d1))
  d2$pair <- d1$pair

  d1$Group <- "ATAC"
  d2$Group <- "RNA"

  d <- rbind(d1,d2)


  require(ggrastr)
  ggplot(umap.df, aes(x = UMAP1, y = UMAP2)) +
    geom_point_rast(color="gray95",size=1.5)+
    geom_point(data=d,aes(UMAP1,UMAP2,color=Group),alpha=0.8,size=pairPoiintSize) +
    scale_color_manual(values=c("cadetblue","darkorange"))+
    geom_line(data = d, aes(x = UMAP1, y = UMAP2, group = pair),alpha=0.4,color="darkgray")+
    theme_void()+
    labs(color="")+
    theme(legend.position = "bottom")+
    guides(colour = guide_legend(override.aes = list(size=3)))
}


# Function that gives you an "up-sampled" dataset given two paired datasets
# Draws cells using density-based upsampling from smaller one till it matches larger one

smartUpSample <- function(mat1,
                          mat2,
                          seed=123,
                          dist.use="pearson"){
  if(nrow(mat1) >= nrow(mat2)){
    bigger <- mat1
    smaller <- mat2 } else {
      bigger <- mat2
      smaller <- mat1
    }

  big.cells <- rownames(bigger)
  small.cells <- rownames(smaller)

  if(is.null(big.cells) | is.null(small.cells))
    stop("Both ATAC and RNA CCA matrices must have rownames (cell IDs) ..\n")

  # Perform Greedy pairing for these cells based on pearson correlation (max)
  cat("Performing initial greedy pairing using ",dist.use, "as metric ..\n")
  # Rows are bigger matrix cells, cols are smaller matrix cells
  cca.cor <- cor(t(bigger),t(smaller),method = dist.use)
  cor.max <- apply(cca.cor,1,function(x) names(which.max(x)))

  stopifnot(all.equal(rownames(cca.cor),names(cor.max)))

  pairing.greedy <- data.frame("bigger"=names(cor.max),"smaller"=as.character(cor.max),stringsAsFactors = FALSE)

  rm(cca.cor)
  gc()

  # E.g. For each RNA cell, see how many ATAC cells it paired to
  numBiggerpairs <- c()
  for(i in small.cells){
    numBiggerpairs <- c(numBiggerpairs,sum(pairing.greedy$smaller %in% i))
  }


  # Smooth the above number based on RNA CCA knn
  cat("Getting kNN of mapped neigbors ..\n")
  smaller.knn <- FNN::get.knn(smaller,k=200)$nn.index
  rownames(smaller.knn) <- small.cells

  numBiggerpairssmoothed <-  unlist(lapply(1:length(numBiggerpairs),function(i){
    mean(numBiggerpairs[smaller.knn[i,]])
  }))

  samplingDensity <- numBiggerpairssmoothed/length(big.cells)

  cat("Up-sampling smaller dataset ..\n")
  set.seed(seed)
  small.upCells <- sample(small.cells,size = length(big.cells),replace = TRUE,prob = samplingDensity)
  cat("Finished\n")


  smaller.up <- smaller[small.upCells,]

  stopifnot(nrow(smaller.up)==nrow(bigger))

  cat("Both datasets now have ",nrow(bigger), " cells\n")

  cat("% smaller dataset coverage up-on up-sampling: ",round((length(unique(small.upCells))/length(small.cells)) *100,2),"\n")

  if(length(intersect(rownames(smaller),rownames(mat1))) > 0) {
    list(smaller.up,bigger)
  } else {
    list(bigger,smaller.up)
  }

}

# Function that evenly (randomly) samples larger dataset to smaller one
# Chunks up and returns list of chunked CCA ATAC/RNA components

chunk_CCA <- function(CCA_1,
                      CCA_2,
                      chunkSize=NULL,
                      seed=123){
  if(nrow(CCA_1) > nrow(CCA_2)){
    bigger <- CCA_1
    smaller <- CCA_2
  } else {
    bigger <- CCA_2
    smaller <- CCA_1
  }

  if(is.null(rownames(bigger)) | is.null(rownames(smaller)))
    stop("CCA matrices must have valid cell barcode names as rownames ..\n")

  # Shuffle them to make sure there's no bias in chunking
  set.seed(seed)
  bigger <- bigger[sample(1:nrow(bigger),nrow(bigger),replace = FALSE),]
  smaller <- smaller[sample(1:nrow(smaller),nrow(smaller),replace = FALSE),]

  cat("Number of cells in bigger dataset: ",nrow(bigger),"\n")
  cat("Number of cells in smaller dataset: ",nrow(smaller),"\n")
  cat("Difference in cell count between 2 datasets: ",nrow(bigger) - nrow(smaller),"\n")

  if(is.null(chunkSize))
     chunkSize <- nrow(smaller)

  cat("Chunking larger dataset to match smaller datset ..\n")
  cat("Chunk size n = ",chunkSize," cells\n")

  frac <- nrow(bigger) %/% chunkSize # Number of iterations of smaller in bigger
  remain <- nrow(bigger) %% chunkSize # Extra amount we have to sample after even chunking

  cat("Total number of chunks: ",length(1:frac) + sum(remain > 0),"\n")

  bigCells <- c()
  smallCells <- c()

  set.seed(seed)
  chunkList <- list()

  for(i in 1:frac){
    # What has already not been sampled
    bigCellsBag <- rownames(bigger)[!rownames(bigger) %in% bigCells]
    smallCellsBag <- rownames(smaller) # Sample from all cells

    # Fetch chunk for each dataset
    bigCellsChunk <- sample(bigCellsBag,chunkSize,replace=FALSE)
    smallCellsChunk <- sample(smallCellsBag,chunkSize,replace=FALSE)

    # Update wha thas been sampled already
    bigCells <- c(bigCells,bigCellsChunk)
    smallCells <- c(smallCells,smallCellsChunk)

    # List of bigger and smaller CCA PCs for chunk
    # NOTE: Each chunk is always returned as Bigger 1st, smaller 2nd (ATAC, then RNA for us)
    chunkList[[i]] <- list(bigger[bigCellsChunk,],
                          smaller[smallCellsChunk,])
  }




  # Add extra cells for last uneven chunk (if any)
  if(remain > 0){
  chunkList[[i+1]] <- list(bigger[rownames(bigger)[!rownames(bigger) %in% bigCells],],
                           smaller[sample(rownames(smaller),remain,replace=FALSE),])
  }
  chunkList
}


umap_knn_graph <-function(x, #Input data on which we compute the knn graph
                          umap_dim = 5, # Dimensionality of the embedded UMAP space
                          umap_n_nerighbors = 30, # Number of neighbors for computing the UMAP
                          k = 5,
                          seed){# Number of neighbors for computing the KNN graph{

  library(uwot)

  # UMAP embedding of data
  set.seed(seed)
  umap_embedding <- uwot::umap(x, n_neighbors = umap_n_nerighbors, n_components = umap_dim)

  # Get KNN for all cells
  all_knn <- get.knn(umap_embedding, k = k)
  n_data <- dim(x)[1]

  # Construct adjacency matrix of knn graph
  knn_adj_mtx <- array(0, dim = c(n_data, n_data))
  for (i in 1:n_data) knn_adj_mtx[i,all_knn$nn.index[i,]] <- 1

  # We make the graph undirected (make the adjacency matrix a symmetric matrix)
  knn_adj_mtx <- knn_adj_mtx + t(knn_adj_mtx)

  # Initialize a graph object using the adjacency matrix
  knn_graph <- graph_from_adjacency_matrix(knn_adj_mtx, mode = "undirected")
  rm(knn_adj_mtx)

  knn_graph
}


# This function takes in the output of the fullmatch function and sort the results into a list of ATAC-RNA pairs
get_pair_list <- function(cell_matches, # The output object of fullmatch or pairmatch in the package optmatch
                          ATAC_barcodes, # Barcode of ATAC cells, must match the order of IDs in cell_matches.
                          #e.g. ATAC_1 in cell matches should correspond to the first element in ATAC_barcodes
                          RNA_barcodes # Barcode of RNA cells, must match the order of IDs in cell_matches.
                          #e.g. RNA_1 in cell matches should correspond to the first element in RNA_barcodes
){

  cell_matches <- sort(cell_matches) # Sort to get ATAC, RNA tuples

  if(length(cell_matches)==0)
    stop("Matches could not be found .. Perhaps try adjusting the constraints to allow optimal matching to be solved?\n")

  if(any(is.na(cell_matches)))
    warning("NA pairs exist ..\n")

  # Currently the result list contain cells paired to multiple other cells
  # If a cell is paired to k cells, we duplicate this cell k times to make k pairs
  # Thus we generate a new list consisting pairs of 1 ATAC - 1 RNA cells
  # We also make sure that in a pair, the first ID is ATAC and the second ID is RNA
  matches <- character()
  pair_ids <- unique(unname(cell_matches))
  for (pair_id in 1:length(pair_ids)){

    new_match <- names(cell_matches[unname(cell_matches)== pair_ids[pair_id]])
    new_match_ATAC <- new_match[splitAndFetch(new_match,"_",1) == "ATAC"]
    new_match_RNA <- new_match[splitAndFetch(new_match,"_",1) == "RNA" ]
    new_match <- vector()
    if (length(new_match_ATAC) > length(new_match_RNA)){
      new_match[seq(1, 2*length(new_match_ATAC), 2)] <- new_match_ATAC
      new_match[seq(2, 2*length(new_match_ATAC), 2)] <- rep(new_match_RNA, length(new_match_ATAC))
    } else {
      new_match[seq(1, 2*length(new_match_RNA), 2)] <- rep(new_match_ATAC, length(new_match_RNA))
      new_match[seq(2, 2*length(new_match_RNA), 2)] <- new_match_RNA
    }
    matches <- c(matches, new_match)
  }

  # Make sure pair groupings (factors) are adjacent
  #stopifnot(all.equal(cell_matches[seq(1,length(cell_matches),2)],cell_matches[seq(2,length(cell_matches),2)],check.attributes=FALSE))

  ATAC_IDs <- matches[seq(1,length(matches),2)] # Names of ATAC cells
  RNA_IDs <- matches[seq(2,length(matches),2)] # Names of RNA cells

  # Checking if ATAC and RNA tuples are concordant
  stopifnot(all(splitAndFetch(ATAC_IDs,"_",1) %in% "ATAC"))
  stopifnot(all(splitAndFetch(RNA_IDs,"_",1) %in% "RNA"))

  # This is just to make sure 1-1, with the names we gave
  # (can still comprise actual doublets from upsampling if any)
  #stopifnot(all.unique(ATACp) & all.unique(RNAp))

  # Get corresponding index relative to input matrix order
  ATAC_inds <- as.numeric(splitAndFetch(ATAC_IDs,"_",2))
  RNA_inds <- as.numeric(splitAndFetch(RNA_IDs,"_",2))

  matches_mat <- matrix(c(ATAC_inds,RNA_inds),ncol=2,byrow = FALSE) # ATAC col1, RNA col2

  cat("Assembling pair list ..\n")
  # Make data frame of matches

  pair_df <- data.frame("ATAC"=ATAC_inds,
                        "RNA"=RNA_inds)

  pair_df <- pair_df %>% arrange(ATAC)

  # Convert to labels if they exist
  pair_df$ATAC <- ATAC_barcodes[pair_df$ATAC] # ATAC cell labels
  pair_df$RNA <- RNA_barcodes[pair_df$RNA] # RNA cell labels

  pair_df
}


# Main matching function.
# By default, will return a data frame of pairs
# Uses names of barcodes if present in input PC matrices (rownames)
cell_pairing <- function(ATACpcs, # Input ATAC single cell PCs obtained from Running CCA dimensionality reduction. Shape: (n_ATAC_cells, n_PCs)
                         RNApcs,# Input ATAC single cell PCs obtained from Running CCA dimensionality reduction. Shape: (n_RNA_cells, n_PCs)
                         mode = "geodesic",# Mode of pairing. One of "geodesic", "greedy", "optmatch"
                         tol = 0.0001, # tol times the number of subjects to be matched specifies the extent to which fullmatch's
                         #output is permitted to differ from an optimal solution to the original problem,
                         # for details on this tol parameter, see the documentation of optmatch. Look for function fullmatch
                         # If mode is specified to be "geodesic", then we need to provide the below parameters.
                         search_range = 0.2, # This determines the size of the search knn. search_range * total number of cells = size of knn.
                         max_multimatch = 5, # Maximum number of cells allowed to be matched to each cell.
                         umap_knn_k = 30,
                         min_subgraph_size=50, # Minimum number of cells (ATAC/RNA each) needed for pairing in a given subgraph. Will skip subgraphs with fewer than these cells
                         cca_umap_df = NULL, # Optional umap of all (ATAC + RNA cells) to visualize UMAP of cells from each chunk being paired, colored by subgraph
                         seed=123, # Random seed for subgraph UMAP and for down-sampling if subgraphs are imbalanced
                         # If mode is specified to be "optmatch", then we need to provide the below parameters.
                         doPermutation=FALSE, # Test FDR on pair?
                         returnDistMat=FALSE, # Return the sparse distance matrix used for matching?
                         k=300, # k-NN parameter used for applying constraints on ATAC-RNA pairs
                         nCores=1,
                         sanityLoop=FALSE,
                         forceSymmetry=TRUE # ATAC and RNA have to be same dimensions
){

  library(igraph)
  library(FNN)
  library(pracma)
  library(uwot)

  # For now we pool ATAC cells and RNA cells to delineate a manifold
  all_pcs <- rbind(ATACpcs, RNApcs)
  n_cells <- dim(all_pcs)[1]

  if (mode == "geodesic"){
    cat("\n\n")
    cat("Constructing KNN graph for computing geodesic distance ..\n")
    # Get UMAP-based KNN graph
    knn_graph <- umap_knn_graph(all_pcs,k=umap_knn_k,seed = seed)

    cat("Computing graph-based geodesic distance ..\n")
    # Compute shortest paths between all cell pairs. We checked that this returns a symmetric matrix
    shortest_paths <- shortest.paths(knn_graph)

    # Find connected subgraphs
    subgraphs <- clusters(knn_graph)
    cat("# KNN subgraphs detected:\n", length(unique(subgraphs$membership)),"\n")

    if(!is.null(cca_umap_df)){
      cat("Visualizing subgraph for cells based on original CCA UMAP coords provided..\n")
      sub_umap <- cca_umap_df[rownames(all_pcs),]
      sub_umap$subgraph_clust <- subgraphs$membership

      gSub <- ggplot(sub_umap,aes(x=UMAP1,y=UMAP2,color=factor(subgraph_clust))) +
        geom_point(size=0.1) + theme_classic() +
        guides(colour = guide_legend(override.aes = list(size=3)))
      print(gSub)
    }


    all_pairs <- NULL
    
    #k_pairing <- 10 + n_cells * search_range
    #size_threshold <- k_pairing # Subgraphs with fewer nodes than this threshold would be skipped

    cat("Skipping subgraphs with either ATAC/RNA cells fewer than: ",min_subgraph_size," ..\n")

    # Go through each subgraph
    for (subgraph_ind in unique(subgraphs$membership)){

      cat("Pairing cells for subgraph No.", subgraph_ind,"\n")

      # Retrieve the subgraph
      subgraph_nodes <- subgraphs$membership == subgraph_ind
      knn_subgraph <- induced_subgraph(knn_graph, which(subgraph_nodes))

      # Use down-sampling to make sure in this subgraph the number of ATAC and RNA cells are balanced
      subgraph_cells <- c(rownames(ATACpcs), rownames(RNApcs))[subgraph_nodes]
      n_ATAC <- sum(subgraph_cells %in% rownames(ATACpcs))
      n_RNA <- sum(subgraph_cells %in% rownames(RNApcs))

      cat("Total ATAC cells in subgraph: ",n_ATAC,"\n")
      cat("Total RNA cells in subgraph: ",n_RNA,"\n")

      subgraph_ATAC_pcs <- ATACpcs[subgraph_nodes[1:dim(ATACpcs)[1]],]
      subgraph_RNA_pcs <- RNApcs[subgraph_nodes[(dim(ATACpcs)[1] + 1):n_cells],]

      if (n_ATAC > n_RNA){
        set.seed(seed)
        subgraph_ATAC_pcs <- subgraph_ATAC_pcs[sample(1:n_ATAC, n_RNA, replace = FALSE),]
      } else if (n_ATAC < n_RNA){
        set.seed(seed)
        subgraph_RNA_pcs <- subgraph_RNA_pcs[sample(1:n_RNA, n_ATAC, replace = FALSE),]
      }

      if(is.null(nrow(subgraph_ATAC_pcs)) | is.null(nrow(subgraph_RNA_pcs))){
        message("Down-sampling within subgraph between assays led to 0 cells in one assay..\n")
        message("Skipping current subgraph\n")
        next
      }

      # Subset the original geodesic distance matrix to get the geodesic distance matrix for the subgraph
      subgraph_geodist <- shortest_paths[match(rownames(subgraph_ATAC_pcs), rownames(all_pcs)),
                                         match(rownames(subgraph_RNA_pcs), rownames(all_pcs))]

      subgraph_size <- dim(subgraph_geodist)[1]

      cat("Subgraph size: ",subgraph_size,"\n")

      # TO AVOID MAJOR SUBGRAPH(S) WERE BEING SKIPPED SOMETIMES
      size_threshold <- ceiling((nrow(subgraph_ATAC_pcs) + nrow(subgraph_RNA_pcs)) * search_range)
      k_pairing <- size_threshold

      cat("Search threshold being used: ",k_pairing,"\n")

      if (subgraph_size < size_threshold | subgraph_size < min_subgraph_size) {
        message("Insufficient number of cells in subgraph. Skipping current subgraph")
        next
      }

      # We also calculate euclidean distance matrix
      subgraph_eucdist <- distmat(subgraph_ATAC_pcs, subgraph_RNA_pcs)

      # Find KNN based on geodesic distances.
      print("Constructing KNN based on geodesic distance to reduce search pairing search space")
      geodist_knn <- array(-1, dim=dim(subgraph_geodist))
      for (i in 1:subgraph_size){

        # Find RNA cells in the KNN of each ATAC cell
        geodist_threshold <- sort(subgraph_geodist[i,])[k_pairing]
        knn_ind <- subgraph_geodist[i,] < geodist_threshold
        geodist_knn[i, knn_ind] <- subgraph_eucdist[i, knn_ind]

        # Find ATAC cells in the KNN of each RNA cell
        geodist_threshold <- sort(subgraph_geodist[,i])[k_pairing]
        knn_ind <- subgraph_geodist[,i] < geodist_threshold
        geodist_knn[knn_ind, i] <- subgraph_eucdist[knn_ind, i]

      }

      # For an ATAC-RNA cell pair, if neither of them are in the KNN of the other, we set their distance to be inf
      geodist_knn[geodist_knn < 0] <- Inf

      # Add rownames to the matrix
      rownames(geodist_knn) <- paste0("ATAC_",1:subgraph_size)
      colnames(geodist_knn) <- paste0("RNA_",1:subgraph_size)

      print(paste("Number of cells being paired:",
                  dim(geodist_knn)[1],
                  "ATAC and",
                  dim(geodist_knn)[1],
                  " RNA cells"))

      cat("Determing pairs through optimized bipartite matching ..\n")
      options("optmatch_max_problem_size" = Inf)
      cell_matches <- suppressWarnings(optmatch::fullmatch(optmatch::as.InfinitySparseMatrix(as.matrix(geodist_knn)),
                                                           tol = tol,
                                                           min.controls = 1 / max_multimatch,
                                                           max.controls = max_multimatch))
      pair_list <- get_pair_list(cell_matches,
                                 rownames(subgraph_ATAC_pcs),
                                 rownames(subgraph_RNA_pcs))

      cat("Finished!\n")
      
      # Append the results for this subgraph to the list of all results
      all_pairs <- rbind(all_pairs, pair_list)
    }
  }

  if (mode == "greedy"){

    n_ATAC <- dim(ATACpcs)[1]
    n_RNA <- dim(RNApcs)[1]
    if (n_ATAC >= n_RNA){
      print("Pairing all RNA cells to nearest ATAC cells")
      pair_knn <- get.knnx(data = ATACpcs, query = RNApcs, k = 1)
      ATAC_paired <- rownames(ATACpcs)[pair_knn$nn.index]
      RNA_paired <- rownames(RNApcs)
    } else{
      print("Pairing all ATAC cells to nearest RNA cells")
      pair_knn <- get.knnx(data = RNApcs, query = ATACpcs, k = 1)
      ATAC_paired <- rownames(RNApcs)[pair_knn$nn.index]
      RNA_paired <- rownames(ATACpcs)
    }
    all_pairs <- data.frame(ATAC=ATAC_paired,RNA=RNA_paired,stringsAsFactors=FALSE)
  }

  return(all_pairs)

}

# Wrapper function for running ATAC/RNA geodesic pairing per chunk
chunk_pair_geo <- function(ATACpcs,
                           RNApcs,
                           ...){ # Addition parameters pased to cell_pairing function (see below)

  if(nrow(ATACpcs)!=nrow(RNApcs)){
    cat("Running geodesic pairing in after chunking data ..\n")

    # Chunk up data (larger to smaller)
    chunkList <- chunk_CCA(CCA_1 = ATACpcs,CCA_2 = RNApcs,seed = 123)

    # Chunk lists always returned with larger assay first, then smaller (list of two matrices)
    # Match to input ATAC/RNA specification
    if(all(rownames(chunkList[[1]][[1]]) %in% rownames(ATACpcs) & !rownames(chunkList[[1]][[2]]) %in% rownames(ATACpcs))){
      ATAC.i <- 1
      RNA.i <- 2
    } else {
      ATAC.i <- 2
      RNA.i <- 1
    }
    
    
    pairList <- list()
    for(i in 1:length(chunkList)){
      cat("\nChunk # ",i,"\n")
      cat("No. cells in chunk: ",nrow(chunkList[[i]][[1]]),"\n\n")

      pairList[[i]] <-  cell_pairing(ATACpcs = chunkList[[i]][[ATAC.i]],
                                     RNApcs = chunkList[[i]][[RNA.i]],
                                     mode = "geodesic",...)
    }

    return(do.call('rbind',pairList))

  } else {
    return(cell_pairing(ATACpcs = ATACpcs,RNApcs = RNApcs,mode = "geodesic",...))
  }
}
