# A pathogenic subpopulation of human glioma associated macrophages linked to glioma progression

This GitHub repository includes all scripts used to perform the analyses in the paper entitled "A pathogenic subpopulation of human glioma associated macrophages linked to glioma progression". For further details please consult the methods section of the manuscript.

## Data Processing and Analysis
The folder "Data Processing and Analysis" contains all the initial quality control steps performed prior to the main analysis, followed by the downstream scRNA, scATAC, MiloR, optMatch, and FigR. This Folder contains all the scripts for importing the cellranger outputs into Seurat, QC for scRNA-seq and scATAC-seq and downstream analyses, MiloR analysis identifying "archetypal" cells, matching scRNA and scATAC cells, and running the optMatch and FigR analysis.

## Figure 1
The folder "Figure 1" contains all the scripts needed to reproduce the different panels from Figure 1: UMAP of single-cell RNA and ATAC cells (Figure 1B), Dot plot to show cell type marker expression levels (Figure 1C), Differential gene expression analysis across wildtype, mutant and normal brain samples (Related to Figure 1E), Gene expression clustering and heatmap showing genes between wildtype, mutant and normal brain myeloid populations (Figure 1F), and gene expression correlation between wildtype, mutant and normal brain samples (Figure 1G).

## Figure 2
The folder "Figure 2" contains all the scripts needed to perform the MiloR and FigR visualizations using the FigR analysis results from "Data Preprocessing/05 - optMatch Wrapper and FigR Analysis". The scripts are related to the following panels from Figure 2: Consensus clustering of proprietary and public datasets (Figure 2A), Alluvial plot showing cell type calling between different public reference sets (Figure 2B), Alluvial plot showing Cell type calling for cells originating from wildtype, mutant and normal brain samples using GBMap reference (Figure 2C), Consensus clustering of proprietary GAMs/myeloid cells with public datasets (Figure 2E), Differential abundance testing using MiloR on consensus dataset (Figure 2F), Optmatch and FigR on differentially abundant populations (Figure 2G), FigR plot showing DORC enrichment between different conditions (Figure 2H).

## Figure 3
The folder "Figure 3" contains all the scripts needed to perform the gene regulatory network visualization using the FigR analysis results from "Data Preprocessing/05 - optMatch Wrapper and FigR Analysis", in addition to the gene pair sensitivity analysis for surface marker pair identification. The scripts are related to the following panels from Figure 3: Gene regulatory network showing DORC relationships in all 3 conditions (Figure 3A), Heatmap showing gene pair sensitivity analysis to identify optimal cell surface marker pairing (Figure 3E), mGAM percent bar plots from scRNA samples, and flow cytometry samples (Figures 3C, 3H).

## Figure 5
The folder "Figure 5" contains all the scripts needed to perform the mtATAC analysis and the corresponding visualizations. The scripts are related to the following panels from Figure 5: mtATAC samples UMAP, FOSL2 regulon accessibility profile violin plot, monocyte and mGAM heteroplasmy correlation plot, and percentage of cells with high-confidence monocyte-derived variants within each population (Figure 5B-F).

## Figure 7
The folder "Figure 7" contains all the scripts needed to perform the bulk RNA-seq analysis for the BT142 co-culture experiment data and the corresponding visualizations. The scripts are related to the following panels from Figure 7: Volcano plot for BT142 versus BT142 co-cultured with mGAMs (Figure 7E), BT142 meta-module scoring (Figure 7F), incoming and outgoing connections in mGAM from CellChat (Figure 7K).