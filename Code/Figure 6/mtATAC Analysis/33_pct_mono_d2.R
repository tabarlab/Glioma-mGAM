library(data.table)
library(BuenColors)
source("call_variants.R")

mo <- readRDS("../data/patient3/KY-2936_010424_Mono_IGO_12437_FF_3_hg38_v20-mtMask.rds"); colnames(mo) <- paste0(colnames(mo), "-MO")
mo_mgatk<- call_mutations_mgatk(mo[,mo$depth > 10])
dn <- readRDS("../data/patient3/KY-2936_010424_DN_IGO_12437_FF_2_hg38_v20-mtMask.rds"); colnames(dn) <- paste0(colnames(dn), "-DN")
dn_mgatk<- call_mutations_mgatk(dn[,dn$depth > 10])
dp <- readRDS("../data/patient3/KY-2936_010424_DP_IGO_12437_FF_1_hg38_v20-mtMask.rds"); colnames(dp) <- paste0(colnames(dp), "-DP")
dp_mgatk<- call_mutations_mgatk(dp[,dp$depth > 10])

df_dp <- data.frame(rowData(readRDS("../output/d3_dp_mgatk_called.rds")))
df_dn <- data.frame(rowData(readRDS("../output/d3_dn_mgatk_called.rds")))
df_mono <- data.frame(rowData(readRDS("../output/d3_mo_mgatk_called.rds")))

dim(df_mono)

vars_dp <- df_dp %>%
  dplyr::filter(n_cells_conf_detected >= 1 & strand_correlation > 0.65 & 
           log10(vmr) > -2 & mean_coverage >= 10 & mean < 0.5) %>% pull(variant)
vars_dn <- df_dn %>%
  dplyr::filter(n_cells_conf_detected >= 1 & strand_correlation > 0.65 & 
           log10(vmr) > -2 & mean_coverage >= 10 & mean < 0.5) %>% pull(variant)
vars_mono <- df_mono %>%
  dplyr::filter(n_cells_conf_detected >= 1 & strand_correlation > 0.65 & 
           log10(vmr) > -2 & mean_coverage >= 10 & mean < 0.5) %>% pull(variant)


mean(colSums(assays(mo_mgatk)[["allele_frequency"]][vars_mono,] >= 0.1) >= 1)
mean(colSums(assays(dn_mgatk)[["allele_frequency"]][vars_mono,] >= 0.1) >= 1)
mean(colSums(assays(dp_mgatk)[["allele_frequency"]][vars_mono,] >= 0.1) >= 1)


