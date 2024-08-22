library(data.table)
library(BuenColors)
source("call_variants.R")

if(FALSE){
  
  dn <- readRDS("../data/patient3/KY-2936_010424_DN_IGO_12437_FF_2_hg38_v20-mtMask.rds"); colnames(dn) <- paste0(colnames(dn), "-DN")
  dn_mgatk<- call_mutations_mgatk(dn[,dn$depth > 10])
  saveRDS(dn_mgatk, file = "../output/d3_dn_mgatk_called.rds")
  rm(dn); rm(dn_mgatk)
  
  dp <- readRDS("../data/patient3/KY-2936_010424_DP_IGO_12437_FF_1_hg38_v20-mtMask.rds"); colnames(dp) <- paste0(colnames(dp), "-DP")
  dp_mgatk<- call_mutations_mgatk(dp[,dp$depth > 10])
  saveRDS(dp_mgatk, file = "../output/d3_dp_mgatk_called.rds")
  rm(dp); rm(dp_mgatk)
  
  mo <- readRDS("../data/patient3/KY-2936_010424_Mono_IGO_12437_FF_3_hg38_v20-mtMask.rds"); colnames(mo) <- paste0(colnames(mo), "-MO")
  mo_mgatk<- call_mutations_mgatk(mo[,mo$depth > 10])
  saveRDS(mo_mgatk, file = "../output/d3_mo_mgatk_called.rds")
  rm(mo); rm(mo_mgatk)
}

df_dp <- data.frame(rowData(readRDS("../output/d3_dp_mgatk_called.rds")))
df_dn <- data.frame(rowData(readRDS("../output/d3_dn_mgatk_called.rds")))
df_mono <- data.frame(rowData(readRDS("../output/d3_mo_mgatk_called.rds")))
vars_dp <- df_dp %>%
  dplyr::filter(n_cells_conf_detected >= 3 & strand_correlation > 0.65 & 
           log10(vmr) > -2 & mean_coverage >= 10) %>% pull(variant)
vars_dn <- df_dn %>%
  dplyr::filter(n_cells_conf_detected >= 3 & strand_correlation > 0.65 & 
           log10(vmr) > -2 & mean_coverage >= 10) %>% pull(variant)
vars_mono <- df_mono %>%
  dplyr::filter(n_cells_conf_detected >= 3 & strand_correlation > 0.65 & 
           log10(vmr) > -2 & mean_coverage >= 10) %>% pull(variant)
varset <- unique(c(vars_dp, vars_dn, vars_mono))
varset <- varset[!(varset == "310T>C")]
df <- data.frame(variant = varset,
                 dp = df_dp[varset,"mean"], 
                 dn = df_dn[varset,"mean"], 
                 mono = df_mono[varset,"mean"]
) 
df %>% arrange(desc(mono)) %>% head()

library(ggpubr)
p1 <- ggplot(df, aes(x = dp, y = dn)) + 
  geom_point() + theme_bw() + labs(x = "DP allele frequency", y = "DN allele frequency") +
  scale_y_log10() + scale_x_log10() +stat_cor(method="spearman")

p2 <- ggplot(df, aes(x = mono, y = dp)) + 
  geom_point() + theme_bw() + labs(x = "Mono allele frequency", y = "DP allele frequency") +
  scale_y_log10() + scale_x_log10()+stat_cor(method="spearman")

p3 <- ggplot(df, aes(x = mono, y = dn)) + 
  geom_point() + theme_bw() + labs(x = "mono allele frequency", y = "DN allele frequency") +
  scale_y_log10() + scale_x_log10()+stat_cor(method="spearman")

cowplot::plot_grid(p1,p2,p3, nrow = 1)


library(pheatmap)
(cor(log10(df[,2:4])))


saveRDS(assays(readRDS("../output/d3_dp_mgatk_called.rds"))[["allele_frequency"]][varset,], file = "../output/d3_dp_afs.rds")
saveRDS(assays(readRDS("../output/d3_dn_mgatk_called.rds"))[["allele_frequency"]][varset,], file = "../output/d3_dn_afs.rds")
saveRDS(assays(readRDS("../output/d3_mo_mgatk_called.rds"))[["allele_frequency"]][varset,], file = "../output/d3_mono_afs.rds")



