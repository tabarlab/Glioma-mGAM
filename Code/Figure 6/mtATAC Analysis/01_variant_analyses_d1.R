library(data.table)
library(BuenColors)
source("call_variants.R")

if(FALSE){
  
  dn <- readRDS("../data/KY-DN_hg38_v20-mtMask.rds"); colnames(dn) <- paste0(colnames(dn), "-DN")
  dn_mgatk<- call_mutations_mgatk(dn[,dn$depth > 10])
  saveRDS(dn_mgatk, file = "../output/dn_mgatk_called.rds")
  rm(dn); rm(dn_mgatk)
  
  dp <- readRDS("../data/KY-DP_hg38_v20-mtMask.rds"); colnames(dp) <- paste0(colnames(dp), "-DP")
  dp_mgatk<- call_mutations_mgatk(dp[,dp$depth > 10])
  saveRDS(dp_mgatk, file = "../output/dp_mgatk_called.rds")
  rm(dp); rm(dp_mgatk)
  
  mo <- readRDS("../data/KY-Mono_hg38_v20-mtMask.rds"); colnames(mo) <- paste0(colnames(mo), "-MO")
  mo_mgatk<- call_mutations_mgatk(mo[,mo$depth > 15])
  saveRDS(mo_mgatk, file = "../output/mo_mgatk_called.rds")
  rm(mo); rm(mo_mgatk)
}
df_dp <- data.frame(rowData(readRDS("../output/dp_mgatk_called.rds")))
df_dn <- data.frame(rowData(readRDS("../output/dn_mgatk_called.rds")))
df_mono <- data.frame(rowData(readRDS("../output/mo_mgatk_called.rds")))
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
  scale_y_log10() + scale_x_log10() 

p2 <- ggplot(df, aes(x = dp*100, y = mono*100, color = variant %in% c("8780T>C"))) + 
  geom_point() + pretty_plot(fontsize = 8) + L_border() + labs(x = "DP allele frequency", y = "mono allele frequency") +
  scale_y_log10() + scale_x_log10()

p2
p3 <- ggplot(df, aes(x = mono, y = dn)) + 
  geom_point() + theme_bw() + labs(x = "monocyte allele frequency", y = "DN allele frequency") +
  scale_y_log10() + scale_x_log10()

cowplot::ggsave2(p2, file = "../figures/correlation_DPmono.pdf",
                   width = 2, height = 2)


library(pheatmap)
(cor(log10(df[,2:4])))

df %>% filter(dn > 0.5)

saveRDS(assays(readRDS("../output/dp_mgatk_called.rds"))[["allele_frequency"]][varset,], file = "../output/dp_afs.rds")
saveRDS(assays(readRDS("../output/dn_mgatk_called.rds"))[["allele_frequency"]][varset,], file = "../output/dn_afs.rds")
saveRDS(assays(readRDS("../output/mo_mgatk_called.rds"))[["allele_frequency"]][varset,], file = "../output/mono_afs.rds")



