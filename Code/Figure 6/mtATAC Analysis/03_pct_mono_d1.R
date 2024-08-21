library(data.table)
library(BuenColors)
source("call_variants.R")

dn <- readRDS("../data/patient1/KY-DN_hg38_v20-mtMask.rds"); colnames(dn) <- paste0(colnames(dn), "-DN")
dn_mgatk<- call_mutations_mgatk(dn[,dn$depth > 10])

dp <- readRDS("../data/patient1/KY-DP_hg38_v20-mtMask.rds"); colnames(dp) <- paste0(colnames(dp), "-DP")
dp_mgatk<- call_mutations_mgatk(dp[,dp$depth > 10])

mo <- readRDS("../data/patient1/KY-Mono_hg38_v20-mtMask.rds"); colnames(mo) <- paste0(colnames(mo), "-MO")
mo_mgatk<- call_mutations_mgatk(mo[,mo$depth > 15])

df_dp <- data.frame(rowData(readRDS("../output/dp_mgatk_called.rds")))
df_dn <- data.frame(rowData(readRDS("../output/dn_mgatk_called.rds")))
df_mono <- data.frame(rowData(readRDS("../output/mo_mgatk_called.rds")))

vars_dp <- df_dp %>%
  dplyr::filter(n_cells_conf_detected >= 1 & strand_correlation > 0.65 & 
           log10(vmr) > -2 & mean_coverage >= 10 & mean < 0.5) %>% pull(variant)
vars_dn <- df_dn %>%
  dplyr::filter(n_cells_conf_detected >= 1 & strand_correlation > 0.65 & 
           log10(vmr) > -2 & mean_coverage >= 10& mean < 0.5) %>% pull(variant)
vars_mono <- df_mono %>%
  dplyr::filter(n_cells_conf_detected >= 1 & strand_correlation > 0.65 & 
           log10(vmr) > -2 & mean_coverage >= 10& mean < 0.5) %>% pull(variant)


vm <- mean(colSums(assays(mo_mgatk)[["allele_frequency"]][vars_mono,] >= 0.1) >= 1)
vdn <- mean(colSums(assays(dn_mgatk)[["allele_frequency"]][vars_mono,] >= 0.1) >= 1)
vdp <- mean(colSums(assays(dp_mgatk)[["allele_frequency"]][vars_mono,] >= 0.1) >= 1)

data.frame(
  value = c(vdn, vdp, vm)*100,
  what = c("DN", "DP", "Mono")
) %>%
  ggplot(aes(x = what, y = value, fill = what)) + 
  geom_bar(stat = "identity", color = "black", width = 0.8) + 
  scale_y_continuous(expand = c(0,0)) + 
  coord_flip() +
  labs(x = "", y = "") + 
  scale_fill_manual(values = c("#6FC276","#B19CD9","#A0D5F6")) +
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none") -> p1
cowplot::ggsave2(p1, file = "../figures/bars_prop.pdf", width = 1.2, height = 2)
