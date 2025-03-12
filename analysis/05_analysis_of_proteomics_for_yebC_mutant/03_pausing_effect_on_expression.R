library(tidyverse)

df_ms <- read_tsv("delta0316_ms.tsv")
df_rnaseq <- read_tsv('input_files/rnaseq_results.tsv')

pause_genes <- c(
  "SPy_0015",
  "SPy_0077",
  "SPy_0145",
  "SPy_0148",
  "SPy_0154",
  "SPy_0155",
  #"SPy_0212",
  "SPy_0259",
  "SPy_0269",
  "SPy_0288",
  "SPy_0296",
  "SPy_0306",
  "SPy_0312",
  "SPy_0385",
  "SPy_0453",
  "SPy_0458",
  "SPy_0460",
  "SPy_0510",
  "SPy_0542",
  "SPy_0543",
  #"SPy_0553",
  "SPy_0653",
  "SPy_0739",
  "SPy_0745",
  "SPy_0758",
  "SPy_0760",
  "SPy_0776",
  "SPy_0794",
  "SPy_0817",
  "SPy_0827",
  "SPy_0831",
  "SPy_0835",
  "SPy_0892",
  "SPy_0899",
  "SPy_0910",
  "SPy_0939",
  "SPy_1029",
  "SPy_1082",
  "SPy_1101",
  "SPy_1155",
  #"SPy_1177",
  "SPy_1246",
  "SPy_1286",
  "SPy_1368",
  "SPy_1370",
  "SPy_1434",
  "SPy_1509",
  "SPy_1547",
  "SPy_1549",
  "SPy_1562",
  "SPy_1568",
  "SPy_1582",
  "SPy_1629",
  "SPy_1637",
  "SPy_1639",
  "SPy_1678",
  "SPy_1704",
  "SPy_1765",
  "SPy_1770",
  "SPy_1781",
  "SPy_1827",
  "SPy_1849",
  "SPy_1916",
  "SPy_1917",
  "SPy_1919",
  "SPy_1973",
  "SPy_1992",
  "SPy_1998",
  "SPy_2019",
  "SPy_2034",
  "SPy_2073",
  "SPy_2096",
  "SPy_2111",
  "SPy_2174",
  "SPy_2191",
  "SPy_2211"
)

ppg_genes <- c(
  'SPy_0015',
  'SPy_0038',
  'SPy_0154',
  'SPy_0532',
  'SPy_0628',
  'SPy_0653',
  'SPy_0700',
  'SPy_0758',
  'SPy_0760',
  'SPy_0923',
  'SPy_1101',
  'SPy_1436',
  'SPy_1509',
  'SPy_1568',
  'SPy_1577',
  'SPy_1765',
  'SPy_1916',
  'SPy_1992',
  'SPy_2073'
)

# Merge RNA-seq and MS data
rnaseq <- df_rnaseq |> 
  select(
    locus_tag,
    old_locus_tag,
    logfc_rl_dw,
    logfc_rl_dc,
    logfc_rs_dw,
    logfc_rs_dc
  )

ms <- df_ms |> 
  select(
    old_locus_tag,
    logfc_lys_del_wt_log,
    logfc_lys_del_comp_log,
    logfc_lys_del_wt_stat,
    logfc_lys_del_comp_stat
  ) |> 
  rename(
    logfc_pl_dw = logfc_lys_del_wt_log,
    logfc_pl_dc = logfc_lys_del_comp_log,
    logfc_ps_dw = logfc_lys_del_wt_stat,
    logfc_ps_dc = logfc_lys_del_comp_stat
  )

df_ms_rnaseq <- inner_join(ms, rnaseq)

# Calculate protein to RNA logfc and find genes with pause sites and PPG/PPP motifs
df_ms_rnaseq <- df_ms_rnaseq |> 
  mutate(
    logfc_tl_dw = logfc_pl_dw - logfc_rl_dw,
    logfc_tl_dc = logfc_pl_dc - logfc_rl_dc,
    logfc_ts_dw = logfc_ps_dw - logfc_rs_dw,
    logfc_ts_dc = logfc_ps_dc - logfc_rs_dc
  ) |>
  mutate(
    pause = if_else(old_locus_tag %in% pause_genes, "pause", "other"),
    ppg = if_else(old_locus_tag %in% ppg_genes, "ppg", "other")
  )

##### Compare protein levels for genes with and without the pause sites or PPG/PPP motifs
df_prot_long <- df_ms_rnaseq |>
  select(
    old_locus_tag,
    logfc_pl_dw:logfc_ps_dc,
    pause,
    ppg
  )|> 
  pivot_longer(
    cols = starts_with("logfc"),
    names_to = "comparison",
    values_to = "logfc"
  ) |>
  mutate(
    comparison = factor(comparison, levels = c("logfc_pl_dw",
                                               "logfc_pl_dc",
                                               "logfc_ps_dw",
                                               "logfc_ps_dc"))
  )

# Wilcoxon test to check the difference between the genes with and without pause or PPG/PPP motif
pval_wilcox_pause <- df_prot_long |> 
  group_by(comparison) |> 
  do(w = wilcox.test(logfc ~ pause, data=., paired=FALSE)) %>% 
  summarise(comparison, wilcox = w$p.value)

pval_wilcox_ppg <- df_prot_long |> 
  group_by(comparison) |> 
  do(w = wilcox.test(logfc ~ ppg, data=., paired=FALSE)) %>% 
  summarise(comparison, wilcox = w$p.value)

# Plot logfc distributions for the proteins with and without pausing sites
ggplot(df_prot_long, aes(x = pause, y = logfc, color = pause)) +
  geom_boxplot(outlier.shape = NA) +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  scale_y_continuous(limits = c(-2, 2)) +
  facet_grid(~comparison) +
  scale_color_manual(values=c("black", "red")) +
  theme_classic()
ggsave('figures/boxplot_prot_DE.pdf', height = 5, width = 8)


##### Compare mRNA levels for genes with and without the pause sites
df_rna_long <- df_ms_rnaseq |>
  select(
    old_locus_tag,
    logfc_rl_dw:logfc_rs_dc,
    pause
  )|> 
  pivot_longer(
    cols = starts_with("logfc"),
    names_to = "comparison",
    values_to = "logfc"
  ) |>
  mutate(
    comparison = factor(comparison, levels = c("logfc_rl_dw",
                                               "logfc_rl_dc",
                                               "logfc_rs_dw",
                                               "logfc_rs_dc"))
  )

# Wilcoxon test to check the difference between the genes with and without pause
pval_wilcox_rna <- df_rna_long |> 
  group_by(comparison) |> 
  do(w = wilcox.test(logfc ~ pause, data=., paired=FALSE)) %>% 
  summarise(comparison, wilcox = w$p.value)

# Plot logfc distributions for the proteins with and without pausing sites
ggplot(df_rna_long, aes(x = pause, y = logfc, color = pause)) +
  geom_boxplot(outlier.shape = NA) +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  scale_y_continuous(limits = c(-2, 2)) +
  facet_grid(~comparison) +
  scale_color_manual(values=c("black", "red")) +
  theme_classic()
ggsave('figures/boxplot_rna_DE.pdf', height = 5, width = 8)

##### Compare translation efficiency for genes with and without the pause sites
df_trans_long <- df_ms_rnaseq |>
  select(
    old_locus_tag,
    logfc_tl_dw:logfc_ts_dc,
    pause
  )|> 
  pivot_longer(
    cols = starts_with("logfc"),
    names_to = "comparison",
    values_to = "logfc"
  ) |>
  mutate(
    comparison = factor(comparison, levels = c("logfc_tl_dw",
                                               "logfc_tl_dc",
                                               "logfc_ts_dw",
                                               "logfc_ts_dc"))
  )

# Wilcoxon test to check the difference between the genes with and without pause
pval_wilcox_trans <- df_trans_long |> 
  group_by(comparison) |> 
  do(w = wilcox.test(logfc ~ pause, data=., paired=FALSE)) %>% 
  summarise(comparison, wilcox = w$p.value)
  
# Plot logfc distributions for the proteins with and without pausing sites
ggplot(df_trans_long, aes(x = pause, y = logfc, color = pause)) +
  geom_boxplot(outlier.shape = NA) +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  scale_y_continuous(limits = c(-2, 2)) +
  facet_grid(~comparison) +
  scale_color_manual(values=c("black", "red")) +
  theme_classic()
ggsave('figures/boxplot_translation_DE.pdf', height = 5, width = 8)


########## End ##########

