library(tidyverse)

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

# Read and modify the data

ifelse(!dir.exists("figures"), dir.create("figures"),
       "Folder 'figures' exists already")

df_expr <- read_delim("input_files/24C003_DVI_df_ProteinLevelData.tsv")
df_comp <- read_delim("input_files/24C003_DVI_df_comparisonResult.tsv")
df_anno <- read_delim("input_files/Proteome annotation.tsv")

df_expr <- df_expr |> 
  select(
    uniprot = Protein,
    sample = GROUP,
    replicate = SUBJECT,
    log2intens = LogIntensities
  ) |> mutate(
    sample = str_to_lower(sample),
    sample = str_replace_all(sample, "lysate", "lys"),
    sample = str_replace_all(sample, "ribo", "rib")
  )

df_comp <- df_comp |> 
  select(
    uniprot = Protein,
    comparison = Label,
    logfc = log2FC,
    pval = pvalue,
    padj = adj.pvalue
  ) |> 
  mutate(
    logfc = as.numeric(logfc),
    comparison = str_to_lower(comparison),
    comparison = str_replace_all(comparison, "lysate", "lys"),
    comparison = str_replace_all(comparison, "ribo", "rib"),
    comparison = str_replace_all(comparison, " vs ", "_vs_")
  )
  
df_anno <- df_anno |> 
  select(
    uniprot = UniProt,
    old_locus_tag = ENSG,
    symbol = 'Gene name',
    anno_interpro = 'Interpro name'
  ) |> 
  mutate(
    uniprot = str_sub(uniprot, 1, 6),
    old_locus_tag = str_replace(old_locus_tag, ";.+", "")
  )

df_expr <- df_expr |> 
  left_join(df_anno) |> 
  relocate(c('old_locus_tag', 'symbol', 'anno_interpro'), .after = uniprot)
df_comp <- df_comp |> 
  left_join(df_anno) |> 
  relocate(c('old_locus_tag', 'symbol', 'anno_interpro'), .after = uniprot)

df_expr <- left_join(df_expr, df_anno)
df_comp <- left_join(df_comp, df_anno)

# Filter the data

df_expr <- df_expr |> 
  filter(!sample %in% c("rib_wt_log", "rib_del_log", "rib_comp_log",
                        "rib_wt_stat", "rib_del_stat", "rib_comp_stat"))

df_comp <- df_comp |> 
  filter(
    comparison %in% c(
      "lys_del_log_vs_lys_wt_log",
      "lys_comp_log_vs_lys_del_log",
      "lys_del_stat_vs_lys_wt_stat",
      "lys_comp_stat_vs_lys_del_stat"
    )
  )

# Merge the data

df_expr_wide_rep <- df_expr |> 
  filter(
    str_starts(sample, 'lys')
  ) |> 
  pivot_wider(
    names_from = c(sample, replicate),
    values_from = log2intens
  )

df_expr_wide <- df_expr |>
  filter(
    str_starts(sample, 'lys')
  ) |> 
  mutate(intens = 2^log2intens) |> 
  group_by(uniprot, old_locus_tag, symbol, anno_interpro, sample) |> 
  summarise(log2intens_avg = log2(mean(intens))) |> 
  ungroup() |> 
  pivot_wider(
    names_from = sample,
    values_from = log2intens_avg
  )

df_comp_wide <- df_comp |>
  pivot_wider(
    id_cols = uniprot:anno_interpro,
    names_from = comparison,
    values_from = c("logfc", "pval", "padj")
  ) |> 
  mutate(
    logfc_lys_del_wt_log = logfc_lys_del_log_vs_lys_wt_log,
    pval_lys_del_wt_log = pval_lys_del_log_vs_lys_wt_log,
    padj_lys_del_wt_log = padj_lys_del_log_vs_lys_wt_log,
    
    logfc_lys_del_comp_log = -logfc_lys_comp_log_vs_lys_del_log,
    pval_lys_del_comp_log = pval_lys_comp_log_vs_lys_del_log,
    padj_lys_del_comp_log = padj_lys_comp_log_vs_lys_del_log,
    
    logfc_lys_del_wt_stat = logfc_lys_del_stat_vs_lys_wt_stat,
    pval_lys_del_wt_stat = pval_lys_del_stat_vs_lys_wt_stat,
    padj_lys_del_wt_stat = padj_lys_del_stat_vs_lys_wt_stat,
    
    logfc_lys_del_comp_stat = -logfc_lys_comp_stat_vs_lys_del_stat,
    pval_lys_del_comp_stat = pval_lys_comp_stat_vs_lys_del_stat,
    padj_lys_del_comp_stat = padj_lys_comp_stat_vs_lys_del_stat,
    
  ) |> 
  select(!old_locus_tag:padj_lys_del_stat_vs_lys_wt_stat)

df_ms <- left_join(df_expr_wide, df_comp_wide, by = "uniprot")
df_ms_rep <- left_join(df_expr_wide_rep, df_comp_wide, by = "uniprot")

# Identify DE genes in Log and Stat
df_ms <- df_ms |> 
  filter(old_locus_tag != 'SPy_0316') |> 
  mutate(
    de_del_wt_log = if_else(padj_lys_del_wt_log < .05 & !between(logfc_lys_del_wt_log, -1, 1), TRUE, FALSE),
    de_del_comp_log = if_else(padj_lys_del_comp_log < .05 & !between(logfc_lys_del_comp_log, -1, 1), TRUE, FALSE),
    de_del_wt_stat = if_else(padj_lys_del_wt_stat < .05 & !between(logfc_lys_del_wt_stat, -1, 1), TRUE, FALSE),
    de_del_comp_stat = if_else(padj_lys_del_comp_stat < .05 & !between(logfc_lys_del_comp_stat, -1, 1), TRUE, FALSE),
    pause = if_else(old_locus_tag %in% pause_genes, TRUE, FALSE)
  )

df_m_reps <- df_ms_rep |> 
  filter(old_locus_tag != 'SPy_0316') |> 
  mutate(
    de_del_wt_log = if_else(padj_lys_del_wt_log < .05 & !between(logfc_lys_del_wt_log, -1, 1), TRUE, FALSE),
    de_del_comp_log = if_else(padj_lys_del_comp_log < .05 & !between(logfc_lys_del_comp_log, -1, 1), TRUE, FALSE),
    de_del_wt_stat = if_else(padj_lys_del_wt_stat < .05 & !between(logfc_lys_del_wt_stat, -1, 1), TRUE, FALSE),
    de_del_comp_stat = if_else(padj_lys_del_comp_stat < .05 & !between(logfc_lys_del_comp_stat, -1, 1), TRUE, FALSE),
    pause = if_else(old_locus_tag %in% pause_genes, TRUE, FALSE)
  )

# Save the combined datasets
write_tsv(df_ms, file = "delta0316_ms.tsv")
write_tsv(df_ms_rep, file = "delta0316_ms_rep.tsv")

########## End ##########

