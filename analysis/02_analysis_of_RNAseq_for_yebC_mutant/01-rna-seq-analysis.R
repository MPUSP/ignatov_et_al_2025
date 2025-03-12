library(tidyverse)
library(ggpubr)
library(ggvenn)
library(ggrepel)

ifelse(!dir.exists("Figures"), dir.create("Figures"), "Folder 'Figures' exists already")

# Read KEGG meta-data
kegg <- read.delim(
  "input_files/gene_go_kegg_meta_info.tsv"
) |> select(
  "locus_tag",
  "symbol",
  "old_locus_tag",
  "product",
  "GO_process"
)

# Read the table with RNA-seq TPM values
tpm <- read.delim(
  "input_files/all_tpm_counts.tsv"
) |> select(
  locus_tag,
  biotype,
  tpm_lw1 = yebC.log_w_1,
  tpm_lw2 = yebC.log_w_2,
  tpm_lw3 = yebC.log_w_3,
  tpm_ld1 = yebC.log_d_1,
  tpm_ld2 = yebC.log_d_2,
  tpm_ld3 = yebC.log_d_3,
  tpm_lc1 = yebC.log_c_1,
  tpm_lc2 = yebC.log_c_2,
  tpm_lc3 = yebC.log_c_3,
  tpm_sw1 = yebC.stat_w_1,
  tpm_sw2 = yebC.stat_w_2,
  tpm_sw3 = yebC.stat_w_3,
  tpm_sd1 = yebC.stat_d_1,
  tpm_sd2 = yebC.stat_d_2,
  tpm_sd3 = yebC.stat_d_3,
  tpm_sc1 = yebC.stat_c_1,
  tpm_sc2 = yebC.stat_c_2,
  tpm_sc3 = yebC.stat_c_3
)

# Read the results of DE analysis
de_rl_dw <- read.delim(
  "input_files/log_d-vs-log_w.tsv"
) |> select(
  locus_tag = gene,
  base_mean_rl_dw = baseMean,
  logfc_rl_dw = log2FoldChange,
  padj_rl_dw = padj
)

de_rl_dc <- read.delim(
  "input_files/log_d-vs-log_c.tsv"
) |> select(
  locus_tag = gene,
  base_mean_rl_dc = baseMean,
  logfc_rl_dc = log2FoldChange,
  padj_rl_dc = padj
)

de_rl_dwc <- read.delim(
  "input_files/log_d-vs-log_wc.tsv"
) |> select(
  locus_tag = gene,
  base_mean_rl_dwc = baseMean,
  logfc_rl_dwc = log2FoldChange,
  padj_rl_dwc = padj
)

de_rs_dw <- read.delim(
  "input_files/stat_d-vs-stat_w.tsv"
) |> select(
  locus_tag = gene,
  base_mean_rs_dw = baseMean,
  logfc_rs_dw = log2FoldChange,
  padj_rs_dw = padj
)

de_rs_dc <- read.delim(
  "input_files/stat_d-vs-stat_c.tsv"
) |> select(
  locus_tag = gene,
  base_mean_rs_dc = baseMean,
  logfc_rs_dc = log2FoldChange,
  padj_rs_dc = padj
)

de_rs_dwc <- read.delim(
  "input_files/stat_d-vs-stat_wc.tsv"
) |> select(
  locus_tag = gene,
  base_mean_rs_dwc = baseMean,
  logfc_rs_dwc = log2FoldChange,
  padj_rs_dwc = padj
)

# Merge the tables and identify DE genes
rnaseq <- left_join(kegg, tpm)
rnaseq <- left_join(rnaseq, de_rl_dw)
rnaseq <- left_join(rnaseq, de_rl_dc)
rnaseq <- left_join(rnaseq, de_rl_dwc)
rnaseq <- left_join(rnaseq, de_rs_dw)
rnaseq <- left_join(rnaseq, de_rs_dc)
rnaseq <- left_join(rnaseq, de_rs_dwc)

rnaseq <- rnaseq |> mutate(
  de_rl_dw = if_else(padj_rl_dw < .05 & !between(logfc_rl_dw, -1, 1), TRUE, FALSE),
  de_rl_dc = if_else(padj_rl_dc < .05 & !between(logfc_rl_dc, -1, 1), TRUE, FALSE),
  de_rl_dwc = if_else(padj_rl_dwc < .05 & !between(logfc_rl_dwc, -1, 1), TRUE, FALSE),
  de_rs_dw = if_else(padj_rs_dw < .05 & !between(logfc_rs_dw, -1, 1), TRUE, FALSE),
  de_rs_dc = if_else(padj_rs_dc < .05 & !between(logfc_rs_dc, -1, 1), TRUE, FALSE),
  de_rs_dwc = if_else(padj_rs_dwc < .05 & !between(logfc_rs_dwc, -1, 1), TRUE, FALSE)
) |> mutate(
  de_cat_rl = case_when(
    de_rl_dw & de_rl_dc & sign(logfc_rl_dw) == sign(logfc_rl_dc) ~ "rlwc",
    de_rl_dw ~ "rlw",
    de_rl_dc ~ "rlc"
  ),
  de_cat_rs = case_when(
    de_rs_dw & de_rs_dc & sign(logfc_rs_dw) == sign(logfc_rs_dc) ~ "rswc",
    de_rs_dw ~ "rsw",
    de_rs_dc ~ "rsc"
  ),
  de_cat_rl = factor(de_cat_rl, levels = c('rlw', 'rlwc', 'rlc')),
  de_cat_rs = factor(de_cat_rs, levels = c('rsw', 'rswc', 'rsc'))
)

# Save the merged results. Only keep the protein-coding genes
rnaseq <- filter(rnaseq, biotype == "mRNA")
rnaseq <- rnaseq[complete.cases(rnaseq[,43:48]),]
write_tsv(rnaseq, "rnaseq_results.tsv")

# Remove Spy_0316
rnaseq <- filter(rnaseq, old_locus_tag != "SPy_0316")

# Select genes with the highest changes in Log and Stat
res_rl <- filter(rnaseq, !is.na(de_cat_rl))
res_rlwc <- filter(res_rl, de_cat_rl == "rlwc")
res_rlwc_max <- res_rlwc %>% 
  arrange(desc(abs(logfc_rl_dw) + abs(logfc_rl_dc))) %>% 
  slice(1:10)

res_rs <- filter(rnaseq, !is.na(de_cat_rs))
res_rswc <- filter(res_rs, de_cat_rs == "rswc")
res_rswc_max <- res_rswc %>% 
  arrange(desc(abs(logfc_rs_dw) + abs(logfc_rs_dc))) %>% 
  slice(1:10)

# Select SpeB
res_rswc_speB <- filter(res_rswc, old_locus_tag == 'SPy_2039')

# Venn diagrams for the pairwise comparisons
p_venn_log <- ggvenn(list(d0316_vs_WT = rnaseq$locus_tag[rnaseq$de_rl_dw],
                          d0316_vs_comp = rnaseq$locus_tag[rnaseq$de_rl_dc]),
             fill_color = c("#A6CEE3", "#B2DF8A"), show_percentage = FALSE)

p_venn_stat <- ggvenn(list(d0316_vs_WT = rnaseq$locus_tag[rnaseq$de_rs_dw],
                           d0316_vs_comp = rnaseq$locus_tag[rnaseq$de_rs_dc]),
             fill_color = c("#A6CEE3", "#B2DF8A"), show_percentage = FALSE)

# Correlation plots for the pairwise comparisons
p_cor_log <- ggplot(res_rl, mapping = aes(x = logfc_rl_dw, y = logfc_rl_dc, color = de_cat_rl)) +
  geom_point() +
  # geom_label_repel(res_rlwc_max, mapping = aes(x = logfc_rl_dw, y = logfc_rl_dc, label = symbol),
  #                 max.overlaps = Inf, show.legend = FALSE) +
  coord_fixed(xlim = c(-3, 3), ylim = c(-3, 3)) +
  scale_color_brewer(name = "Diff expr genes", labels = c("d0316 vs. WT", "both","d0316 vs. comp"), palette="Paired") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position="none") +
  labs(title = "DE genes in d0316 relative to WT and comp in Log growth phase",
       x = "d0316 vs. WT", y = "d0316 vs. comp")

p_cor_stat <- ggplot(res_rs, mapping = aes(x = logfc_rs_dw, y = logfc_rs_dc, color = de_cat_rs)) +
  geom_point() +
  # geom_label_repel(res_rswc_max, mapping = aes(x = logfc_rs_dw, y = logfc_rs_dc, label = symbol),
  #                 max.overlaps = Inf, show.legend = FALSE) +
  coord_fixed(xlim = c(-3, 3), ylim = c(-3, 3)) +
  scale_color_brewer(name = "Diff expr genes", labels = c("d0316 vs. WT", "both","d0316 vs. comp"), palette="Paired") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position="none") +
  labs(title = "DE genes in d0316 relative to WT and comp in Stat growth phase",
       x = "d0316 vs. WT", y = "d0316 vs. comp")

# MA plots for the pairwise comparisons
p_ma_log <- ggplot(rnaseq) +
  geom_point(aes(x = log2(base_mean_rl_dw), y = logfc_rl_dw),
             alpha = 0.1, size = 2) +
  geom_point(data = res_rlwc,
             aes(x = log2(base_mean_rl_dw), y = logfc_rl_dw),
             color = "dodgerblue3", size = 2) + 
  geom_label_repel(data = res_rlwc_max,
             aes(x = log2(base_mean_rl_dw), y = logfc_rl_dw, label = symbol),
             max.overlaps = Inf) +
  ggtitle(label = "DE genes in d0316 relative to WT in Log growth phase") +
  scale_x_continuous(name = "Mean expression, log2",
                     limits = c(0, 20),
                     breaks = seq(0, 20, by = 2)) +
  scale_y_continuous(name = "Fold change, log2",
                     limits = c(-3, 3)) +
  theme_bw() +
  theme(aspect.ratio=0.6,
        panel.grid.minor = element_blank())
  

p_ma_stat <- ggplot(rnaseq) +
  geom_point(aes(x = log2(base_mean_rs_dw), y = logfc_rs_dw),
             alpha = 0.1, size = 2) +
  geom_point(data = res_rswc,
             aes(x = log2(base_mean_rs_dw), y = logfc_rs_dw),
             color = "dodgerblue3", size = 2) + 
  geom_label_repel(data = res_rswc_max,
                  aes(x = log2(base_mean_rs_dw), y = logfc_rs_dw, label = symbol),
                  max.overlaps = Inf) +
  ggtitle(label = "DE genes in d0316 relative to WT in Stat growth phase") +
  scale_x_continuous(name = "Mean expression, log2",
                     limits = c(0, 20),
                     breaks = seq(0, 20, by = 2)) +
  scale_y_continuous(name = "Fold change, log2",
                     limits = c(-3, 3)) +
  theme_bw() +
  theme(aspect.ratio=0.6,
        panel.grid.minor = element_blank())

pdf(file = "./Figures/rnaseq_main.pdf",
    width = 12, height = 15)
ggarrange(plotlist = list(
  p_venn_log,
  p_venn_stat,
  p_cor_log,
  p_cor_stat,
  p_ma_log,
  p_ma_stat
  ), nrow = 3, ncol = 2)
dev.off()

pdf(file = "./Figures/rnaseq_main_2.pdf",
    width = 18, height = 12)
ggarrange(plotlist = list(
  p_venn_log, p_cor_log, p_ma_log,
  p_venn_stat, p_cor_stat, p_ma_stat
), nrow = 2, ncol = 3)
dev.off()

# MA plot with speB highlighted
p_ma_stat_speB <- ggplot(rnaseq) +
  geom_point(aes(x = log2(base_mean_rs_dw), y = logfc_rs_dw),
             alpha = 0.1, size = 2) +
  geom_point(data = res_rswc,
             aes(x = log2(base_mean_rs_dw), y = logfc_rs_dw),
             color = "dodgerblue3", size = 2) + 
  geom_text_repel(data = res_rswc_speB,
                   aes(x = log2(base_mean_rs_dw), y = logfc_rs_dw, label = symbol),
                   max.overlaps = Inf) +
  ggtitle(label = "DE genes in d0316 relative to WT in Stat growth phase") +
  scale_x_continuous(name = "Mean expression, log2",
                     limits = c(0, 20),
                     breaks = seq(0, 20, by = 2)) +
  scale_y_continuous(name = "Fold change, log2",
                     limits = c(-3, 3)) +
  theme_bw() +
  theme(aspect.ratio=0.6,
        panel.grid.minor = element_blank())

# Save all plots separately

ggsave(filename = "./Figures/cor_log.pdf",
       plot = p_cor_log,
       width = 3, height = 3)
ggsave(filename = "./Figures/cor_stat.pdf",
       plot = p_cor_stat,
       width = 3, height = 3)


ggsave(filename = "./Figures/ma_log.pdf",
       plot = p_ma_log,
       width = 5, height = 5)
ggsave(filename = "./Figures/ma_stat.pdf",
       plot = p_ma_stat,
       width = 5, height = 5)

ggsave(filename = "./Figures/ma_stat_speB.pdf",
       plot = p_ma_stat_speB,
       width = 5, height = 5)
