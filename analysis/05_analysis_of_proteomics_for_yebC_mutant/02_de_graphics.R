library(tidyverse)
library(ggpubr)
library(ggvenn)
library(ggrepel)

df_ms <- read_tsv("delta0316_ms.tsv")

# MA plots for the pairwise comparisons
ggplot(df_ms) +
  geom_point(aes(x = lys_wt_log, y = logfc_lys_del_wt_log),
             alpha = 0.1, size = 2) +
  geom_point(data = filter(df_ms, de_del_wt_log),
             aes(x = lys_wt_log, y = logfc_lys_del_wt_log),
             color = "dodgerblue3", size = 2) +
  geom_label_repel(data = filter(df_ms, de_del_wt_log),
             aes(x = lys_wt_log, y = logfc_lys_del_wt_log, label = symbol),
             color = "dodgerblue3") +
  ggtitle(label = "d0316 vs. WT in Log") +
  scale_x_continuous(name = "Mean expression, log2",
                     limits = c(0, 20),
                     breaks = seq(0, 20, by = 2)) +
  scale_y_continuous(name = "Fold change, log2",
                     limits = c(-4, 4)) +
  theme_bw() +
  theme(aspect.ratio=0.6,
        panel.grid.minor = element_blank())

ggsave("figures/ma_plot_log.pdf")


ggplot(df_ms) +
  geom_point(aes(x = lys_wt_stat, y = logfc_lys_del_wt_stat),
             alpha = 0.1, size = 2) +
  geom_point(data = filter(df_ms, de_del_wt_stat),
             aes(x = lys_wt_stat, y = logfc_lys_del_wt_stat),
             color = "dodgerblue3", size = 2) +
  geom_label_repel(data = filter(df_ms, de_del_wt_stat),
             aes(x = lys_wt_stat, y = logfc_lys_del_wt_stat, label = symbol),
             color = "dodgerblue3") +
  ggtitle(label = "d0316 vs. WT in Stat") +
  scale_x_continuous(name = "Mean expression, log2",
                     limits = c(0, 20),
                     breaks = seq(0, 20, by = 2)) +
  scale_y_continuous(name = "Fold change, log2",
                     limits = c(-4, 4)) +
  theme_bw() +
  theme(aspect.ratio=0.6,
        panel.grid.minor = element_blank())

ggsave("figures/ma_plot_stat.pdf")


ggplot(df_ms) +
  geom_point(aes(x = lys_comp_log, y = logfc_lys_del_comp_log),
             alpha = 0.1, size = 2) +
  ggtitle(label = "d0316 vs. comp in Log") +
  scale_x_continuous(name = "Mean expression, log2",
                     limits = c(0, 20),
                     breaks = seq(0, 20, by = 2)) +
  scale_y_continuous(name = "Fold change, log2",
                     limits = c(-4, 4)) +
  theme_bw() +
  theme(aspect.ratio=0.6,
        panel.grid.minor = element_blank())

ggplot(df_ms) +
  geom_point(aes(x = lys_comp_stat, y = logfc_lys_del_comp_stat),
             alpha = 0.1, size = 2) +
  ggtitle(label = "d0316 vs. comp in Stat") +
  scale_x_continuous(name = "Mean expression, log2",
                     limits = c(0, 20),
                     breaks = seq(0, 20, by = 2)) +
  scale_y_continuous(name = "Fold change, log2",
                     limits = c(-5, 5)) +
  theme_bw() +
  theme(aspect.ratio=0.6,
        panel.grid.minor = element_blank())


# MA plots for the pairwise comparisons indicating genes with YebC pauses
ggplot(df_ms) +
  geom_point(aes(x = lys_wt_log, y = logfc_lys_del_wt_log),
             alpha = 0.1, size = 2) +
  geom_point(data = filter(df_ms, pause),
             aes(x = lys_wt_log, y = logfc_lys_del_wt_log),
             color = "dodgerblue3", size = 2)

ggplot(df_ms) +
  geom_point(aes(x = lys_comp_log, y = logfc_lys_del_comp_log),
             alpha = 0.1, size = 2) +
  geom_point(data = filter(df_ms, pause),
             aes(x = lys_comp_log, y = logfc_lys_del_comp_log),
             color = "dodgerblue3", size = 2)

ggplot(df_ms) +
  geom_point(aes(x = lys_wt_stat, y = logfc_lys_del_wt_stat),
             alpha = 0.1, size = 2) +
  geom_point(data = filter(df_ms, pause),
             aes(x = lys_wt_stat, y = logfc_lys_del_wt_stat),
             color = "dodgerblue3", size = 2)

ggplot(df_ms) +
  geom_point(aes(x = lys_comp_stat, y = logfc_lys_del_comp_stat),
             alpha = 0.1, size = 2) +
  geom_point(data = filter(df_ms, pause),
             aes(x = lys_comp_stat, y = logfc_lys_del_comp_stat),
             color = "dodgerblue3", size = 2)


