library(tidyverse)

ifelse(!dir.exists("figures"), dir.create("figures"), "Folder 'figures' exists already")

rsd <- read_csv("riboseq.csv")

# Calculate the total coverage around start and stop codons
cvg_start <- rsd |>
  group_by(nt_locus_coord) |> 
  summarise(cvg_by_pos = sum(cvg_wt_1, na.rm = TRUE))
cvg_stop <- rsd |>
  group_by(nt_locus_coord_rev) |> 
  summarise(cvg_by_pos = sum(cvg_wt_1, na.rm = TRUE))

# Plot the cumulative coverage for the first an last 75 nt
ggplot(cvg_start, aes(nt_locus_coord, cvg_by_pos)) +
  geom_line(col = 'dodgerblue3') +
  scale_x_continuous(limits = c(-15, 60),
                     breaks = seq(-15, 60, 15),
                     minor_breaks = NULL) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(cvg_stop, aes(-nt_locus_coord_rev, cvg_by_pos)) +
  geom_line(col = 'dodgerblue3') +
  scale_x_continuous(limits = c(-60, 15),
                     breaks = seq(-60, 15, 15),
                     minor_breaks = NULL) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Plot the cumulative coverage around start and stop codons
ggplot(cvg_start, aes(nt_locus_coord, cvg_by_pos)) +
  geom_col(fill = "dodgerblue3",
           color = "black",
           width = 1) +
  scale_x_continuous(limits = c(-15, 24),
                     breaks = seq(-15, 24, 1),
                     minor_breaks = NULL) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("figures/coverage_start_codon.pdf",
       width = 8, height = 8, device = cairo_pdf)

ggplot(cvg_stop, aes(-nt_locus_coord_rev, cvg_by_pos)) +
  geom_col(fill = "dodgerblue3",
           color = "black",
           width = 1) +
  scale_x_continuous(limits = c(-24, 15),
                     breaks = seq(-24, 15, 1),
                     minor_breaks = NULL) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())