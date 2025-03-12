library(tidyverse)

ifelse(!dir.exists("gene_pausing_figures"), dir.create("gene_pausing_figures"),
       "Folder 'gene_pausing_figures' exists already")

cds <- read_csv("riboseq_cds.csv")
cds_dp <- filter(cds, d0316_up | d0316_dw)# cds_dp <- filter(cds, fdr < 1e-04 & !between(logFC, -2, 2) &
#                    score_mean_wt + score_mean_d0316 + score_mean_comp + score_mean_m6  > 4)

# A function to plot the pausing scores of a selected locus
plot_scores <- function(site_description = "SPY_RS06520-SPy_1568-N-46"){
  site_description_v <- str_split_1(site_description, "-")
  pausing_coord <- as.integer(site_description_v[4])
  print(site_description_v)
  cds_sel <- cds |> 
    filter(locus_tag == site_description_v[1]) |>
    mutate(de_pos = if_else(d0316_up | d0316_dw, "*", "", missing = "")) |> 
    unite("labels", c(aa, de_pos), sep = "\n", remove = FALSE)
  
  if (pausing_coord < 20){
    region = c(1, 40)
  } else if (pausing_coord > nrow(cds_sel) - 20){
    region = c(nrow(cds_sel) - 40, nrow(cds_sel))
  } else {
    region = c(pausing_coord - 19, pausing_coord + 20)
  }
  
  cds_sel <- cds_sel |>
    filter(between(aa_locus_coord, region[1], region[2])) |> 
    arrange(aa_locus_coord)|> 
    select(locus_tag, aa_locus_coord, labels, score_mean_wt:score_sd_m6) 
  
  cds_plot <- cds_sel |>
    pivot_longer(
      cols = score_mean_wt:score_sd_m6,
      names_pattern = "(.+_.+)_(.+)",
      names_to = c(".value", "sample"),
    )
  
  cds_plot$sample <- factor(cds_plot$sample, levels = c("wt", "d0316", "comp", "m6"))
 
  ggplot(cds_plot, aes(x = aa_locus_coord, y = score_mean + 0.1)) +
    geom_errorbar(aes(ymin = score_mean, ymax = score_mean + score_sd),
                  linewidth = 0.5,
                  width = 0.4) +
    geom_col(fill = "red",
             color = "black",
             width = 1) +
    facet_wrap(~sample,
               nrow = 4,
               scales = "free_x",
               #strip.position = "right",
               labeller = labeller(sample = 
                                     c("wt" = "WT",
                                       "d0316" = "\u0394yebC",
                                       "comp" = "\u0394yebC + yebC",
                                       "m6" = "\u0394yebC + yebC_M2"))) +
    scale_x_continuous(name = sprintf("%s (%s) (%d..%d)",
                                      site_description_v[1],
                                      site_description_v[2],
                                      region[1],
                                      region[2]), 
                       breaks = seq(region[1], region[2]),
                       labels = cds_sel$labels,
                       minor_breaks = FALSE,
                       expand = c(0, 0)) +
    scale_y_continuous(name = "pausing score",
                       minor_breaks = FALSE,
                       expand = expansion(mult = c(0, 0.05))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #strip.background = element_blank(),
          strip.text = element_text(face = "italic"))
  ggsave(
    filename = paste("gene_pausing_figures/pausing-", site_description, ".pdf", sep = ""),
    width = 6, height = 8, device = cairo_pdf
  )
}
plot_scores()

plot_scores("SPY_RS00060-SPy_0015-T-231")
plot_scores("SPY_RS00460-SPy_0077-P-10")
plot_scores("SPY_RS00715-NA-P-281")
plot_scores("SPY_RS00745-NA-S-68")
plot_scores("SPY_RS00795-SPy_0145-P-7")
plot_scores("SPY_RS00810-SPy_0148-W-341")
plot_scores("SPY_RS00835-SPy_0154-D-397")
plot_scores("SPY_RS00835-SPy_0154-I-446")
plot_scores("SPY_RS00840-SPy_0155-I-317")
plot_scores("SPY_RS01065-SPy_0212-P-103")
plot_scores("SPY_RS01215-SPy_0259-F-264")
plot_scores("SPY_RS01265-SPy_0269-D-329")
plot_scores("SPY_RS01335-SPy_0288-A-30")
plot_scores("SPY_RS01360-NA-T-589")
plot_scores("SPY_RS01375-SPy_0296-P-308")
plot_scores("SPY_RS01435-SPy_0306-T-24")
plot_scores("SPY_RS01460-SPy_0312-F-125")
plot_scores("SPY_RS01730-SPy_0385-V-58")
plot_scores("SPY_RS01930-SPy_0453-P-67")
plot_scores("SPY_RS01950-SPy_0458-G-421")
plot_scores("SPY_RS01960-SPy_0460-A-18")
plot_scores("SPY_RS02160-SPy_0510-E-103")
plot_scores("SPY_RS02305-SPy_0542-E-122")
plot_scores("SPY_RS02310-SPy_0543-A-310")
plot_scores("SPY_RS02355-SPy_0553-I-23")
plot_scores("SPY_RS02720-SPy_0653-D-56")
plot_scores("SPY_RS03075-SPy_0739-P-16")
plot_scores("SPY_RS03105-SPy_0745-P-348")
plot_scores("SPY_RS03155-SPy_0758-A-455")
plot_scores("SPY_RS03165-SPy_0760-M-210")
plot_scores("SPY_RS03230-SPy_0776-G-106")
plot_scores("SPY_RS03310-SPy_0794-D-181")
plot_scores("SPY_RS03415-SPy_0817-I-204")
plot_scores("SPY_RS03450-SPy_0827-S-186")
plot_scores("SPY_RS03460-SPy_0831-P-293")
plot_scores("SPY_RS03475-SPy_0835-V-998")
plot_scores("SPY_RS03725-SPy_0892-W-53")
plot_scores("SPY_RS03745-SPy_0899-D-147")
plot_scores("SPY_RS03795-SPy_0910-N-104")
plot_scores("SPY_RS03945-SPy_0939-D-232")
plot_scores("SPY_RS04285-SPy_1029-P-109")
plot_scores("SPY_RS04520-SPy_1082-A-440")
plot_scores("SPY_RS04605-SPy_1101-H-225")
plot_scores("SPY_RS04845-SPy_1155-S-62")
plot_scores("SPY_RS04935-SPy_1177-I-153")
plot_scores("SPY_RS05220-SPy_1246-W-56")
plot_scores("SPY_RS05385-SPy_1286-V-158")
plot_scores("SPY_RS05695-SPy_1368-K-92")
plot_scores("SPY_RS05705-SPy_1370-T-256")
plot_scores("SPY_RS05935-SPy_1434-T-231")
plot_scores("SPY_RS06270-SPy_1509-E-546")
plot_scores("SPY_RS06435-SPy_1547-A-50")
plot_scores("SPY_RS06445-SPy_1549-F-153")
plot_scores("SPY_RS06490-SPy_1562-F-339")
plot_scores("SPY_RS06520-SPy_1568-N-46")
plot_scores("SPY_RS06580-SPy_1582-L-141")
plot_scores("SPY_RS06800-SPy_1629-Y-698")
plot_scores("SPY_RS06825-SPy_1637-A-357")
plot_scores("SPY_RS06835-SPy_1639-F-91")
plot_scores("SPY_RS06990-SPy_1678-T-138")
plot_scores("SPY_RS07085-SPy_1704-L-164")
plot_scores("SPY_RS09415-NA-L-65")
plot_scores("SPY_RS07355-SPy_1765-E-162")
plot_scores("SPY_RS07375-SPy_1770-L-279")
plot_scores("SPY_RS07415-SPy_1781-Q-51")
plot_scores("SPY_RS07575-SPy_1827-Y-76")
plot_scores("SPY_RS07665-SPy_1849-F-89")
plot_scores("SPY_RS07985-SPy_1916-Y-180")
plot_scores("SPY_RS07990-SPy_1917-E-270")
plot_scores("SPY_RS08000-SPy_1919-L-163")
plot_scores("SPY_RS08200-SPy_1973-H-215")
plot_scores("SPY_RS08275-SPy_1992-I-50")
plot_scores("SPY_RS08300-SPy_1998-Q-144")
plot_scores("SPY_RS08375-SPy_2019-P-526")
plot_scores("SPY_RS08420-SPy_2034-L-68")
plot_scores("SPY_RS08565-SPy_2073-Y-604")
plot_scores("SPY_RS08660-SPy_2096-R-57")
plot_scores("SPY_RS08710-SPy_2111-R-22")
plot_scores("SPY_RS08975-SPy_2174-E-244")
plot_scores("SPY_RS09045-SPy_2191-S-84")
plot_scores("SPY_RS09135-SPy_2211-P-568")


########## End ##########
