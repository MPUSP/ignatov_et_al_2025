library(tidyverse)

ifelse(!dir.exists("figures"), dir.create("figures"), "Folder 'figures' exists already")

iclip_df <- read.csv('iclip.csv')
pl <- c("grey", '#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f')

# Plot the 23S rRNA coverage
cvg_23S_df <- iclip_df |> 
  filter(locus_tag == 'SPY_RS00080') |> 
  select(nt_locus_coord,
         cvg_yebc_uv_1_plus,
         cvg_smi_1_plus,
         cvg_wt_uv_plus,
         cvg_yebc_no_uv_plus) |> 
  pivot_longer(cols = -1, names_to = 'sample_name', values_to = 'rpm') |> 
  mutate(sample_name = factor(sample_name, levels = c('cvg_yebc_uv_1_plus',
                                                      'cvg_smi_1_plus',
                                                      'cvg_wt_uv_plus',
                                                      'cvg_yebc_no_uv_plus')))

ggplot(cvg_23S_df, aes(x = nt_locus_coord, y = rpm, color = sample_name)) +
  geom_line() +
  facet_wrap(~sample_name, nrow = 4) +
  theme_bw() + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
ggsave('figures/23S_coverage.pdf', height = 4, width = 10)

cl_23S_df <- iclip_df |> 
  filter(locus_tag == 'SPY_RS00080') |> 
  select(nt_locus_coord,
         cl_yebc_uv_1_plus,
         cl_smi_1_plus,
         cl_wt_uv_plus,
         cl_yebc_no_uv_plus) |> 
  pivot_longer(cols = -1, names_to = 'sample_name', values_to = 'rpm') |> 
  mutate(sample_name = factor(sample_name, levels = c('cl_yebc_uv_1_plus',
                                                      'cl_smi_1_plus',
                                                      'cl_wt_uv_plus',
                                                      'cl_yebc_no_uv_plus')))

ggplot(cl_23S_df, aes(x = nt_locus_coord, y = rpm, color = sample_name)) +
  geom_col() +
  xlim(c(2450, 2470)) +
  facet_wrap(~sample_name, nrow = 4) +
  theme_bw() + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

# Count the number of CL sites by category
iclip_df_mod <- iclip_df |> 
  rowwise() |> 
  mutate(
    cl_yebc_uv_1 = sum(cl_yebc_uv_1_plus, cl_yebc_uv_1_minus),
    cl_yebc_uv_2 = sum(cl_yebc_uv_2_plus, cl_yebc_uv_2_minus),
    cl_smi_1 = sum(cl_smi_1_plus, cl_smi_1_minus),
    cl_smi_2 = sum(cl_smi_2_plus, cl_smi_2_minus),
    cl_yebc_no_uv = sum(cl_yebc_no_uv_plus, cl_yebc_no_uv_minus),
    cl_wt_uv = sum(cl_wt_uv_plus, cl_wt_uv_minus)
  )

cl_biotype <- iclip_df_mod |> 
  select(
    biotype,
    cl_yebc_uv_1,
    cl_yebc_uv_2,
    cl_smi_1,
    cl_smi_2,
    cl_yebc_no_uv,
    cl_wt_uv
  ) |> 
  group_by(biotype) |> 
  summarise(
    cl_yebc_uv_1 = sum(cl_yebc_uv_1),
    cl_yebc_uv_2 = sum(cl_yebc_uv_2),
    cl_smi_1 = sum(cl_smi_1),
    cl_smi_2 = sum(cl_smi_2),
    cl_yebc_no_uv = sum(cl_yebc_no_uv),
    cl_wt_uv = sum(cl_wt_uv)
  ) |>
  pivot_longer(!biotype, names_to = 'sample_name', values_to = 'cl_sum')

# Data frame to store the clusters
cluster_df <- data.frame(
  sample_name = character(),
  locus_tag = character(),
  biotype = character(),
  start = numeric(),
  end = numeric(),
  cl_sum = numeric()
)

# Iterate over the CL sites for each sample and strand

for(sample_name in c(
  'cl_yebc_uv_1_plus', 'cl_yebc_uv_1_minus',
  'cl_yebc_uv_2_plus', 'cl_yebc_uv_2_minus',
  'cl_smi_1_plus', 'cl_smi_1_minus',
  'cl_smi_2_plus', 'cl_smi_2_minus',
  'cl_yebc_no_uv_plus', 'cl_yebc_no_uv_minus',
  'cl_wt_uv_plus', 'cl_wt_uv_minus'
)){
  # Map the clusters using the sliding window approach 
  cl_dens <- iclip_df[, sample_name]
  cluster_dens <- rep(0, length(cl_dens))
  for(i in 1:(length(cl_dens) - 90)){
    up_dens = cl_dens[i:(i+30)]
    up_sum = sum(up_dens)
    win_dens = cl_dens[(i+31):(i+60)]
    win_sum = sum(win_dens)
    dw_dens = cl_dens[(i+61):(i+90)]
    dw_sum = sum(dw_dens)
    if(win_sum / (up_sum + 1E-10) > 5 & win_sum / (dw_sum + 1E-10) > 5 & win_sum > 3){
      cluster_dens[i+46] = win_sum
    }
  }
  
  # Refine the clusters and submit them to the data frame
  cluster_forming = FALSE
  cluster_coord = numeric()
  
  for(i in 1:length(cluster_dens)){
    if(cluster_dens[i] == 0){
      if(cluster_forming){
        cluster_center = round(mean(cluster_coord))
        cluster_df <- cluster_df |> 
          add_row(
            sample_name = sample_name,
            locus_tag = iclip_df[cluster_center, 'locus_tag'],
            biotype = iclip_df[cluster_center, 'biotype'],
            start = cluster_center - 15,
            end = cluster_center + 15,
            cl_sum = sum(cl_dens[start:end])
          )
        cluster_coord = numeric()
        cluster_forming = FALSE
      }
    } else {
      cluster_forming = TRUE
      cluster_coord = c(cluster_coord, i)
    }
  }
}

# Combine clusters on plus and minus strands and modify the data frame
cluster_df <- cluster_df |> 
  mutate(
    strand = if_else(str_detect(sample_name, 'plus'), '+', '-'),
    #sample_name = str_replace(sample_name, 'cl_', ''),
    sample_name = str_replace(sample_name, '_plus', ''),
    sample_name = str_replace(sample_name, '_minus', ''),
    cluster_name = sprintf("%s_(%d..%d)%s", locus_tag, start, end, strand),
    cluster_name = str_replace(cluster_name, 'NA_', 'ncRNA_')
  ) |> 
  relocate(cluster_name, .before = sample_name) |> 
  relocate(strand, .before = start)

# Save the table with the clusters of CL sites
write_csv(cluster_df, file = "cl_clusters.csv")
write_tsv(cluster_df, file = "cl_clusters.tsv")

# Select top clusters
top_cluster_df <- cluster_df |>
  group_by(sample_name) |> 
  slice_max(cl_sum, n = 1) |> 
  ungroup() |> 
  select(
    sample_name,
    cl_sum
  ) |> 
  mutate(biotype = 'top_cluster')

# Add top clusters to the CL statistics and modify it
cl_stat <- cl_biotype |> 
  bind_rows(top_cluster_df)

cl_stat <- mutate(cl_stat, n = 1:nrow(cl_stat))

cl_stat[37, 'cl_sum'] = cl_stat[37, 'cl_sum'] - cl_stat[59, 'cl_sum']
cl_stat[38, 'cl_sum'] = cl_stat[38, 'cl_sum'] - cl_stat[60, 'cl_sum']
cl_stat[39, 'cl_sum'] = cl_stat[39, 'cl_sum'] - cl_stat[55, 'cl_sum']
cl_stat[40, 'cl_sum'] = cl_stat[40, 'cl_sum'] - cl_stat[56, 'cl_sum']
cl_stat[41, 'cl_sum'] = cl_stat[41, 'cl_sum'] - cl_stat[58, 'cl_sum']
cl_stat[42, 'cl_sum'] = cl_stat[42, 'cl_sum'] - cl_stat[57, 'cl_sum']

cl_stat <- cl_stat |> 
  filter(
    biotype %in% c(
      'top_cluster',
      'rRNA',
      'tRNA',
      'protein_coding',
      'intergenic'
    )
  ) |> 
  mutate(
    biotype = factor(biotype, levels = c('top_cluster',
                                         'rRNA',
                                         'tRNA',
                                         'protein_coding',
                                         'intergenic')),
    sample_name = factor(sample_name, levels = rev(c('cl_yebc_uv_1',
                                                 'cl_yebc_uv_2',
                                                 'cl_smi_1',
                                                 'cl_smi_2',
                                                 'cl_wt_uv',
                                                 'cl_yebc_no_uv')))
  )

ggplot(cl_stat, aes(x = sample_name, y = cl_sum, fill = biotype)) +
  geom_col(col = 'black') +
  scale_fill_manual(values=c('#1B9E77', '#66c2a5','#fc8d62','#8da0cb','#e78ac3')) +
  coord_flip() +
  theme_void()

ggsave('figures/barplot_cl_stat.pdf', width = 7, height = 5)

########## End ##########
