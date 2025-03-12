library(tidyverse)
library(heatmaply)
library(Biostrings)
library(ggseqlogo)

cds <- read_csv("riboseq_cds.csv")

ifelse(!dir.exists("gene_pausing_figures"), dir.create("gene_pausing_figures"),
       "Folder 'gene_pausing_figures' exists already")
ifelse(!dir.exists("figures"), dir.create("figures"),
       "Folder 'figures' exists already")

pl_dark2 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D")

# Filter out the lowly expressed genes
cds_sel <- cds |> 
  rowwise() |>
  mutate(cvg_average = mean(cvg_wt_1:cvg_m6_3)) |>
  group_by(locus_tag) |> 
  mutate(cvg_average_gene = mean(cvg_average)) |>
  ungroup() |> 
  filter(cvg_average_gene >= 10) |> 
  mutate(cvg_average = NULL, cvg_average_gene = NULL)

# Filter out the first and the last 10 amino acids
cds_sel <- cds_sel |> 
  group_by(locus_tag) |> 
  mutate(cds_length = max(aa_locus_coord)) |> 
  ungroup() |> 
  filter(between(aa_locus_coord, 10, cds_length - 10))

# Extract aa sequences of proteins
genome_seq <- readDNAStringSet("../../data/genome/NC_002737.2.fa")[[1]]

extract_aa <- function(genome_seq, locus_start, locus_stop, locus_strand){
  dna_seq = subseq(genome_seq, locus_start, locus_stop)
  if(locus_strand == '-'){dna_seq = reverseComplement(dna_seq)}
  as.character(translate(dna_seq))
}

protein_df <- cds |> 
  group_by(locus_tag) |> 
  slice_head(n = 1) |> 
  select(locus_tag:locus_strand) |> 
  mutate(prot_seq = extract_aa(genome_seq, locus_start, locus_stop, locus_strand))

protein_df <- as.data.frame(protein_df)

########## Consensus around the differential pausing sites ##########

# Select the pausing sites up-regulated in d0316 and m6
cds_d0316_up <- filter(cds, d0316_up)

# Extract aa sequences around the pausing sites (proteins should be longer than 20 aa)
pausing_seq_df <- cds_d0316_up |> 
  select(!cvg_wt_1:score_m6_3) |> 
  left_join(protein_df) |> 
  arrange(locus_tag, aa_locus_coord)
pausing_seq_df <- as.data.frame(pausing_seq_df)

extract_pausing_aa <- function(aa_locus_coord, prot_seq){
  aa_start = aa_locus_coord - 9
  aa_end = aa_locus_coord + 10
  prot_end = nchar(prot_seq) - 1
  if(aa_start < 1){
    pausing_seq = strrep(" ", 1 - aa_start)
    pausing_seq = paste(pausing_seq, subseq(prot_seq, 1, aa_end), sep = "")
  } else if(aa_end > prot_end){
    pausing_seq = strrep(" ", aa_end - prot_end)
    pausing_seq = paste(subseq(prot_seq, aa_start, prot_end), pausing_seq, sep ="")
  } else {
    pausing_seq = subseq(prot_seq, aa_start, aa_end)
  }
  pausing_seq
}

pausing_seq_vector = character()
for(i in 1:nrow(pausing_seq_df)){
  pausing_seq = extract_pausing_aa(pausing_seq_df[i, "aa_locus_coord"],
                                   pausing_seq_df[i, "prot_seq"])
  pausing_seq_vector = c(pausing_seq_vector, pausing_seq)
}

pausing_seq_df$pausing_seq = pausing_seq_vector

# Create the WebLogo plots and export the neighboring regions as FASTA 
ggseqlogo(pausing_seq_vector)
ggsave("figures/weblogo_bits.pdf")

ggseqlogo(pausing_seq_vector, method = 'prob')
ggsave("figures/weblogo_prob.pdf")

names_vector <- paste(1:nrow(pausing_seq_df), pausing_seq_df$site_description, sep = "_")
names_vector <- paste(">", names_vector, sep = "")
write_lines(paste(names_vector, pausing_seq_vector, sep = "\n"),
           "pausing_seq.fa")

########## Motifs at the clustered pausing sites ##########

# Cluster the neighboring pausing sites
pausing_cluster_df <- data.frame(
  locus_tag = character(),
  name = character(),
  old_locus_tag = character(),
  aa = character(),
  pausing_seq = character()
)

for(locus_tag_i in unique(pausing_seq_df$locus_tag)){
  pausing_seq_df_i =  filter(pausing_seq_df, locus_tag == locus_tag_i)
  aa_i = paste0(pausing_seq_df_i$aa, pausing_seq_df_i$aa_locus_coord, sep = "", collapse = ", ")
  pausing_cluster_df_i <- data.frame(
    locus_tag = locus_tag_i,
    name = pausing_seq_df_i[1, "name"],
    old_locus_tag = pausing_seq_df_i[1, "old_locus_tag"],
    aa = paste0(pausing_seq_df_i$aa, pausing_seq_df_i$aa_locus_coord, sep = "", collapse = ", "),
    pausing_seq = pausing_seq_df_i[1, "pausing_seq"]
  )
  pausing_cluster_df = rbind(pausing_cluster_df, pausing_cluster_df_i)
}

# Search for the D.P, P.P and PP. motifs
pausing_cluster_df <- pausing_cluster_df |> 
  mutate(
    pausing_motif = case_when(
      str_detect(pausing_seq, "PPG") ~ "PPG",
      str_detect(pausing_seq, "PP.") ~ "PPX",
      str_detect(pausing_seq, "PIP") ~ "PIP",
      str_detect(pausing_seq, "P.P") ~ "PXP",
      str_detect(pausing_seq, "DIP") ~ "DIP",
      str_detect(pausing_seq, "D.P") ~ "DXP",
      .default = NA
    )
  )

# Save the pausing sites
write_tsv(pausing_cluster_df, 'pausing_clusters.tsv')

# Plot the number of the motifs in the clustered pausing sites
pausing_motif <- pausing_cluster_df |> 
  dplyr::group_by(pausing_motif) |> 
  dplyr::count() |> 
  mutate(
    pausing_motif = factor(pausing_motif,
                              levels = c("PPG", "PPX", "PIP", "PXP", "DIP", "DXP", "NA"))
  )

ggplot(pausing_motif, aes(x="", y=n, fill=pausing_motif)) +
  geom_bar(stat="identity", width=1, color = 'black') +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c('#1B9E77','#66c2a5','#7570B3','#8da0cb','#D95F02','#fc8d62','grey')) +
  coord_polar("y", start=0) +
  theme_void()
ggsave('figures/pausing_motifs_pie_chart.pdf', width = 6, height = 5)

########## Plot the pausing score for each amino acid ##########

aa_score <- cds_sel |> 
  select(aa, score_mean_wt, score_mean_d0316, score_mean_comp, score_mean_m6) |> 
  pivot_longer(
    cols = starts_with("score_mean"),
    names_to = "strain",
    values_to = "score_mean"
  ) |> 
  mutate(
    strain = factor(strain, levels = c("score_mean_wt", "score_mean_d0316", "score_mean_comp", "score_mean_m6"))
  )

ggplot(aa_score, aes(x = aa, y = score_mean, fill = strain)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#66c2a5','#fc8d62','#8da0cb','#e78ac3')) +
  scale_y_continuous(limits = c(0, 5)) +
  theme_bw()
ggsave('figures/scores_amino_acids.pdf', height = 6, width = 12)

########## Average pausing scores around PPX, PXP, DXP etc. ##########

# Collect data for different motifs
seq_score <- data.frame(
  locus_tag = character(),
  old_locus_tag = character(),
  pattern_seq = factor(),
  pattern_exact_seq = factor(),
  coord_rel = numeric(),
  coord_abs = numeric()
)

patterns <- c("PP.", "P.P", "D.P") # The patterns should consist of 3 amino acids

for (pattern_seq in patterns){
  protein_m_df <- protein_df |> 
    filter(str_detect(prot_seq, pattern_seq))
  
  for(i in 1:nrow(protein_m_df)){
    prot_seq = protein_m_df[i, "prot_seq"]
    m = str_locate_all(prot_seq, pattern_seq)[[1]]
    for(j in 1:nrow(m)){
      pattern_exact_seq = str_sub(prot_seq, m[j, 1], m[j, 2])
      if(m[j, 1] > 14 & m[j, 2] < nchar(prot_seq) - 16){
        seq_score_i <- data.frame(
          locus_tag = rep(protein_m_df[i, "locus_tag"], 21),
          old_locus_tag = rep(protein_m_df[i, "old_locus_tag"], 21),
          pattern_seq = rep(pattern_seq, 21),
          pattern_exact_seq = rep(pattern_exact_seq, 21),
          coord_rel = -9:11,
          coord_abs = m[j, 1] + -9:11
        )
        seq_score = rbind(seq_score, seq_score_i)
      }
    }
  }
}

seq_score <- left_join(
  seq_score,
  select(cds_sel, locus_tag, aa, aa_locus_coord,
         score_mean_wt, score_mean_d0316, score_mean_comp, score_mean_m6),
  by = c("locus_tag" = "locus_tag", "coord_abs" = "aa_locus_coord")
)

seq_score_long <- seq_score |>
  pivot_longer(
    cols = starts_with("score_mean"),
    names_to = "strain",
    values_to = "score_mean"
  ) |> 
  mutate(
    strain = factor(strain, levels = c("score_mean_wt", "score_mean_d0316", "score_mean_comp", "score_mean_m6")),
    coord_rel = factor(coord_rel)
  )

# Calculate and plot the average scores per position for different motifs
seq_score_sum <- seq_score_long |> 
  group_by(pattern_seq, pattern_exact_seq, coord_rel, strain) |> 
  summarize(
    median_score = median(score_mean, na.rm = TRUE),
    mean_score = mean(score_mean, na.rm = TRUE),
    n = n()
  ) |>
  unite("subtitle", c("pattern_exact_seq", "n"), sep = " - ", remove = FALSE) |> 
  filter(n >= 10)

for(motif in unique(seq_score_sum$pattern_seq)){
  p <- seq_score_sum |> 
    filter(strain %in% c("score_mean_wt", "score_mean_d0316")) |> 
    filter(pattern_seq == motif) |> 
    ggplot(aes(x = coord_rel, y = mean_score, group = strain)) +
    geom_step(aes(color = strain),
              alpha = 0.5,
              #position = position_dodge(width = 0.2),
              position = position_nudge(x = -0.5)
              ) +
    scale_x_discrete(labels = c(9:1, str_split_1(motif, ""), 1:9),
                     expand = c(0, 0)) +
    #scale_y_continuous(limits = c(0, 15)) +
    scale_color_manual(values=c("black", "red")) +
    ylab("Average pausing score") +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    facet_wrap(~subtitle,
               scales = "free_x",
               drop = FALSE,
               ncol = 6)
  print(p)
}

p <- seq_score_sum |> 
  filter(strain %in% c("score_mean_wt", "score_mean_d0316")) |> 
  filter(pattern_exact_seq %in% c("DIP", "PIP", "PPG")) |> 
  ggplot(aes(x = coord_rel, y = mean_score, group = strain)) +
  geom_step(aes(color = strain),
            alpha = 0.5,
            #position = position_dodge(width = 0.2),
            position = position_nudge(x = -0.5)
  ) +
  scale_x_discrete(labels = c(9:1, str_split_1("XXX", ""), 1:9),
                   expand = c(0, 0)) +
  #scale_y_continuous(limits = c(0, 15)) +
  scale_color_manual(values=c("black", "dodgerblue3")) +
  ylab("Average pausing score") +
  theme_bw() +
  theme(legend.position = "bottom",
        #axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~subtitle,
             scales = "free_x",
             drop = FALSE,
             ncol = 6)
print(p)
ggsave('figures/DIP_PIP_PPG_scores.pdf', width = 12, height = 3.5)

# Extract the genes with PPG or PPP motifs
ppg_genes <- seq_score |> 
  filter(
    pattern_exact_seq %in% c('PPP', 'PPG')
  ) |> 
  dplyr::count(old_locus_tag)
ppg_genes <- unique(filter(seq_score, pattern_exact_seq %in% c('PPP', 'PPG'))$old_locus_tag)

########## End ##########
