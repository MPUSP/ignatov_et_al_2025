library(tidyverse)
library(edgeR)

rsd <- read_csv("riboseq.csv")

rsd <- filter(rsd, locus_tag != "SPY_RS01470")

cds <- data.frame(
  locus_tag = character(),
  name = character(),
  old_locus_tag = character(),
  locus_start = numeric(),
  locus_stop = numeric(),
  locus_strand = character(),
  aa = character(),
  aa_locus_coord = numeric(),
  
  cvg_wt_1 = numeric(),
  cvg_wt_2 = numeric(),
  cvg_wt_3 = numeric(),
  cvg_d0316_1 = numeric(),
  cvg_d0316_2 = numeric(),
  cvg_d0316_3 = numeric(),
  cvg_comp_1 = numeric(),
  cvg_comp_2 = numeric(),
  cvg_comp_3 = numeric(),
  cvg_m6_1 = numeric(),
  cvg_m6_2 = numeric(),
  cvg_m6_3 = numeric(),
  
  score_wt_1 = numeric(),
  score_wt_2 = numeric(),
  score_wt_3 = numeric(),
  score_d0316_1 = numeric(),
  score_d0316_2 = numeric(),
  score_d0316_3 = numeric(),
  score_comp_1 = numeric(),
  score_comp_2 = numeric(),
  score_comp_3 = numeric(),
  score_m6_1 = numeric(),
  score_m6_2 = numeric(),
  score_m6_3 = numeric(),
  
  score_mean_wt = numeric(),
  score_mean_d0316 = numeric(),
  score_mean_comp = numeric(),
  score_mean_m6 = numeric(),
  
  score_sd_wt = numeric(),
  score_sd_d0316 = numeric(),
  score_sd_comp = numeric(),
  score_sd_m6 = numeric(),

  logFC = numeric(),
  PValue = numeric()
)

# Iterate over CDSs
locus_tag_ids = unique(rsd$locus_tag)

for(i in 1:length(locus_tag_ids)){
  locus_tag_id = locus_tag_ids[i]
  print(locus_tag_id)
  
  rsd_sel <- rsd |> 
    filter(locus_type == "CDS" & locus_tag == locus_tag_id)
  
  cds_sel <- rsd_sel|> 
    #filter(between(nt_locus_coord, 16, max(nt_locus_coord) - 15)) |> 
    group_by(aa_locus_coord) |> 
    summarize(
      cvg_wt_1 = sum(cvg_wt_1),
      cvg_wt_2 = sum(cvg_wt_2),
      cvg_wt_3 = sum(cvg_wt_3),
      cvg_d0316_1 = sum(cvg_d0316_1),
      cvg_d0316_2 = sum(cvg_d0316_2),
      cvg_d0316_3 = sum(cvg_d0316_3),
      cvg_comp_1 = sum(cvg_comp_1),
      cvg_comp_2 = sum(cvg_comp_2),
      cvg_comp_3 = sum(cvg_comp_3),
      cvg_m6_1 = sum(cvg_m6_1),
      cvg_m6_2 = sum(cvg_m6_2),
      cvg_m6_3 = sum(cvg_m6_3)
    )
  
  m <- as.matrix(cds_sel[,-1])
  row.names(m) <- cds_sel$aa_locus_coord
  
  # Test for the expression level
  if(sum(m)/nrow(m) >= 4){
    
    # Calculate pausing scores
    cds_sel <- cds_sel |> 
      mutate(
        score_wt_1 = cvg_wt_1 * nrow(cds_sel) / sum(cvg_wt_1),
        score_wt_2 = cvg_wt_2 * nrow(cds_sel) / sum(cvg_wt_2),
        score_wt_3 = cvg_wt_3 * nrow(cds_sel) / sum(cvg_wt_3),
        score_d0316_1 = cvg_d0316_1 * nrow(cds_sel) / sum(cvg_d0316_1),
        score_d0316_2 = cvg_d0316_2 * nrow(cds_sel) / sum(cvg_d0316_2),
        score_d0316_3 = cvg_d0316_3 * nrow(cds_sel) / sum(cvg_d0316_3),
        score_comp_1 = cvg_comp_1 * nrow(cds_sel) / sum(cvg_comp_1),
        score_comp_2 = cvg_comp_2 * nrow(cds_sel) / sum(cvg_comp_2),
        score_comp_3 = cvg_comp_3 * nrow(cds_sel) / sum(cvg_comp_3),
        score_m6_1 = cvg_m6_1 * nrow(cds_sel) / sum(cvg_m6_1),
        score_m6_2 = cvg_m6_2 * nrow(cds_sel) / sum(cvg_m6_2),
        score_m6_3 = cvg_m6_3 * nrow(cds_sel) / sum(cvg_m6_3)
      ) |> rowwise() |> 
      mutate(
        score_mean_wt = mean(c(score_wt_1, score_wt_2, score_wt_3)),
        score_sd_wt = sd(c(score_wt_1, score_wt_2, score_wt_3)),
        score_mean_d0316 = mean(c(score_d0316_1, score_d0316_2, score_d0316_3)),
        score_sd_d0316 = sd(c(score_d0316_1, score_d0316_2, score_d0316_3)),
        score_mean_comp = mean(c(score_comp_1, score_comp_2, score_comp_3)),
        score_sd_comp = sd(c(score_comp_1, score_comp_2, score_comp_3)),
        score_mean_m6 = mean(c(score_m6_1, score_m6_2, score_m6_3)),
        score_sd_m6 = sd(c(score_m6_1, score_m6_2, score_m6_3))
      ) |> ungroup()
    
    # Use edgeR to find DE positions
    sample_groups <- rep(c("W", "D", "W", "D"), each = 3)
    d <- DGEList(counts = m, group = sample_groups)
    d <- calcNormFactors(d)
    d <- estimateCommonDisp(d)
    d <- estimateTagwiseDisp(d)
    et <- exactTest(d, pair=c(2,1))
    cds_sel$logFC <- et$table$logFC
    cds_sel$PValue <- et$table$PValue
    
  }else{
    
    # For the lowly expressed genes do not do further calculations
    cds_sel <- cds_sel |> 
      mutate(
        score_wt_1 = rep(NA, nrow(cds_sel)),
        score_wt_2 = rep(NA, nrow(cds_sel)),
        score_wt_3 = rep(NA, nrow(cds_sel)),
        score_d0316_1 = rep(NA, nrow(cds_sel)),
        score_d0316_2 = rep(NA, nrow(cds_sel)),
        score_d0316_3 = rep(NA, nrow(cds_sel)),
        score_comp_1 = rep(NA, nrow(cds_sel)),
        score_comp_2 = rep(NA, nrow(cds_sel)),
        score_comp_3 = rep(NA, nrow(cds_sel)),
        score_m6_1 = rep(NA, nrow(cds_sel)),
        score_m6_2 = rep(NA, nrow(cds_sel)),
        score_m6_3 = rep(NA, nrow(cds_sel)),
        
        score_mean_wt = rep(NA, nrow(cds_sel)),
        score_mean_d0316 = rep(NA, nrow(cds_sel)),
        score_mean_comp = rep(NA, nrow(cds_sel)),
        score_mean_m6 = rep(NA, nrow(cds_sel)),
        score_sd_wt = rep(NA, nrow(cds_sel)),
        score_sd_d0316 = rep(NA, nrow(cds_sel)),
        score_sd_comp = rep(NA, nrow(cds_sel)),
        score_sd_m6 = rep(NA, nrow(cds_sel)),
        
        logFC = rep(NA, nrow(cds_sel)),
        PValue = rep(NA, nrow(cds_sel))
      )
  }
  
  cds_sel <- rsd_sel |>
    distinct(aa_locus_coord, .keep_all = TRUE) |> 
    select(locus_tag, name, old_locus_tag, locus_start, locus_stop, locus_strand,
           aa, aa_locus_coord) |> 
    left_join(cds_sel, join_by(aa_locus_coord))
  
  cds <- rbind(cds, cds_sel)
}

cds <- cds |> 
  mutate(
    fdr = p.adjust(PValue, method = "fdr"),
    d0316_up = fdr < 1e-03 & logFC > 2 & score_mean_d0316 + score_mean_m6 > 4,
    d0316_dw = fdr < 1e-03 & logFC < -2 & score_mean_wt + score_mean_comp > 4
  ) |> 
  unite("site_description", locus_tag, old_locus_tag, aa, aa_locus_coord,
         sep = "-", remove = FALSE)

cds_diff <- filter(cds, d0316_up | d0316_dw)
write_csv(cds, file = "riboseq_cds.csv")
write_csv(cds_diff, file = "riboseq_cds_diff.csv")
write_tsv(cds_diff, file = "riboseq_cds_diff.tsv")

########## End ##########
