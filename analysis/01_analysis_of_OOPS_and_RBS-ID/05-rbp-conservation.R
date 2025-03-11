library(janitor)
library(tidyverse)

pal_set2 <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
pal_dark_2 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")

rbp <- read_tsv("Raw data/241120_RAB_spy_ignatov_other_studies.tsv")

rbp[rbp == "-"] <- NA
rbp <- clean_names(rbp)
rbp <- rbp |> rename(uniprot = uni_prot)

# Add different definitions of RBPs
rbp <- rbp |>
  mutate(
    rbp_oops = if_else(oops_p_adj < 0.05, 'RBP', 'no_RBP'),
    rbp_gt2 = if_else(rbp_oops == 'RBP' & oops_mean > 1, 'RBP', 'no_RBP'),
    rbp_rbs = if_else(rbp_gt2 == 'RBP' & rna_binding_sites != '-', 'RBP', 'no_RBP'),
    
    rbp_oops = replace_na(rbp_oops, 'no_RBP'),
    rbp_gt2 = replace_na(rbp_gt2, 'no_RBP'),
    rbp_rbs = replace_na(rbp_rbs, 'no_RBP')
  )

# Export OOPS and RBS-ID results with information from previous studies
oops_rbsid <- rbp |> 
  dplyr::select(
    locus_tag = old_locus_tag,
    gene_name,
    uniprot,
    pfam,
    interpro,
    annot_rbp_type = rbp_type,
    oops_rep1:rna_binding_sites,
    rbp = rbp_rbs,
    novel_rbp,
    eco_shchepachev:sau_chu
  )

write.table(oops_rbsid, "OOPS_RBSID_orthologs.tsv", row.names = FALSE,
            quote = FALSE, sep = "\t")

# Plot the number of RBPs identified in different studies

rbp_annot <-  rbp |> 
  dplyr::select(locus_tag, eco_shchepachev:sau_chu, rbp_rbs) |> 
  pivot_longer(
    cols = eco_shchepachev:rbp_rbs,
    names_to = 'study',
    values_to = 'annotation'
  ) |> 
  mutate(
    study = factor(study, levels = c('rbp_rbs',
                                     'eco_queiroz',
                                     'eco_shchepachev',
                                     'eco_monti',
                                     'eco_stenum',
                                     'sal_urdaneta',
                                     'sau_chu')),
    annotation = factor(annotation, levels = c('no_homolog',
                                               'no_RBP',
                                               'RBP',
                                               NA))
  )

ggplot(dplyr::select(rbp_annot, 2:3), aes(x = study, fill = annotation)) +
  geom_bar(color = 'black') +
  scale_fill_manual(values = c("white", "#B3B3B3", "#8DA0CB")) +
  scale_y_continuous(expand = c(0, 150)) +
  theme_light()
ggsave("Figures/Number of RBPs in different studies.pdf")

annot_stat <- dplyr::select(rbp_annot, 2:3) |> 
  count(annotation, study)

# Plot the number of RBPs from this and previous studies that identified the annotated RBPs

an_rbp <- rbp |> 
  mutate(
    an_oops = if_else(rbp_oops == 'RBP', 'RBP', NA),
    an_oops = if_else(!is.na(rbp_type) & an_oops == 'RBP', 'an_RBP', an_oops),
    
    an_gt2 = if_else(rbp_gt2 == 'RBP', 'RBP', NA),
    an_gt2 = if_else(!is.na(rbp_type) & an_gt2 == 'RBP', 'an_RBP', an_gt2),
    
    an_rbs = if_else(rbp_rbs == 'RBP', 'RBP', NA),
    an_rbs = if_else(!is.na(rbp_type) & an_rbs == 'RBP', 'an_RBP', an_rbs),
    
    an_eco_shchepachev = if_else(eco_shchepachev == 'RBP', 'RBP', NA),
    an_eco_shchepachev = if_else(!is.na(rbp_type) & eco_shchepachev == 'RBP', 'an_RBP', an_eco_shchepachev),
    
    an_eco_queiroz = if_else(eco_queiroz == 'RBP', 'RBP', NA),
    an_eco_queiroz = if_else(!is.na(rbp_type) & eco_shchepachev == 'RBP', 'an_RBP', an_eco_queiroz),
    
    an_eco_monti = if_else(eco_monti == 'RBP', 'RBP', NA),
    an_eco_monti = if_else(!is.na(rbp_type) & eco_monti == 'RBP', 'an_RBP', an_eco_monti),
    
    an_eco_stenum = if_else(eco_stenum == 'RBP', 'RBP', NA),
    an_eco_stenum = if_else(rbp_gt2 == 'RBP' & eco_stenum == 'RBP', 'an_RBP', an_eco_stenum),
    
    an_sal_urdaneta = if_else(sal_urdaneta == 'RBP', 'RBP', NA),
    an_sal_urdaneta = if_else(!is.na(rbp_type) & sal_urdaneta == 'RBP', 'an_RBP', an_sal_urdaneta),
    
    an_sau_chu = if_else(sau_chu == 'RBP', 'RBP', NA),
    an_sau_chu = if_else(!is.na(rbp_type) & sau_chu == 'RBP', 'an_RBP', an_sau_chu),
  )

an_rbp_annot <-  an_rbp |> 
  dplyr::select(locus_tag, an_oops:an_sau_chu) |> 
  pivot_longer(
    cols = an_oops:an_sau_chu,
    names_to = 'study',
    values_to = 'annotation'
  ) |> 
  mutate(
    study = factor(study, levels = c('an_rbs',
                                     'an_eco_queiroz',
                                     'an_eco_shchepachev',
                                     'an_eco_monti',
                                     'an_eco_stenum',
                                     'an_sal_urdaneta',
                                     'an_sau_chu')),
    annotation = factor(annotation, levels = c('RBP',
                                               'an_RBP',
                                               NA))
  )

an_annot_stat <- dplyr::select(an_rbp_annot, 2:3) |> 
  count(annotation, study) |> 
  filter(!is.na(annotation))

ggplot(an_annot_stat, aes(x = study, y = n, fill = annotation)) +
  geom_col(color = 'black') +
  scale_fill_manual(values = c("#c2cce3", "#5874b3")) +
  scale_y_continuous(expand = c(0, 30)) +
  theme_light()

ggsave("Figures/Annotated RBPs in different studies.pdf")


####### End ########
