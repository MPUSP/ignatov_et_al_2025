library(tidyverse)
library(VennDiagram)
library(beeswarm)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(plotROC)

oops <- read.delim("Raw data/200519_20C011_results.tsv", skip = 1)
clp <- read.delim("RBS-ID - CL proteins.tsv")
cls <- read.delim("RBS-ID - CL sites.tsv")
pa <- read.delim("Proteome annotation with RBP.tsv")

pal_set2 <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
pal_dark_2 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")

# Fix the table with OOPS results

oops$RBP <- NULL
oops$Gene.name <- NULL
oops <- within(oops, {
  Ensembl.gene <- str_replace(Ensembl.gene, "M5005_Spy\\d+;", "")
})
oops[oops$Accession == "P68999", "Ensembl.gene"] <- "Spy_0209"

# Identify significantly enriched proteins in Log and Stat

oops <- within(oops, {
  enriched.log <- p.adj.RBP.log < .05 & mean.ratio.RBP.log > 1
  enriched.log[is.na(enriched.log)] <- FALSE
  enriched.stat <- p.adj.RBP.stat < .05 & mean.ratio.RBP.stat > 1
  enriched.stat[is.na(enriched.stat)] <- FALSE
  #enriched.log_stat <- enriched.log | enriched.stat
})

# Merge Annotation, OOPS and RBS-ID data

rbp <- pa[, c("UniProt", "Pfam", "ENSG", "Gene.name", "Interpro.name",
              "rbp",
              "translation",
              "rnase",
              "rna.modif",
              "ribosome.biogenesis",
              "rna.binding",
              "rna.pol",
              "ribosome",
              "aatrna.synthesis",
              "glycolysis",
              "rbp.type")]

rbp <- merge(rbp,
             oops[, c("Accession",
                      "RBP_log_rep1", "RBP_log_rep2", "RBP_log_rep3", "RBP_log_rep4", 
                      "p.value.RBP.log", "p.adj.RBP.log", "mean.ratio.RBP.log", "enriched.log",
                      "RBP_stat_rep1", "RBP_stat_rep2", "RBP_stat_rep3", "RBP_stat_rep4", 
                      "p.value.RBP.stat", "p.adj.RBP.stat", "mean.ratio.RBP.stat", "enriched.stat")],
             by.x = "UniProt", by.y = "Accession", all.x= TRUE)
rbp$enriched.log[is.na(rbp$enriched.log)] <- FALSE
rbp$enriched.stat[is.na(rbp$enriched.stat)] <- FALSE

clp$CL.identified <- rep(TRUE, nrow(clp))
rbp <- merge(rbp,
             clp[, c("UniProt", "Sequence", "CL.identified", "Rep_1", "Rep_2", "Rep_3",
                     "Modifications")],
             by = "UniProt", all.x= TRUE)
rbp <- within(rbp, {
  CL.identified[is.na(CL.identified)] <- FALSE
  Rep_1[is.na(Rep_1)] <- FALSE
  Rep_2[is.na(Rep_2)] <- FALSE
  Rep_3[is.na(Rep_3)] <- FALSE
})

# Identify candidates to novel RBPs

rbp$novel.rbp <- ifelse(!rbp$rbp & rbp$enriched.log & rbp$mean.ratio.RBP.log > 1 & (rbp$Rep_1 + rbp$Rep_2 + rbp$Rep_3 > 0), 
                        TRUE, FALSE)

# Save the RBP table

write.table(rbp, "Candidate RBPs.tsv", row.names = FALSE,
            quote = FALSE, sep = "\t")

# Calculate the number of CL sites in domains

cls_domain <- cls |> 
  mutate(Func_cat = str_replace(Pfam, " \\(.+", "")) |> 
  count(Func_cat) |> 
  arrange(desc(n))

write.table(cls_domain, "CL sites in domains.tsv", row.names = FALSE,
            quote = FALSE, sep = "\t")

cls_domain_top <- slice(cls_domain, 2:16)

ggplot(cls_domain_top, aes(y = reorder(Func_cat, n), x = n)) + 
  geom_col() +
  scale_x_continuous(expand = c(0, 2)) +
  theme_light()
ggsave("Figures/CL sites in top domains.pdf")

# Calculate the number of CL sites in domain groups

cls_domain_group <- cls |> 
  mutate(Func_cat = str_replace(Pfam, " \\(.+", "")) |>
  mutate(Func_group = "other") |> 
  mutate(
    Func_group = case_when(
      str_detect(Func_cat, regex("ribosom", ignore_case = TRUE)) ~ "ribosome",
      str_detect(Func_cat, regex("RNA")) ~ "RNA",
      str_detect(Func_cat, regex("S1")) ~ "RNA",
      str_detect(Func_cat, regex("KH")) ~ "RNA",
      str_detect(Func_cat, regex("K Homology")) ~ "RNA",
      str_detect(Func_cat, regex("NusB")) ~ "RNA",
      str_detect(Func_cat, regex("Cold-shock protein")) ~ "RNA",
      is.na(Func_cat) ~ "no domain",
      .default = "other"
    )
  ) |> 
  count(Func_group) |> 
  arrange(desc(n))

ggplot(cls_domain_group, aes(y = reorder(Func_group, c(4, 2, 1, 3)), x = n)) + 
  geom_col() +
  scale_x_continuous(expand = c(0, 20)) +
  theme_light()
ggsave("Figures/CL sites in domain groups.pdf")

# Create a dataframe for plotting statistics for RBP types

rbp.sel <- subset(rbp, #!is.na(mean.ratio.RBP.log),
                 select = c(UniProt, enriched.log,
                            mean.ratio.RBP.log, p.adj.RBP.log, CL.identified,
                            rbp, rnase, rna.pol, ribosome, aatrna.synthesis,
                            translation, rna.modif, ribosome.biogenesis,
                            glycolysis))
rbp.rbp <- rbp.sel |> 
  mutate(
    group = ifelse(rbp, "RBP", "non-RBP")
  )
rbp.func.cat <- rbp.sel |> 
  add_column(
    group = NA
  ) |> 
  mutate(
    group = ifelse(rnase, "RNases", group),
    group = ifelse(rna.pol, "RNA polymerase", group),
    group = ifelse(ribosome, "Ribosomal proteins", group),
    group = ifelse(aatrna.synthesis, "tRNA ligases", group),
    group = ifelse(translation, "Translation factors", group),
    group = ifelse(rna.modif, "RNA modification", group),
    group = ifelse(ribosome.biogenesis, "Ribosome biogenesis", group),
    group = ifelse(glycolysis, "Glycolysis", group)
  ) |> 
  mutate(
    group = ifelse(rbp & is.na(group), "Other RBP", group)
  )
d <- rbind(rbp.rbp, rbp.func.cat[!is.na(rbp.func.cat$group),])
d$group <- factor(d$group, levels = c("non-RBP",
                                      "RBP",
                                      "RNA polymerase",
                                      "RNA modification",
                                      "RNases",
                                      "Ribosomal proteins",
                                      "Ribosome biogenesis",
                                      "tRNA ligases",
                                      "Translation factors",
                                      "Other RBP",
                                      "Glycolysis"))
d.enr <- d[d$p.adj.RBP.log < 0.05 & d$mean.ratio.RBP.log > 1,]
d.nenr <- d[d$p.adj.RBP.log >= 0.05 | d$mean.ratio.RBP.log <= 1,]

# Venn diagram for OOPS, RBS-ID and RBP annotation

setwd("./Figures")
venn <- venn.diagram(list(rbp = rbp$UniProt[rbp$rbp],
                          oops = rbp.sel$UniProt[rbp.sel$enriched.log],
                          rbsid = rbp$UniProt[rbp$CL.identified]),
                     category = c("RBPs", "OOPS", "RBS-ID"),
                     fontfamily = "sans",
                     cat.fontfamily = "sans",
                     filename = "Venn - RBP, OOPS, RBS-ID.tiff",
                     fill = c("#66c2a5", "#fc8d62", "#8da0cb"))
setwd("..")
rm(venn)

# Plot the number of proteins with CLs for RBP types

df_cl <- d |> 
  group_by(group) |> 
  summarise(
    prot = n(),
    cl_plus = sum(CL.identified)
  ) |> 
  mutate(cl_minus = prot - cl_plus) |> 
  pivot_longer(
    cols = c("cl_plus", "cl_minus"),
    names_to = "cl_status",
    values_to = "number"
  )

ggplot(df_cl, aes(x = group, y = number, fill = cl_status)) +
  geom_bar(stat = "identity", position="fill") +
  geom_text(aes(label = number), position = "fill", hjust = 3) +
  scale_x_discrete(limits = rev(levels(df_cl$group))) +
  coord_flip()

pdf(file = "./Figures/Number of CL in RBP types.pdf",
    width = 5, height = 2)
ggplot(df_cl, aes(x = group, y = number, fill = cl_status, color = group)) +
  geom_bar(stat = "identity",
           position="fill",
           size = 1,
           width = 0.75) +
  scale_fill_manual(values = c("white","grey")) +
  scale_color_manual(values = c("grey",'#66c2a5','#e78ac3','#e78ac3','#e78ac3',
                                '#A6761D','#A6761D','#A6761D','#A6761D','#66c2a5',"grey")) +
  geom_text(aes(label = number), position = "fill", color = 'black', vjust = 3) +
  theme_void() +
  theme(legend.position = 'none')
dev.off()

# Plot Beeswarm plots of OOPS enrichment for RBP types
d <- filter(d, !is.na(mean.ratio.RBP.log))
d.enr <- filter(d.enr, !is.na(mean.ratio.RBP.log))
d.nenr <- filter(d.nenr, !is.na(mean.ratio.RBP.log))

pdf(file = "./Figures/Beeswarm - OOPS enrichment of RBPs and CL - 2.pdf",
     width = 6, height = 4.5)
par(mar=c(10,5,1,1))
beeswarm(d.enr$mean.ratio.RBP.log ~ d.enr$group,
         col = c("grey",'#66c2a5','#e78ac3','#e78ac3','#e78ac3',
                 '#A6761D','#A6761D','#A6761D','#A6761D','#66c2a5',"grey"),
         corral = "wrap",
         pch = 16,  
         ylim = c(min(d$mean.ratio.RBP.log), max(d$mean.ratio.RBP.log)),
         xlab = "", ylab = "OOPS enrichment, log2", 
         las = 2, xaxt = "n")
beeswarm(d.nenr$mean.ratio.RBP.log ~ d.nenr$group,
         add = TRUE,
         col = "black", corral = "wrap",
         pch = 16)
abline(h = 0, lty = 3)
axis(1, at = 1:11, las = 2, c("non-RBP",
                              "RBP",
                              "RNA polymerase",
                              "RNA modification",
                              "RNases",
                              "Ribosomal proteins",
                              "Ribosome biogenesis",
                              "tRNA ligases",
                              "Translation factors",
                              "Other RBPs",
                              "Glycolysis"))

dev.off()
#rm(d, d.enr, d.nenr, rbp.sel, rbp.rbp, rbp.func.cat)

# Density plots of OOPS scores for proteins, RBPs and CL proteins

p <- ggplot(data = filter(rbp.sel, !rbp), aes(x = mean.ratio.RBP.log)) +
  geom_density() +
  geom_density(data = filter(rbp.sel, rbp), aes(x = mean.ratio.RBP.log), color = '#66c2a5') +
  geom_density(data = filter(rbp.sel, CL.identified), aes(x = mean.ratio.RBP.log), color = '#8da0cb') +
  ggtitle(label = "Density plots of OOPS scores") +
  scale_x_continuous(name = "OOPS enrichment, log2",
                     breaks = seq(-2, 7, by = 1),
                     minor_breaks = NULL) +
  scale_y_continuous(minor_breaks = NULL) +
  theme_bw()
pdf(file = "./Figures/Density plots of OOPS scores.pdf",
    width = 5, height = 5)
p
dev.off()

# Create volcano plot for Log showing RBPs

d <- subset(rbp, !is.na(mean.ratio.RBP.log),
            select = c(ENSG, rbp, Gene.name, mean.ratio.RBP.log, p.adj.RBP.log, CL.identified,
                       rbp, rnase, rna.pol, ribosome, aatrna.synthesis))
d.sel <- filter(d, rbp)

p <- ggplot(data = d)  +
  geom_point(mapping = aes(x = mean.ratio.RBP.log,
                           y = -log10(p.adj.RBP.log)),
             alpha = 0.1, size = 2) +
  geom_point(data = d.sel,
             mapping = aes(x = mean.ratio.RBP.log,
                           y = -log10(p.adj.RBP.log)),
             pch = 21, color = 'black', fill = '#66c2a5', size = 3) +
  ggtitle(label = "Candidate RBPs - Log") +
  scale_x_continuous(name = "OOPS enrichment, log2",
                     breaks = seq(-2, 7, by = 1),
                     minor_breaks = NULL) +
  scale_y_continuous(name = "FDR adjusted p value, -log10",
                     minor_breaks = NULL) +
  theme_bw()
pdf(file = "./Figures/Volcano - OOPS enrichment of RBPs.pdf",
    width = 5, height = 5)
p
dev.off()
rm(d, d.sel, p)

# Create volcano plots for Log showing different RBP categories

d <- subset(rbp, !is.na(mean.ratio.RBP.log),
            select = c(Gene.name, mean.ratio.RBP.log, p.adj.RBP.log, CL.identified,
                       rbp, rnase, rna.pol, ribosome, aatrna.synthesis))

plots <- vector("list", length = 4)
categories <- c("ribosome", "aatrna.synthesis", "rna.pol", "rnase")
titles <- c("Ribosomal proteins - Log", "tRNA ligases", "RNA polymerase", "Ribonucleases")

for(i in 1:length(categories)){
  d.sel <- d[d[,categories[i]],]
  p <- ggplot(data = d)  +
    geom_point(mapping = aes(x = mean.ratio.RBP.log,
                             y = -log10(p.adj.RBP.log)),
               alpha = 0.1, size = 2) +
    geom_point(data = d.sel,
               mapping = aes(x = mean.ratio.RBP.log,
                             y = -log10(p.adj.RBP.log)),
               color = "red", size = 3) +
    geom_text_repel(data = d.sel,
                    mapping = aes(x = mean.ratio.RBP.log,
                                  y = -log10(p.adj.RBP.log),
                                  label = Gene.name),
                    max.overlaps = Inf) +
    ggtitle(label = titles[i]) +
    scale_x_continuous(name = "OOPS enrichment, log2",
                       breaks = seq(-2, 7, by = 1)) +
    scale_y_continuous(name = "FDR adjusted p value, -log10") +
    theme_bw() +
    theme(panel.grid.minor=element_blank(),
          plot.title = element_text(hjust = 0.5))
  plots[[i]] <- p
}
pdf(file = "./Figures/Volcano - OOPS enrichment of RBP types.pdf",
    width = 12, height = 10)
ggarrange(plotlist = plots, nrow = 2, ncol = 2)
dev.off()
rm(d, d.sel, p, plots, categories, i, titles)

# Volcano plot for Log showing the candidate RBPs

d <- subset(rbp, !is.na(mean.ratio.RBP.log),
            select = c(ENSG, Gene.name, mean.ratio.RBP.log, p.adj.RBP.log, CL.identified,
                       rbp, rnase, rna.pol, ribosome, aatrna.synthesis))
d.sel <- filter(d, ENSG %in% c("SPy_0267", "SPy_1371", "SPy_0316", "SPy_0471",
                               "SPy_0539", "SPy_1124", "SPy_1608"))
d.sel[d.sel$ENSG == "SPy_0316", "Gene.name"] <- "yebC"
d.sel[d.sel$ENSG == "SPy_1124", "Gene.name"] <- "yjbK"
d.sel[d.sel$ENSG == "SPy_0539", "Gene.name"] <- "thuC"
d.sel[d.sel$ENSG == "SPy_1608", "Gene.name"] <- "ygaC"

p <- ggplot(data = d)  +
  geom_point(mapping = aes(x = mean.ratio.RBP.log,
                           y = -log10(p.adj.RBP.log)),
             alpha = 0.1, size = 2) +
  geom_point(data = d.sel,
             mapping = aes(x = mean.ratio.RBP.log,
                           y = -log10(p.adj.RBP.log)),
             color = "red", size = 3) +
  geom_label_repel(data = d.sel,
                  mapping = aes(x = mean.ratio.RBP.log,
                                y = -log10(p.adj.RBP.log),
                                label = Gene.name),
                  max.overlaps = Inf) +
  ggtitle(label = "Candidate RBPs - Log") +
  scale_x_continuous(name = "OOPS enrichment, log2",
                     breaks = seq(-2, 7, by = 1),
                     minor_breaks = NULL) +
  scale_y_continuous(name = "FDR adjusted p value, -log10",
                     minor_breaks = NULL) +
  theme_bw()
pdf(file = "./Figures/Volcano - OOPS enrichment for PNK assay.pdf",
    width = 5, height = 5)
p
dev.off()
rm(d, d.sel, p)

# Pie chart for functional categories of the candidate RBPs

pdf(file = "./Figures/Pie chart - func cat.pdf",
    width = 7, height = 7)
rbp_cat <- data.frame(group = c("Unknown",
                                "Signalling",
                                "Transport",
                                "Protein metabolism",
                                "DNA binding",
                                "Metabolic enzyme"),
                      value = c(10, 1, 2, 2, 9, 6))
ggplot(rbp_cat, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1, color = 'black') +
  geom_text(aes(label = value),
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c('#66c2a5','#8da0cb','#fc8d62','#e78ac3','#a6d854', 'grey')) +
  coord_polar("y", start=0) +
  theme_void()
dev.off()

# ROC curves for RBP identification (all proteins)

pdf(file = "./Figures/ROC curve - all proteins.pdf",
    width = 7, height = 7)
ggplot(rbp) +
  geom_roc(
    aes(d = as.numeric(rbp), m = mean.ratio.RBP.log),
    cutoffs.at = c(0.585, 1, 2), labels = FALSE, pointsize = 1,
    ) +
  geom_roc(
    aes(d = as.numeric(rbp & (Rep_1 + Rep_2 + Rep_3 >= 1)), m = mean.ratio.RBP.log),
    color = '#fc8d62', labels = FALSE, pointsize = 1,
    cutoffs.at = c(0.585, 1, 2)
    ) +
  geom_abline(slope = 1, linetype = "dashed") +
  coord_fixed(expand = FALSE) +
  scale_x_continuous(
    name = "FPR", breaks = seq(0, 1, 0.2), minor_breaks = FALSE
    ) +
  scale_y_continuous(
    name = "TPR", breaks = seq(0, 1, 0.2), minor_breaks = FALSE
    ) +
  theme_bw()
dev.off()

# ROC curves for RBP identification (only OOPS enriched proteins)

pdf(file = "./Figures/ROC curve - OOPS enriched proteins.pdf",
    width = 7, height = 7)
rbp_oops <- filter(rbp, enriched.log)
ggplot(rbp_oops) +
  geom_roc(
    aes(d = as.numeric(rbp), m = mean.ratio.RBP.log),
    cutoffs.at = c(0.585, 1, 2), labels = FALSE, pointsize = 1,
  ) +
  geom_roc(
    aes(d = as.numeric(rbp & (Rep_1 + Rep_2 + Rep_3 >= 1)), m = mean.ratio.RBP.log),
    color = '#fc8d62', labels = FALSE, pointsize = 1,
    cutoffs.at = c(0.585, 1, 2)
  ) +
  geom_abline(slope = 1, linetype = "dashed") +
  coord_fixed(expand = FALSE) +
  scale_x_continuous(
    name = "FPR", breaks = seq(0, 1, 0.2), minor_breaks = FALSE
  ) +
  scale_y_continuous(
    name = "TPR", breaks = seq(0, 1, 0.2), minor_breaks = FALSE
  ) +
  theme_bw()
dev.off()