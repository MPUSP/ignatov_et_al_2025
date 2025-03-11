library(tidyverse)
library(stringr)
library(ggplot2)
library(drawProteins)
library(ggrepel)


ifelse(!dir.exists("CL sites for PNK assay"),
       dir.create("CL sites for PNK assay"),
       "Folder 'CL sites for PNK assay' exists already")

ifelse(!dir.exists("CL sites for novel RBPs"),
       dir.create("CL sites for novel RBPs"),
       "Folder 'CL sites for novel RBPs' exists already")

rbp <- read.delim("Candidate RBPs.tsv")
ipr <- read.delim("Raw data/interpro.STRP1.tsv")
df_domains <- read_csv2("Raw data/Domains_PNK_assay.csv")

# Plot CL sites for proteins in the PNK assay
uniprot_ids  <-  c(
  "Q9A1H5", # Spy_0267
  "Q99Z67", # Spy_1371
  "P67188", # Spy_0316
  "Q9A145", # Spy_0471
  "Q9A0Z7", # Spy_0539
  "Q99ZQ9", # Spy_1124
  "Q99YP1"  # Spy_1608
)
uniprot_ids  <- rev(uniprot_ids)

df <- data.frame(
  accession = character(),
  type = character(),
  description = character(),
  entryName = character(),
  begin = integer(),
  end = integer(),
  length = integer(),
  taxid = integer(),
  order = numeric()
)

for(i in 1:length(uniprot_ids)) {
  accession = uniprot_ids[i]
  description = rbp[rbp$UniProt == accession, "ENSG"]
  entryName = description
  
  df_chain = data.frame(
    accession,
    type = "CHAIN",
    description = description,
    entryName = entryName,
    begin = 1,
    end = nchar(rbp[rbp$UniProt == accession, "Sequence"]),
    length = nchar(rbp[rbp$UniProt == accession, "Sequence"]),
    taxid = 301447,
    order = i
  )
  
  coords = as.numeric(str_match_all(rbp[rbp$UniProt == accession, "Modifications"], "(\\d+)\\w")[[1]][, 2])
  if (length(coords) > 0) {
    df_coord = data.frame(
      type = rep("CL site", length(coords)),
      description = rep(description, length(coords)),
      begin = coords,
      end = coords,
      length = 1,
      accession = rep(accession, length(coords)),
      entryName = rep(entryName, length(coords)),
      taxid = rep(301447, length(coords)),
      order = i
    )
    df = bind_rows(df, df_chain, df_coord)
  }
}

df = bind_rows(df, df_domains)

p <- draw_canvas(df)
p <- p + geom_segment(data = df[df$type == "CL site",],
                      aes(x=begin, y=order + 0.4, xend=begin, yend=order))
p <- p + geom_point(data = df[df$type == "CL site",], 
                    aes(x = begin, y = order + 0.4), shape = 23, 
                    colour = "black", fill = "yellow", size = 3, show.legend = FALSE)
p <- draw_chains(p, df)
p <- draw_domains(p, df)
p <- p + theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        panel.border = element_blank())
pdf(file = "./Figures/CL sites for PNK assay.pdf", width = 10, height = 6)
print(p)
dev.off()

# Functions to plot chains, domains and CL sites
draw_cl <- function(p, data = data, size = 5, fill = "yellow", show.legend = FALSE){
  begin = end = description = NULL
  p <- p + geom_segment(data = data[data$type == "CL site",],
                        aes(x=begin, y=order + 0.2, xend=begin, yend=order))
  p <- p + geom_point(data = data[data$type == "CL site",], 
                      aes(x = begin, y = order + 0.2), shape = 21, 
                      colour = "black", fill = fill, size = size, show.legend = show.legend)
  # p <- p + geom_label_repel(data = data[data$type == "CL site",],
  #                           aes(x = begin, y = order + 0.2, label = begin),
  #                           max.overlaps = 100,
  #                           nudge_y = 0.1)
  return(p)
}

draw_chains <- function (p, data = data, outline = "black", fill = "grey", 
                         label_chains = TRUE, 
                         labels = data[data$type == "CHAIN",]$entryName,
                         size = 0.5, label_size = 4){
  begin = end = NULL
  p <- p + geom_rect(data = data[data$type == "CHAIN",], 
                     mapping = aes(xmin = begin, xmax = end,
                                   ymin = order - 0.1, ymax = order + 0.1),
                     colour = outline, fill = fill, size = size)
  if(label_chains == TRUE){
    p <- p + annotate("text", x = -10, y = data[data$type == "CHAIN", ]$order,
                      label = labels, hjust = 1, size = label_size)
  }
  return(p)
}

draw_domains <- function (p, data = data, label_domains = TRUE, label_size = 4, 
                          show.legend = TRUE, type = "DOMAIN"){
  begin = end = description = NULL
  p <- p + geom_rect(data = data[data$type == type,],
                     mapping = aes(xmin = begin, xmax = end,
                                   ymin = order - 0.1, ymax = order + 0.1, fill = description), 
                     color = "black", show.legend = show.legend)
  if (label_domains == TRUE) {
    p <- p + geom_label(data = data[data$type == type, ],
                        aes(x = begin + (end - begin)/2, y = order, label = description),
                        size = label_size)
  }
  return(p)
}

# Function to drawing the protein with domains and CL sites
draw_prot <- function(up = "Q9A102", folder = "./CL sites for novel RBPs",
                      db.select = c("cd", "G3DSA", "MF_", "PF", "PIRSF", "PR",
                                    "PS", "PTHR", "SM", "SSF", "TIGR")){
  if(!(up %in% rbp[rbp$CL.identified, "UniProt"])){
    print("The protein does not have detected cross-links")
    return(NULL)
  }
  gene.name <- ifelse(rbp[rbp$UniProt == up, "Gene.name"] == "",
                      rbp[rbp$UniProt == up, "ENSG"],
                      rbp[rbp$UniProt == up, "Gene.name"])
  df.chain <- data.frame(type = "CHAIN",
                         description = gene.name,
                         begin = 1,
                         end = nchar(rbp[rbp$UniProt == up, "Sequence"]),
                         length = nchar(rbp[rbp$UniProt == up, "Sequence"]),
                         accession = up,
                         entryName = gene.name,
                         taxid = 301447,
                         order = 1)
  coords <- as.numeric(str_match_all(rbp[rbp$UniProt == up, "Modifications"],
                                     "(\\d+)\\w")[[1]][,2])
  df.coord <- data.frame(type = rep("CL site", length(coords)),
                         description = rep(gene.name, length(coords)),
                         begin = coords,
                         end = coords,
                         length = 1,
                         accession = rep(up, length(coords)),
                         entryName = rep(gene.name, length(coords)),
                         taxid = rep(301447, length(coords)),
                         order = 1)
  df.chain.coord <- rbind(df.chain, df.coord)
  ipr.up <- ipr[ipr$Uniprot == up,]
  if(nrow(ipr.up) == 0){
    df <- df.chain.coord
  }else{
    db.ids <- c("CDD", "CATH", "HAMAP", "Pfam", "PIR",
                "PRINTS", "Prosite", "Panther", "SMART", "Superfamily", "TIGR")
    names(db.ids) <- c("cd", "G3DSA", "MF_", "PF", "PIRSF",
                       "PR", "PS", "PTHR", "SM", "SSF", "TIGR")
    db.ids <-  db.ids[db.select]
    df <- data.frame(type = character(),
                     description = character(),
                     begin = numeric(),
                     end = numeric(),
                     length = numeric(),
                     accession = character(),
                     entryName = character(),
                     taxid = numeric(),
                     order = numeric())
    df.domains <- cbind(df)

    for(i in 1:nrow(ipr.up)){
      domain.name <- ipr.up[i, "Domain"]
      domain.id <- ipr.up[i, "Database.ID"]
      begin <- ipr.up[i, "Start"]
      end <- ipr.up[i, "End"]
      db.id <- db.ids[str_detect(domain.id, names(db.ids))]
      if(any(str_detect(domain.id, names(db.ids)))){
        df.line <- data.frame(type = db.id,
                              description = domain.name,
                              begin = begin,
                              end = end,
                              length = end - begin + 1,
                              accession = up,
                              entryName = gene.name,
                              taxid = 301447,
                              order = 1)
        df.domains <- rbind(df.domains, df.line)
      }
    }
    for(i in 1:length(unique(df.domains$type))){
      db.type = unique(df.domains$type)[i]
      df.db <- rbind(df.domains[df.domains$type == db.type,])
      df.db$type = "DOMAIN"
      df.db <- rbind(df.chain.coord, df.db)
      df.db$entryName <- db.type
      df.db$order <- i
      df <- rbind(df, df.db)
    }
  }
  p <- draw_canvas(df)
  p <- draw_cl(p, df)
  p <- draw_chains(p, df)
  p <- draw_domains(p, df, label_domains = FALSE)
  title <- paste(gene.name, rbp[rbp$UniProt == up, "Modifications"], sep = " - ")
  p <- p + labs(title = gene.name,
                subtitle = rbp[rbp$UniProt == up, "Modifications"])
  p <- p + theme_bw() + theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank())
  file_path <- paste(folder, "/", gene.name, ".pdf", sep = "")
  pdf(file = file_path, width = 10, height = 6)
  print(p)
  dev.off()
}

# Proteins for PNK assay
draw_prot("Q9A1H5", db.select = c("PF"), folder = "./CL sites for PNK assay")   # Spy_0267
draw_prot("Q99Z67", db.select = c("PF"), folder = "./CL sites for PNK assay")   # Spy_1371
draw_prot("P67188", db.select = c("PF"), folder = "./CL sites for PNK assay")   # Spy_0316
draw_prot("Q9A145", db.select = c("PF"), folder = "./CL sites for PNK assay")   # Spy_0471
draw_prot("Q9A0Z7", db.select = c("PF"), folder = "./CL sites for PNK assay")   # Spy_0539
draw_prot("Q99ZQ9", db.select = c("PF"), folder = "./CL sites for PNK assay")    # Spy_1124
draw_prot("Q99YP1", db.select = c("PF"), folder = "./CL sites for PNK assay")    # Spy_1608

# Novel RBPs
for(up in rbp[rbp$novel.rbp, "UniProt"]){
  print(up)
  draw_prot(up)
}
