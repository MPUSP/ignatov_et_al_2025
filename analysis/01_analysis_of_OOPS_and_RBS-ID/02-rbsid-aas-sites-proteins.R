library(tidyverse)
library(UniProt.ws)
library(stringr)
library(VennDiagram)

rbsid <- read.delim("Raw data/RBS-ID results.tsv")
rbsid$Replicate <- as.factor(rbsid$Replicate)
ipr <- read.delim("Raw data/interpro.STRP1.tsv")

# Detect cross-links and select peptides with exactly one cross-link
proper_cl <- rep(NA, nrow(rbsid))
for(i in 1:nrow(rbsid)){
  mods <- str_match_all(rbsid[i,"Assigned.Modifications"],
                        "(\\d+)(\\w)\\(((?:244|226|227).\\d+)\\)")[[1]]
  proper_cl[i] <- ifelse(nrow(mods) == 1, TRUE, FALSE)
}
rbsid <- rbsid[proper_cl,]
rm(proper_cl)

# Remove peptides with spectral count equal to zero
rbsid <- rbsid[rbsid$Spectral.Count != 0,]

# Identify cross-linked proteins in each replicate
clp <- data.frame(UniProt = unique(rbsid$Protein.ID))
clp <- within(clp, {
  Rep_1 <- is.element(UniProt, rbsid$Protein.ID[rbsid$Replicate == "Rep_1"])
  Rep_2 <- is.element(UniProt, rbsid$Protein.ID[rbsid$Replicate == "Rep_2"])
  Rep_3 <- is.element(UniProt, rbsid$Protein.ID[rbsid$Replicate == "Rep_3"])
})

# Extract the protein sequences from Uniprot and remove entries without sequence
up <- UniProt.ws(taxId=301447)
up_output <- select(up, keys = clp$UniProt, columns = c("sequence"))
clp$Sequence <- up_output$Sequence
clp <- clp[!is.na(clp$Sequence),]
rm(up, up_output)

# Remove peptides from proteins without sequence in Uniprot
rbsid <- rbsid[is.element(rbsid$Protein.ID, clp$UniProt),]

# Identify the positions and amino acids of CL within proteins
rbsid$Coord.prot <- rep(NA, nrow(rbsid))
rbsid$AA <- rep(NA, nrow(rbsid))
rbsid$Mass.shift <- rep(NA, nrow(rbsid))

for(i in 1:nrow(rbsid)){
  uniprot <- rbsid[i, "Protein.ID"]
  pep_seq <- rbsid[i, "Peptide.Sequence"]
  mods <- str_match(rbsid[i,"Assigned.Modifications"],
                    "(\\d+)(\\w)\\(((?:244|226|227).\\d+)\\)")
  coord <- as.integer(mods[2])
  aa <- mods[3]
  
  prot_seq <- clp[clp$UniProt == uniprot, "Sequence"]
  prot_length <- nchar(prot_seq)
  pep_coords <- str_locate(prot_seq, pep_seq)
  
  rbsid$Coord.prot[i] <- pep_coords[1] + coord - 1
  rbsid$AA[i] <- mods[3]
  rbsid$Mass.shift[i] <- mods[4]
}
rm(i, uniprot, pep_seq, mods, coord, aa, prot_seq, prot_length, pep_coords)

# For each replicate detect CL sites and sum up their spectral counts
cls <- data.frame(Protein.ID = character(),
                  Coord.prot = integer(),
                  Rep_1 = integer(),
                  Rep_2 = integer(),
                  Rep_3 = integer())
for(up in clp$UniProt){
  rbsid.up <- rbsid[rbsid$Protein.ID == up,]
  coord.up <- sort(unique(rbsid.up$Coord.prot))
  cls.up <- data.frame(Protein.ID = up,
                       Coord.prot = coord.up,
                       Rep_1 = 0, Rep_2 = 0, Rep_3 = 0)
  for(coord in coord.up){
    rbsid.up.coord <- rbsid.up[rbsid.up$Coord.prot == coord,]
    sc.coord <- tapply(rbsid.up.coord$Spectral.Count, rbsid.up.coord$Replicate, sum)
    cls.up[cls.up$Coord.prot == coord, c("Rep_1", "Rep_2", "Rep_3")] <- sc.coord
  }
  cls <- rbind(cls, cls.up)
}
cls$CL.site <- sprintf("%s - %d", cls$Protein.ID, cls$Coord.prot)
rm(up, rbsid.up, coord.up, cls.up, coord, rbsid.up.coord, sc.coord)

# Identify the CL amino acid for each CL site
cls$AA <- rep(NA, nrow(cls))
for(i in 1:nrow(cls)){
  cls$AA[i] <- rbsid$AA[rbsid$Protein.ID == cls[i, "Protein.ID"]
                        & rbsid$Coord.prot == cls[i, "Coord.prot"]][1]
}
rm(i)

# Add the CL amino acids and their positions to CL proteins
clp$Modifications <- rep(NA, nrow(clp))
for(i in 1:nrow(clp)){
  up <- clp[i, "UniProt"]
  cls.sel <- cls[cls$Protein.ID == up,]
  mods <- sprintf("%d%s", cls.sel$Coord.prot, cls.sel$AA)
  clp[i, "Modifications"] <- paste(mods, collapse = '; ')
}
rm(i, up, cls.sel, mods)

# Add the domain information for CL sites
db.ids <- c("CDD", "CATH", "HAMAP", "Pfam", "PIR",
            "PRINTS", "Prosite", "Panther", "SMART", "Superfamily", "TIGR")
names(db.ids) <- c("cd", "G3DSA", "MF_", "PF", "PIRSF",
                   "PR", "PS", "PTHR", "SM", "SSF", "TIGR")
for(db.id in db.ids){
  cls[db.id] <- rep(NA, nrow(cls))
}
for(i in 1:nrow(cls)){
  up <- cls[i, "Protein.ID"]
  coord <- cls[i, "Coord.prot"]
  ipr.up <- ipr[ipr$Uniprot == up,]
  if(nrow(ipr.up) > 0){
    for(j in 1:nrow(ipr.up)){
      domain.name <- ipr.up[j, "Domain"]
      domain.id <- ipr.up[j, "Database.ID"]
      start <- ipr.up[j, "Start"]
      end <- ipr.up[j, "End"]
      if(coord >= start & coord <= end){
        for(db.id.name in names(db.ids)){
          if(str_detect(domain.id, db.id.name)){
            cls[i, db.ids[db.id.name]] <- sprintf("%s (%d..%d)",
                                                  domain.name, start, end)
          }
        }
      }
    }
  }
}
rm(db.ids, db.id, i, up, coord, ipr.up, j, domain.name, domain.id, start, end, 
   db.id.name)

# Draw Venn diagrams for CL sites and proteins in different RBS-ID replicates
ifelse(!dir.exists("Figures"), dir.create("Figures"), "Folder 'Figures' exists already")
setwd("./Figures")
venn.diagram(list(rep_1 = cls[!is.na(cls$Rep_1), "CL.site"],
                  rep_2 = cls[!is.na(cls$Rep_2), "CL.site"],
                  rep_3 = cls[!is.na(cls$Rep_3), "CL.site"]),
             filename = "Venn - CL sites in RBS-ID replicates.tiff",
             category = c("Replicate 1", "Replicate 2", "Replicate 3"),
             fontfamily = "sans",
             cat.fontfamily = "sans",
             fill = c("#c2cce3", "#8da0cb", "#6A83bb"))
venn.diagram(list(rep_1 = clp[clp$Rep_1, "UniProt"],
                  rep_2 = clp[clp$Rep_2, "UniProt"],
                  rep_3 = clp[clp$Rep_3, "UniProt"]),
             filename = "Venn - CL proteins in RBS-ID replicates.tiff",
             category = c("Replicate 1", "Replicate 2", "Replicate 3"),
             fontfamily = "sans",
             cat.fontfamily = "sans",
             fill = c("#c2cce3", "#8da0cb", "#6A83bb"))
setwd("..")

# Calculate the CL frequency for different amino acids
aa.frac.prot <- function(rep.name = ""){
  if(rep.name == ""){
    concat.seq <- paste(clp[,"Sequence"], collapse = "")
  }else{
    concat.seq <- paste(clp[clp[,rep.name], "Sequence"], collapse = "")
  }
  aa.count <- table(strsplit(concat.seq, ""))
  aa.count <- as.vector(aa.count[aa_stat$aa])
  aa.frac <- aa.count / sum(aa.count)
}
aa_stat <- data.frame(aa = c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S",
                             "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T"))
aa_stat <- within(aa_stat, {
  cl = table(cls$AA)[aa]
  cl_1 = table(cls$AA[!is.na(cls$Rep_1)])[aa]
  cl_2 = table(cls$AA[!is.na(cls$Rep_2)])[aa]
  cl_3 = table(cls$AA[!is.na(cls$Rep_3)])[aa]
  cl.frac <- cl / sum(cl)
  cl_1.frac <- cl_1 / sum(cl_1)
  cl_2.frac <- cl_2 / sum(cl_2)
  cl_3.frac <- cl_3 / sum(cl_3)
  cl.freq <- cl.frac / aa.frac.prot()
  cl_1.freq <- cl_1.frac / aa.frac.prot("Rep_1")
  cl_2.freq <- cl_2.frac / aa.frac.prot("Rep_2")
  cl_3.freq <- cl_3.frac / aa.frac.prot("Rep_3")
})

# Plot the CL amino acids counts and frequencies in total and for replicates
plot.cl.aa <- function(name.rep = "cl"){
  name.freq <- paste(name.rep, ".freq", sep = "")
  title <- switch(name.rep,
                  cl = "Crosslinked amino acids",
                  cl_1 = "Crosslinked amino acids in replicate 1",
                  cl_2 = "Crosslinked amino acids in replicate 2",
                  cl_3 = "Crosslinked amino acids in replicate 3",)
  file.path <- paste("Figures/", title, ".pdf", sep = "")
  aa_stat <- aa_stat[order(aa_stat[, name.freq], decreasing = TRUE),]
  pdf(file = file.path, width = 9, height = 5.5)
  par(mfrow=c(2,1),
      oma = c(5,4,1,1) + 1,
      mar = c(1,4,2,2) + 0.1)
  barplot(aa_stat[, name.freq],
          names.arg = aa_stat$aa,
          main = title,
          ylab = "Frequency",
          ylim = c(0, max(aa_stat[, name.freq]) + 2),
          las = 1,
          col = "blue4")
  barplot(aa_stat[, name.rep],
          ylab = "RBS counts",
          ylim = c(max(aa_stat[, name.rep]) + 50, 0),
          las = 1,
          col = "#0099ff")
  par(mfrow=c(1,1))
  dev.off()
}
for(name.rep in c("cl", "cl_1", "cl_2", "cl_3")) plot.cl.aa(name.rep)
rm(aa.frac.prot, plot.cl.aa, name.rep)

# Statistics of RNA-binding sites and domains
func_cat <- cls |>
  transform(Func_cat = str_replace(Pfam, " \\(.+", "")) |> 
  transform(Func_cat = str_replace(Func_cat, "Ribosomal protein.*", "Ribosomal protein")) |> 
  count(Func_cat) |> 
  arrange(desc(n)) |> 
  slice_head(n = 10)

func_cat

# Write the tables of CL sites and proteins to files
write.table(cls, "RBS-ID - CL sites.tsv", row.names = FALSE,
            quote = FALSE, sep = "\t")
write.table(clp, "RBS-ID - CL proteins.tsv", row.names = FALSE,
            quote = FALSE, sep = "\t")