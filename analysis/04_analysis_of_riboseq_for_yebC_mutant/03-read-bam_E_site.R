library(tidyverse)
library(GenomicAlignments)

# Read the GFF file
gff <- read_table(
  "../../data/genome/220215_NC_002737.2_refseq_gene.gff",
  skip = 2,
  col_names = c(
    "seqid",
    "source",
    "type",
    "start",
    "stop",
    "score",
    "strand",
    "phase",
    "attr"
  )
)
gff <- as.data.frame(gff)

# Function to extract attributes from the "attr" column of gff data frame
extract_attribute <- function(attr_column, attr_name) {
  attr_fields <- unlist(strsplit(attr_column, ";"))
  attr_name_re <- paste("^", attr_name, sep = "")
  extr_attr <- grep(attr_name_re, attr_fields, value = TRUE)
  substr(extr_attr, nchar(attr_name) + 2, nchar(extr_attr))
}

# Extract attributes and add them as columns
gff <- within(gff, {
  name <- extract_attribute(attr, "Name")
  biotype <- extract_attribute(attr, "gene_biotype")
  id <- extract_attribute(attr, "ID")
  locus_tag <- extract_attribute(attr, "locus_tag")
})

gff <- gff |> filter(
  biotype == "protein_coding",
  start < stop,
  !(abs(start - stop) + 1) %% 3
)

# Add a column with old locus tags
gff <- gff |> mutate(
  old_locus_tag = ifelse(str_detect(attr, "SPy_"),
                         str_extract(attr, "SPy_...."),
                         NA)
)

# Extract nt and aa sequences and add them as columns
genome_seq <- readDNAStringSet("../../data/genome/NC_002737.2.fa")[[1]]
gff$nt = rep(NA, nrow(gff))
gff$aa = rep(NA, nrow(gff))

for(i in 1:nrow(gff)) {
  nt_seq = subseq(genome_seq, gff[i, "start"], gff[i, "stop"])
  if(gff[i, "strand"] == "-") {nt_seq = reverseComplement(nt_seq)}
  aa_seq = translate(nt_seq)
  gff[i, "nt"] = as.character(nt_seq)
  gff[i, "aa"] = as.character(aa_seq)
}

# Create a data.frame to store the Ribo-seq data
rsd <- data.frame(
  locus_tag = character(),
  name = character(),
  old_locus_tag = character(),
  locus_type = character(),
  locus_start = integer(),
  locus_stop = integer(),
  locus_strand = character(),
  nt = character(),
  nt_locus_coord = integer(),
  nt_genome_coord = integer(),
  aa = character(),
  aa_locus_coord = integer()
)

# Add CDS positions
for(i in 1:nrow(gff)){
  nt = str_split_1(gff[i, "nt"], "")
  cds_length = length(nt)
  locus_tag = rep(gff[i, "locus_tag"], cds_length)
  name = rep(gff[i, "name"], cds_length)
  old_locus_tag = rep(gff[i, "old_locus_tag"], cds_length)
  locus_type = rep("CDS", cds_length)
  locus_start = rep(gff[i, "start"], cds_length)
  locus_stop = rep(gff[i, "stop"], cds_length)
  locus_strand = rep(gff[i, "strand"], cds_length)
  nt_locus_coord = seq(cds_length)
  aa_single = str_split_1(gff[i, "aa"], "")
  aa = rep(aa_single, each = 3)
  aa_locus_coord = rep(seq(length(aa_single)), each = 3)
  if(gff[i, "strand"] == '+'){
    nt_genome_coord = gff[i, "start"] + nt_locus_coord - 1
  }else{
    nt = rev(nt)
    nt_locus_coord = rev(nt_locus_coord)
    aa = rev(aa)
    aa_locus_coord = rev(aa_locus_coord)
    nt_genome_coord = gff[i, "stop"] - nt_locus_coord + 1
  }
  rsd_cds = data.frame(
    locus_tag,
    name,
    old_locus_tag,
    locus_type,
    locus_start,
    locus_stop,
    locus_strand,
    nt,
    nt_locus_coord,
    nt_genome_coord,
    aa,
    aa_locus_coord
  )
  rsd = rbind(rsd, rsd_cds)
}


# Add 30 nt 5 and 3UTRs 
rsd_utr <- data.frame(
  locus_tag = character(),
  name = character(),
  old_locus_tag = character(),
  locus_type = character(),
  locus_start = integer(),
  locus_stop = integer(),
  locus_strand = character(),
  nt = character(),
  nt_locus_coord = integer(),
  nt_genome_coord = integer(),
  aa = character(),
  aa_locus_coord = integer()
)

for(i in 1:nrow(gff)){
  utr_length = 30
  nt = rep("N", utr_length)
  locus_tag = rep(gff[i, "locus_tag"], utr_length)
  name = rep(gff[i, "name"], utr_length)
  old_locus_tag = rep(gff[i, "old_locus_tag"], utr_length)
  locus_type = rep("5UTR", utr_length)
  locus_strand = rep(gff[i, "strand"], utr_length)
  nt_locus_coord = seq(- utr_length + 1, 0)
  aa = rep(NA, utr_length)
  aa_locus_coord = rep(NA, utr_length)
  if(gff[i, "strand"] == '+'){
    locus_start = rep(gff[i, "start"] - utr_length, utr_length)
    locus_stop = rep(gff[i, "start"] - 1, utr_length)
    nt_genome_coord = gff[i, "start"] + nt_locus_coord - 1
  }else{
    locus_start = rep(gff[i, "stop"] + 1, utr_length)
    locus_stop = rep(gff[i, "stop"] + utr_length, utr_length)
    nt_locus_coord = rev(nt_locus_coord)
    nt_genome_coord = gff[i, "stop"] - nt_locus_coord + 1
  }
  rsd_5utr = data.frame(
    locus_tag,
    name,
    old_locus_tag,
    locus_type,
    locus_start,
    locus_stop,
    locus_strand,
    nt,
    nt_locus_coord,
    nt_genome_coord,
    aa,
    aa_locus_coord
  )
  rsd_utr = rbind(rsd_utr, rsd_5utr)
}

# Add 30 nt 3UTRs 
for(i in 1:nrow(gff)){
  cds_length = gff[i, "stop"] - gff[i, "start"] + 1
  utr_length = 30
  nt = rep("N", utr_length)
  locus_tag = rep(gff[i, "locus_tag"], utr_length)
  name = rep(gff[i, "name"], utr_length)
  old_locus_tag = rep(gff[i, "old_locus_tag"], utr_length)
  locus_type = rep("3UTR", utr_length)
  locus_strand = rep(gff[i, "strand"], utr_length)
  nt_locus_coord = seq(cds_length + 1, cds_length + utr_length)
  aa = rep(NA, utr_length)
  aa_locus_coord = rep(NA, utr_length)
  if(gff[i, "strand"] == '+'){
    locus_start = rep(gff[i, "stop"] + 1, utr_length)
    locus_stop = rep(gff[i, "stop"] + utr_length, utr_length)
    nt_genome_coord = gff[i, "start"] + nt_locus_coord - 1
  }else{
    locus_start = rep(gff[i, "start"] - utr_length, utr_length)
    locus_stop = rep(gff[i, "start"] - 1, utr_length)
    nt_locus_coord = rev(nt_locus_coord)
    nt_genome_coord = gff[i, "stop"] - nt_locus_coord + 1
  }
  rsd_3utr = data.frame(
    locus_tag,
    name,
    old_locus_tag,
    locus_type,
    locus_start,
    locus_stop,
    locus_strand,
    nt,
    nt_locus_coord,
    nt_genome_coord,
    aa,
    aa_locus_coord
  )
  rsd_utr = rbind(rsd_utr, rsd_3utr)
}

rsd <- rbind(rsd, rsd_utr)

# Add reverse numeration for each locus
rsd <- rsd |>
  group_by(locus_tag) |>
  mutate(nt_locus_coord_rev = n() - nt_locus_coord - 59)

# Function to read BAM files and add to rsd data.frame
read_bam <- function(sample_names = c("wt_1"), read_min = 27, read_max = 40){
  
  for(i in 1:length(sample_names)){
    sample_name <-  sample_names[i]
    
    # Read BAM file as GAlignments object and extract reads of certain lengths
    bam_file <- paste("../../data/riboseq/riboseq-", sample_name, ".bam", sep = "")
    bam <- readGAlignments(bam_file)
    bam <- bam[between(width(bam), read_min, read_max)]
    
    # Calculate coverage and add to the cvg_df
    bam_plus <- bam[strand(bam) == "+"]
    bam_minus <- bam[strand(bam) == "-"]
    pos_plus <- end(bam_plus) - 18
    pos_minus <- start(bam_minus) + 18
    cvg_plus <- tabulate(pos_plus, nbins = length(genome_seq))
    cvg_minus <- tabulate(pos_minus, nbins = length(genome_seq))
    cvg_df <- data.frame(
      nt_genome_coord = seq(1, length(genome_seq)),
      cvg_plus,
      cvg_minus
    )
    
    # Merge with rsd data.frame
    rsd <- rsd |> left_join(cvg_df)
    rsd <- rsd |> mutate(
      cvg = ifelse(locus_strand == "+",
                   cvg_plus, cvg_minus))
    rsd <- rsd |> rename(
      cvg_plus = paste("cvg_", sample_name, "_plus", sep = ""),
      cvg_minus = paste("cvg_", sample_name, "_minus", sep = ""),
      cvg = paste("cvg_", sample_name, sep = "")
    )
    
  }
  rsd
}

# Run the function and save results to CSV file
rsd <- read_bam(c("wt_1", "d0316_1", "comp_1", "m6_1",
                  "wt_2", "d0316_2", "comp_2", "m6_2",
                  "wt_3", "d0316_3", "comp_3", "m6_3"))
rsd <- rsd |> arrange(nt_genome_coord)
write_csv(rsd, file = "riboseq_E_site.csv")


##################### End #####################
