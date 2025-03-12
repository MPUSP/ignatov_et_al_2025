library(tidyverse)
library(GenomicAlignments)

ifelse(!dir.exists("igv"), dir.create("igv"), "Folder 'igv' exists already")
ifelse(!dir.exists("figures"), dir.create("figures"), "Folder 'figures' exists already")

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

# Extract nt sequences and add them as columns
genome_seq <- readDNAStringSet("../../data/genome/NC_002737.2.fa")[[1]]
gff$nt = rep(NA, nrow(gff))

for(i in 1:nrow(gff)) {
  nt_seq = subseq(genome_seq, gff[i, "start"], gff[i, "stop"])
  if(gff[i, "strand"] == "-") {nt_seq = reverseComplement(nt_seq)}
  gff[i, "nt"] = as.character(nt_seq)
}

# Create a data.frame to store data
iclip_df <- data.frame(
  nt_genome_coord = seq(1, length(genome_seq)),
  nt_locus_coord = rep(NA, length(genome_seq)),
  locus_tag = rep(NA, length(genome_seq)),
  biotype = rep('intergenic', length(genome_seq))
)

for(i in 1:nrow(gff)){
  locus_tag = gff[i, 'locus_tag']
  biotype = gff[i, 'biotype']
  start = gff[i, 'start']
  stop = gff[i, 'stop']
  strand = gff[i, 'strand']
  gene_length = stop - start + 1
  
  iclip_df[start:stop, "locus_tag"] = locus_tag
  iclip_df[start:stop, "biotype"] = biotype
  
  if(strand == '+'){
    iclip_df[start:stop, "nt_locus_coord"] = 1:gene_length
  } else {
    iclip_df[start:stop, "nt_locus_coord"] = rev(1:gene_length)
  }
}

# Read BAM files and collect the data to iclip_df
for(sample_name in c(
  'yebc_uv_1',
  'yebc_uv_2',
  'smi_1',
  'smi_2',
  'yebc_no_uv',
  'wt_uv'
)){
  #sample_name <- 'yebc_uv_1'
  bam_file <- paste("../../data/iclip/iCLIP-0316-1124_", sample_name, ".bam", sep = "")
  bam <- readGAlignments(bam_file)
  
  bam_plus <- bam[strand(bam) == "+"]
  bam_minus <- bam[strand(bam) == "-"]
  
  # Extract CL sites and calculate coverage
  cl_sites_plus <- start(bam_plus) - 1
  cl_sites_minus <- end(bam_minus) + 1
  cl_plus <- tabulate(cl_sites_plus, nbins = length(genome_seq))
  cl_minus <- tabulate(cl_sites_minus, nbins = length(genome_seq))
  
  cvg_plus <- coverage(bam_plus)
  cvg_minus <- coverage(bam_minus)
  cvg_plus <- as.numeric(cvg_plus$NC_002737.2)
  cvg_minus <- as.numeric(cvg_minus$NC_002737.2)
  
  # Collapse CL sites of rRNAs to the first operon
  df_16S <- data.frame(
    cl_SPY_RS00070 = cl_plus[17065:18613],
    cl_SPY_RS00135 = cl_plus[23066:24614],
    cl_SPY_RS00480 = cl_plus[79284:80832],
    cl_SPY_RS01390 = cl_plus[264349:265898][-419],
    cl_SPY_RS06725 = rev(cl_minus[1334635:1336183]),
    cl_SPY_RS07920 = rev(cl_minus[1581569:1583117]),
    
    cvg_SPY_RS00070 = cvg_plus[17065:18613],
    cvg_SPY_RS00135 = cvg_plus[23066:24614],
    cvg_SPY_RS00480 = cvg_plus[79284:80832],
    cvg_SPY_RS01390 = cvg_plus[264349:265898][-419],
    cvg_SPY_RS06725 = rev(cvg_minus[1334635:1336183]),
    cvg_SPY_RS07920 = rev(cvg_minus[1581569:1583117])
  ) |>
    rowwise() |> 
    mutate(
      cl_sum = sum(c_across(starts_with("cl"))),
      cvg_sum = sum(c_across(starts_with("cvg")))
    )
  
  cl_plus[17065:18613] <- df_16S$cl_sum
  cl_plus[23066:24614] <- 0
  cl_plus[79284:80832] <- 0
  cl_plus[264349:265898] <- 0
  cl_minus[1334635:1336183] <- 0
  cl_minus[1581569:1583117] <- 0
  
  cvg_plus[17065:18613] <- df_16S$cvg_sum
  cvg_plus[23066:24614] <- 0
  cvg_plus[79284:80832] <- 0
  cvg_plus[264349:265898] <- 0
  cvg_minus[1334635:1336183] <- 0
  cvg_minus[1581569:1583117] <- 0
  
  df_23S <- data.frame(
    cl_SPY_RS00080 = cl_plus[19036:21938],
    cl_SPY_RS00145 = cl_plus[25037:27939],
    cl_SPY_RS00490 = cl_plus[81255:84158][-2399],
    cl_SPY_RS01400 = cl_plus[266321:269226][-c(2104, 2401, 2447)],
    cl_SPY_RS06715 = rev(cl_minus[1331308:1334212])[-c(2104, 2401)],
    cl_SPY_RS07910 = rev(cl_minus[1578235:1581146])[-c(2104, 2538, 2548, 2578, 2581, 2602, 2618, 2624, 2629)],
    
    cvg_SPY_RS00080 = cvg_plus[19036:21938],
    cvg_SPY_RS00145 = cvg_plus[25037:27939],
    cvg_SPY_RS00490 = cvg_plus[81255:84158][-2399],
    cvg_SPY_RS01400 = cvg_plus[266321:269226][-c(2104, 2401, 2447)],
    cvg_SPY_RS06715 = rev(cvg_minus[1331308:1334212])[-c(2104, 2401)],
    cvg_SPY_RS07910 = rev(cvg_minus[1578235:1581146])[-c(2104, 2538, 2548, 2578, 2581, 2602, 2618, 2624, 2629)]
  ) |>
    rowwise() |> 
    mutate(
      cl_sum = sum(c_across(starts_with("cl"))),
      cvg_sum = sum(c_across(starts_with("cvg")))
    )
  
  cl_plus[19036:21938] <- df_23S$cl_sum
  cl_plus[25037:27939] <- 0
  cl_plus[81255:84158] <- 0
  cl_plus[266321:269226] <- 0
  cl_minus[1331308:1334212] <- 0
  cl_minus[1578235:1581146] <- 0
  
  cvg_plus[19036:21938] <- df_23S$cvg_sum
  cvg_plus[25037:27939] <- 0
  cvg_plus[81255:84158] <- 0
  cvg_plus[266321:269226] <- 0
  cvg_minus[1331308:1334212] <- 0
  cvg_minus[1578235:1581146] <- 0
  
  df_5S <- data.frame(
    cl_SPY_RS00085 = cl_plus[22024:22139],
    cl_SPY_RS00150 = cl_plus[28025:28140],
    cl_SPY_RS00495 = cl_plus[84244:84359],
    cl_SPY_RS01405 = cl_plus[269312:269427],
    cl_SPY_RS06710 = rev(cl_minus[1331107:1331222]),
    cl_SPY_RS07905 = rev(cl_minus[1578034:1578149]),
    
    cvg_SPY_RS00085 = cvg_plus[22024:22139],
    cvg_SPY_RS00150 = cvg_plus[28025:28140],
    cvg_SPY_RS00495 = cvg_plus[84244:84359],
    cvg_SPY_RS01405 = cvg_plus[269312:269427],
    cvg_SPY_RS06710 = rev(cvg_minus[1331107:1331222]),
    cvg_SPY_RS07905 = rev(cvg_minus[1578034:1578149])
  ) |>
    rowwise() |> 
    mutate(
      cl_sum = sum(c_across(starts_with("cl"))),
      cvg_sum = sum(c_across(starts_with("cvg")))
    )
  
  cl_plus[22024:22139] <- df_5S$cl_sum
  cl_plus[28025:28140] <- 0
  cl_plus[84244:84359] <- 0
  cl_plus[269312:269427] <- 0
  cl_minus[1331107:1331222] <- 0
  cl_minus[1578034:1578149] <- 0
  
  cvg_plus[22024:22139] <- df_5S$cvg_sum
  cvg_plus[28025:28140] <- 0
  cvg_plus[84244:84359] <- 0
  cvg_plus[269312:269427] <- 0
  cvg_minus[1331107:1331222] <- 0
  cvg_minus[1578034:1578149] <- 0
  
  # Combine data for CL sites and coverage and normalize to the library size
  bam_df <- data.frame(
    nt_genome_coord = seq(1, length(genome_seq)),
    cvg_plus = cvg_plus / length(bam) * 1E6,
    cvg_minus = cvg_minus / length(bam) * 1E6,
    cl_plus = cl_plus / length(bam) * 1E6,
    cl_minus = cl_minus / length(bam) * 1E6
  )
  
  # Write the coverage on plus and minus strands to WIG files
  wig_plus_file <- paste("igv/cvg_", sample_name, "_plus.wig", sep = "")
  wig_plus_header <- paste('track type=wiggle_0 name=', sample_name, '_plus color=255,0,0',
                           sep = "")
  wig_plus <- bam_df |>
    select(nt_genome_coord, cvg_plus) |> 
    filter(cvg_plus != 0)
  write_lines(c(wig_plus_header, 'variableStep chrom=NC_002737.2'), wig_plus_file)
  write_delim(wig_plus, wig_plus_file, append = TRUE, col_names = FALSE)
  
  wig_minus_file <- paste("igv/cvg_", sample_name, "_minus.wig", sep = "")
  wig_minus_header <- paste('track type=wiggle_0 name=', sample_name, '_minus color=255,0,0',
                            sep = "")
  wig_minus <- bam_df |>
    select(nt_genome_coord, cvg_minus) |> 
    filter(cvg_minus != 0)
  write_lines(c(wig_minus_header, 'variableStep chrom=NC_002737.2'), wig_minus_file)
  write_delim(wig_minus, wig_minus_file, append = TRUE, col_names = FALSE)
  
  # Write the cross-links on plus and minus strands to WIG files
  wig_plus_file <- paste("igv/cl_", sample_name, "_plus.wig", sep = "")
  wig_plus_header <- paste('track type=wiggle_0 name=', sample_name, '_plus color=255,0,0',
                           sep = "")
  wig_plus <- bam_df |>
    select(nt_genome_coord, cl_plus) |> 
    filter(cl_plus != 0)
  write_lines(c(wig_plus_header, 'variableStep chrom=NC_002737.2'), wig_plus_file)
  write_delim(wig_plus, wig_plus_file, append = TRUE, col_names = FALSE)
  
  wig_minus_file <- paste("igv/cl_", sample_name, "_minus.wig", sep = "")
  wig_minus_header <- paste('track type=wiggle_0 name=', sample_name, '_minus color=255,0,0',
                            sep = "")
  wig_minus <- bam_df |>
    select(nt_genome_coord, cl_minus) |> 
    filter(cl_minus != 0)
  write_lines(c(wig_minus_header, 'variableStep chrom=NC_002737.2'), wig_minus_file)
  write_delim(wig_minus, wig_minus_file, append = TRUE, col_names = FALSE)
  
  # Add the data to the final data.frame
  bam_df <- bam_df |> 
    rename(
      cvg_plus = paste("cvg_", sample_name, "_plus", sep = ""),
      cvg_minus = paste("cvg_", sample_name, "_minus", sep = ""),
      cl_plus = paste("cl_", sample_name, "_plus", sep = ""),
      cl_minus = paste("cl_", sample_name, "_minus", sep = "")
    )
  
  iclip_df <- left_join(iclip_df, bam_df)
}

# Save the output data.frame
write_csv(iclip_df, file = "iclip.csv")
write_tsv(iclip_df, file = "iclip.tsv")
