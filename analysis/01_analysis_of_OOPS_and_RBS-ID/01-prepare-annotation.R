library(stringr)

# Read the modified Proteome annotation file
pa <- read.delim("Raw data/Proteome annotation.tsv")
pa <- within(pa, {
  UniProt <- str_sub(UniProt, 1, 6)
  ENSG <- str_replace(ENSG, ";.+", "")
})

# Annotate functional categories
pa <- within(pa, {
  glycolysis <- str_detect(KEGG.name, "Glycolysis")
  #glycolysis[UniProt == "Q99YC4"] <- FALSE
  #glycolysis[UniProt == "Q9A0V5"] <- FALSE
  glycolysis[UniProt == "Q99XX0"] <- FALSE
  glycolysis[UniProt == "Q99ZS0"] <- FALSE
  glycolysis[UniProt == "Q99ZX5"] <- FALSE
  glycolysis[UniProt == "Q99ZX6"] <- FALSE
  glycolysis[UniProt == "Q99ZX7"] <- FALSE
  glycolysis[UniProt == "Q99ZX8"] <- FALSE
  glycolysis[UniProt == "Q9A0X3"] <- FALSE
  glycolysis[UniProt == "Q9A1X7"] <- FALSE
  
  rna.binding <- str_detect(GOMF.name, "RNA binding")
  rna.binding[UniProt == "P0C0F1"] <- TRUE    # cspA
  rna.binding[UniProt == "P68893"] <- TRUE    # nusG antiterminator
  rna.binding[UniProt == "P64285"] <- TRUE    # greA
  rna.binding[UniProt == "Q9A1J2"] <- TRUE    # jag
  rna.binding[UniProt == "Q9A0C1"] <- TRUE   # Putative KH domain protein
  rna.binding[UniProt == "Q9A204"] <- TRUE   # Putative S4 domain protein
  rna.binding[UniProt == "Q9A206"] <- TRUE    # Peptidyl-tRNA hydrolase
  rna.binding[UniProt == "P63999"] <- TRUE    # D-aminoacyl-tRNA deacylase
  
  aatrna.synthesis <- str_detect(KEGG.name, "Aminoacyl-tRNA biosynthesis")
  
  ribosome <- str_detect(KEGG.name, "Ribosome")
  ribosome[UniProt == "Q99YN9"] <- TRUE    # Ribosome hibernation promoting factor hpf
  ribosome[UniProt == "P68903"] <- TRUE    # rrf
  
  ribosome.biogenesis <- str_detect(GOBP.name, "ribonucleoprotein complex biogenesis")
  ribosome.biogenesis[UniProt == "Q9A1F4"] <- TRUE    # GTPase yqeH
  ribosome.biogenesis[UniProt == "Q99Z94"] <- TRUE    # GTPase obg
  ribosome.biogenesis[UniProt == "Q9A088"] <- TRUE    # GTPase engB
  ribosome.biogenesis[UniProt == "Q9A1H9"] <- TRUE    # GTPase rsgA
  ribosome.biogenesis[UniProt == "P65971"] <- TRUE    # rbfA
  
  rna.modif <- str_detect(GOMF.name, "RNA methyltransferase activity")
  rna.modif[UniProt == "Q99ZH2"] <- TRUE    # SPy_1232 rRNA methyltransferase
  rna.modif[UniProt == "Q9A1A5"] <- TRUE    # rRNA methylase
  rna.modif[UniProt == "Q99YU1"] <- TRUE    # 16S rRNA methyltransferase
  rna.modif[UniProt == "Q99XI8"] <- TRUE    # tRNA modification enzyme MnmG
  rna.modif[UniProt == "Q7DAN3"] <- TRUE    # tRNA-dihydrouridine synthase
  rna.modif[UniProt == "Q9A1E8"] <- TRUE    # tmcAL
  rna.modif[UniProt == "Q99Y44"] <- TRUE    # tRNA N6-adenosine threonylcarbamoyltransferase
  rna.modif[UniProt == "Q99YL5"] <- TRUE    # putative N6-adenine-specific DNA methylase
  rna.modif[UniProt == "P58075"] <- TRUE    # mnmA
  rna.modif[UniProt == "Q99XI8"] <- TRUE    # mnmG
  rna.modif[UniProt == "Q9A0A5"] <- TRUE    # CCA-adding_enz_firmicutes
  rna.modif[UniProt == "Q9A0D8"] <- TRUE    # thiI
  rna.modif[UniProt == "Q9A1E8"] <- TRUE    # tRNA(Met) cytidine acetate ligase
  rna.modif[UniProt == "Q7DAN3"] <- TRUE    # tRNA-dihydrouridine synthase
  rna.modif[UniProt == "Q99XS3"] <- TRUE    # PseudoU_synth_2
  rna.modif[UniProt == "P65850"] <- TRUE    # PseudoU_synth_1
  rna.modif[UniProt == "Q99YX6"] <- TRUE    # 23S rRNA (cytidine1920-2'-O)-methyltransferase
  rna.modif[UniProt == "P66021"] <- TRUE    # PseudoU_synth_2;S4
  rna.modif[UniProt == "Q99Z92"] <- TRUE    # PseudoU_synth_2;S4
  rna.modif[UniProt == "P65858"] <- TRUE    # tRNA_psdUridine_synth_TruB
  rna.modif[UniProt == "Q99ZP3"] <- TRUE    # L-threonylcarbamoyladenylate synthase
  rna.modif[UniProt == "Q99ZQ6"] <- TRUE    # PseudoU_synth_2
  rna.modif[UniProt == "Q9A0D1"] <- TRUE    # PseudoU_synth_2
  rna.modif[UniProt == "Q9A1B0"] <- TRUE    # PseudoU_synth_2
  rna.modif[UniProt == "Q99Z51"] <- TRUE    # PsdUridine_synth_RsuA/RluD
  rna.modif[UniProt == ""] <- TRUE    # 
  rna.modif[UniProt == ""] <- TRUE    # 
  rna.modif[UniProt == ""] <- TRUE    #
  rna.modif[UniProt == ""] <- TRUE    # 
  rna.modif[UniProt == ""] <- TRUE    # 
  rna.modif[UniProt == ""] <- TRUE    # 
  rna.modif[UniProt == ""] <- TRUE    # 
  rna.modif[UniProt == ""] <- TRUE    #
  
  rna.pol <- str_detect(KEGG.name, "RNA polymerase")
  
  rnase <- str_detect(GOMF.name, "ribonuclease activity")
  rnase[UniProt == "Q9A1H5"] <- TRUE # 3'-5' exonuclease YhaM
  rnase[UniProt == "P66670"] <- TRUE # Ribonuclease 3
  rnase[Gene.name == "SPy_1985"] <- FALSE
  #rnase[Gene.name == "cas9"] <- FALSE
  rnase[Gene.name == "hsdM"] <- FALSE
  rnase[Gene.name == "hsdR"] <- FALSE
  rnase[Gene.name == "uvrA"] <- FALSE
  rnase[Gene.name == "uvrB"] <- FALSE
  rnase[Gene.name == "uvrC"] <- FALSE
  rnase[Gene.name == "xseA"] <- FALSE
  rnase[Gene.name == "xseB"] <- FALSE
  
  translation <- rep(FALSE, nrow(pa))
  translation[UniProt == "P68903"] <- TRUE    # RRF
  translation[UniProt == "Q99XQ7"] <- TRUE    # EF-TS
  translation[UniProt == "P68773"] <- TRUE    # EF-P
  translation[UniProt == "Q99YG1"] <- TRUE    # IF-2
  translation[UniProt == "P66021"] <- TRUE    # PrfC
  translation[UniProt == "Q99ZK1"] <- TRUE    # Signal recognition particle protein ffh
  translation[UniProt == "Q99ZP5"] <- TRUE    # PrfA
  translation[UniProt == "Q99ZV8"] <- TRUE    # Elongation factor 4
  translation[UniProt == "P65146"] <- TRUE    # IF-3
  translation[UniProt == "Q9A0S5"] <- TRUE    # PrfB
  translation[UniProt == "P69952"] <- TRUE    # EF-Tu
  translation[UniProt == "Q9A126"] <- TRUE    # SmpB
  translation[UniProt == "P69946"] <- TRUE    # EF-G
  translation[UniProt == "P65123"] <- TRUE    # IF-1
  translation[UniProt == "Q9A206"] <- TRUE    # Peptidyl-tRNA hydrolase pth
  translation[UniProt == ""] <- TRUE    #
  translation[UniProt == ""] <- TRUE    #
  translation[UniProt == ""] <- TRUE    #
  translation[UniProt == ""] <- TRUE    #
  
  rbp <- rna.binding | aatrna.synthesis | ribosome |  
    ribosome.biogenesis | rna.modif | rna.pol | rnase | translation
})

pa$rbp.type <- "-"
pa <- within(pa, {
  rbp.type[rbp] <- "Other RBP"
  rbp.type[rnase] <- "RNase"
  rbp.type[rna.pol] <- "RNA polymerase"
  rbp.type[rna.modif] <- "RNA modification"
  rbp.type[ribosome.biogenesis] <- "Ribosome biogenesis"
  rbp.type[ribosome] <- "Ribosomal protein"
  rbp.type[aatrna.synthesis] <- "tRNA ligase"
  rbp.type[translation] <- "Translation"
})

# Fix some Gene.names
pa <- within(pa, {
  Gene.name[UniProt == "Q9A1H5"] <- "yhaM"
  Gene.name[UniProt == "Q9A066"] <- "rpsA"
  Gene.name[UniProt == "Q9A0Y2"] <- "asnS"
  Gene.name[UniProt == "Q99Y00"] <- "mrnC"
  Gene.name[UniProt == "Q99Y42"] <- "rnj"
  Gene.name[UniProt == "Q99ZY2"] <- "rnj"
  Gene.name[UniProt == "Q9A1I2"] <- "tatD"
})

write.table(pa, "Proteome annotation with RBP.tsv", row.names = FALSE,
            quote = FALSE, sep = "\t")