library(dplyr)
library(readr)

#lineages <- read_csv("inputs/gtdb-rs202.taxonomy.v2.csv")
lineages <- read_csv(snakemake@input[['lineages']])

gather <- read_csv(snakemake@input[['gather']]) %>%
#gather <- read_csv("outputs/genbank/gather_gtdb-rs202-genomic.x.genbank.gather.csv") %>%
  mutate(accession = gsub(" .*", "", name)) %>%
  left_join(lineages, by = c("accession" = "ident")) %>%
  select(accession, species) %>%
  #mutate(species = gsub(" ", "_", species),
  #       database = paste0("gtdb-rs202.", species, ".protein-k10.nodegraph"))
  mutate(species = gsub(" ", "_", species))

write_csv(gather, snakemake@output[['csv']])
