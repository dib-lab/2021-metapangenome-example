library(dplyr)
library(readr)
library(purrr)
library(tidyr)

gtdb_lineages <- read_csv(snakemake@input[['gtdb_lineages']])

charcoal_lineages <- Sys.glob(paste0("outputs/metabat2_gather/", snakemake@wildcards[['sample']], "_bin.*_gather_gtdb-rs202-genomic.csv")) %>%
#charcoal_lineages <- unlist(snakemake@input[['gather']]) %>% 
  map_dfr(read_csv, col_types = "dddddlllcccddddcccd") %>%
  filter(gather_result_rank == 0) %>%
  select(query_name, name) %>%
  mutate(query_name = gsub(".*\\.", "", query_name),
         query_name = gsub("bin", "bin.", query_name),
         query_name = paste0(query_name, ".fa"),
         accession = gsub(" .*", "", name)) %>%
  left_join(gtdb_lineages, by = c("accession" = "ident")) %>%
  select(-name, -accession)

write_csv(charcoal_lineages, snakemake@output[['charcoal_lineages']])
