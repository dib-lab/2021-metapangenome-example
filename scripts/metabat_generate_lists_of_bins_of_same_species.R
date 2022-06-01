library(readr)
library(dplyr)
library(tidyr)
library(purrr)

gtdb_lineages <- read_csv(snakemake@input[['lineages']])

# this specified as a sysglob to not mess up the sample wildcard in the workflow
gather_bins <- Sys.glob("outputs/metabat2_gather/*csv") %>%
#gather_bins <- unlist(snakemake@input[['gather']]) %>%
  map_dfr(read_csv, col_types = c("dddddlllcccddddcccd")) %>%
  separate(col = name, into = c("accession"), remove = F, sep = " ", extra = "drop") %>%
  mutate(sample = basename(query_filename)) %>%
  mutate(sample = gsub("-.*", "", sample)) %>%
  left_join(gtdb_lineages, by = c("accession" = "ident"))

species_of_interest <- read_csv(unlist(snakemake@input[['acc_to_db']])) %>%
  #filter(!species %in% c("s__Ruminococcus_B_gnavus", "s__Clostridium_Q_symbiosum", "s__Roseburia_intestinalis")) %>%
  mutate(species = gsub("(s__.*?)_", "\\1 ", species))

acc_to_db <- read_csv(unlist(snakemake@input[['acc_to_db']])) %>%
  #filter(!species %in% c("s__Ruminococcus_B_gnavus", "s__Clostridium_Q_symbiosum", "s__Roseburia_intestinalis")) %>%
  mutate(acc_to_db = paste0(accession, "-", species)) %>%
  mutate(species = gsub("(s__.*?)_", "\\1 ", species))

# only consider best match for the bin
gather_bins <- gather_bins %>%
  filter(gather_result_rank == 0)

named_bins <- gather_bins %>% 
  left_join(acc_to_db, by = "species") %>%
  select(acc_to_db, query_filename) %>%
  mutate(sig_filename = gsub("metabat2", "metabat2_prokka_sigs", query_filename),
         sig_filename = gsub("/bin", "_bin", sig_filename),
         sig_filename = gsub(".fa", ".sig", sig_filename)) %>%
  select(acc_to_db, sig_filename)

# filter to single acc_db using wildcard
# and only write out sig path
tmp <- named_bins %>%
  filter(acc_to_db %in% snakemake@wildcards[["acc_db"]]) %>%
  select(sig_filename)

write_tsv(tmp, snakemake@output[['txt']], col_names = FALSE)

#named_bins %>%
#  group_by(acc_to_db) %>%
#  group_walk(~write_tsv(.x, paste0("outputs/metabat2_gather_labelled_bins/", .y$acc_to_db, ".txt"), col_names = FALSE))

