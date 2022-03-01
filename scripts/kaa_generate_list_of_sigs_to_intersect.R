library(readr)
library(dplyr)

#path <- "outputs/nbhd_sketch_tables_species/GCA_000162535.1-s__Parabacteroides_distasonis_long.csv"
#sketch_table <- read_csv(path)

sketch_table <- read_csv(snakemake@input[['csv']])
# filter samples that don't have enough k-mers for the species
# (using 10000 as threshold)
sketch_table_grp <- sketch_table %>%
  group_by(acc) %>%
  tally()

keep <- sketch_table_grp %>%
  filter(n > 10000)

pg_sigs <- sketch_table %>%
  filter(acc %in% keep$acc) %>%
  select(acc) %>%
  distinct() %>%
  mutate(acc = gsub("\\.G", "-G", acc), 
         acc = gsub("\\.s", "-s", acc),
         acc = paste0("outputs/nbhd_sigs_species/", acc, ".sig"))

write_tsv(pg_sigs, snakemake@output[['lst']], col_names = FALSE)
