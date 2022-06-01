library(dplyr)
library(readr)

read_names <- read_tsv(snakemake@input[['names']], col_names = c("names")) %>%
  mutate(names = gsub("^@[0-9]*:[0-5]:[0-5]:", "", names))

write_tsv(read_names, snakemake@output[["lst"]], col_names = FALSE)
