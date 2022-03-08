library(readr)
library(dplyr)

files <- read_tsv(snakemake@input[['txt']], col_names = "bins")
files <- files %>%
  mutate(bins = gsub("_sigs", "", bins),
         bins = gsub("\\.sig", "\\.ffn", bins))

write_tsv(files, snakemake@output[['txt']], col_names = F)
