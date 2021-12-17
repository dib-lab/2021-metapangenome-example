library(pagoo)
library(readr)
library(ggplot2)

read_long_sketch_table_as_pagoo <- function(path, threshold = 2000){
  sketch_table <- read_csv(path)
  
  # filter samples that don't have enough k-mers for the species
  sketch_table_grp <- sketch_table %>%
    group_by(acc) %>%
    tally()

  print(sketch_table_grp)
  
  keep <- sketch_table_grp %>%
    filter(n > threshold)
  
  sketch_table <- sketch_table %>%
    filter(acc %in% keep$acc) %>%
    select(gene    = minhash, 
           org     = acc,
           cluster = minhash) 
  p <- pagoo(data = as.data.frame(sketch_table))
  return(p)
}

species_string1 <- snakemake@wildcards[['acc_db']]
species_string2 <- gsub("-", ".", species_string1)
species_string3 <- gsub(".*s__", "", species_string2)
species_string3 <- gsyb("_", " ", species_string3)

pg <- read_long_sketch_table_as_pagoo(snakemake@input[['csv']], threshold = 0)

pg_organisms <- gsub(paste0(".", species_string2), "", as.character(pg$organisms$org))

destfile <- "inputs/hmp2_metadata.csv"
url <- "https://ibdmdb.org/tunnel/products/HMP2/Metadata/hmp2_metadata.csv"
if (!file.exists(destfile)) {
  download.file(url, destfile, method="auto")
}
hmp_metadata <- read_csv(destfile)

h4017 <- hmp_metadata %>%
  select(data_id = "External ID", data_type,
         week_num, diagnosis, antibiotics = "Antibiotics") %>%
  filter(data_type == "metagenomics") %>%
  filter(data_id %in% c('HSM67VF9', 'HSM67VFD', 'HSM67VFJ', 'HSM6XRQB',
                        'HSM6XRQI', 'HSM6XRQK', 'HSM6XRQM', 'HSM6XRQO',
                        'HSM7CYY7', 'HSM7CYY9', 'HSM7CYYB', 'HSM7CYYD'))

h4017 <- h4017 %>%
  mutate(org = paste0(data_id, ".", species_string2)) %>%
  filter(data_id %in% pg_organisms)

pg$add_metadata(map = "org", as.data.frame(h4017))

pdf(snakemake@output[['pca']], )
pg$gg_pca() +
  geom_point(aes(color = antibiotics)) +
  theme_minimal() +
  theme(plot.title = element_text(face = "italic")) +
  labs(title = species_string3) +
  scale_color_manual(values = c("#1F78B4", "#33A02C"))
dev.off()

# try color binmap --------------------------------------------------------

tpm <- t(pg$pan_matrix)
tpm[which(tpm > 0, arr.ind = TRUE)] <- 1L
bm <- as.data.frame(tpm)
or <- order(rowSums(bm), decreasing = TRUE)
lvls <- rownames(bm)[or]
bm$Cluster <- factor(rownames(bm), levels = lvls)
bm <- reshape2::melt(bm, 'Cluster')
bm$value <- factor(bm$value, levels = c(1, 0))
bm <- left_join(bm, h4017, by = c("variable" = "org"))
bm$value_abx <- paste0(bm$value, "_", bm$antibiotics)
bm$value_abx <- factor(bm$value_abx, levels = c("0_No", "1_No", "0_Yes",  "1_Yes"))
colnames(bm)[which(colnames(bm) == 'variable')] <- "Organism"

pdf(snakemake@output[['binmap']], height = 3, width = 4)
ggplot(bm, aes(Cluster, as.factor(week_num), fill=value_abx)) +
  geom_raster() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(face = "italic")) +
  #scale_fill_grey(start = .2, end = .9)
  scale_fill_brewer(palette = "Paired", labels = c("0, No", "1, No", "0, Yes", "1, Yes")) +
  labs(x = "protein k-mer", y = "week number", fill = "presence\n& antibiotics",
       title = species_string3)
dev.off()

# Summary of pagoo results
# Phocaeicola vulgatus
# core genome is stable across abx administration, but set of genes disappears with second round of abx.
# third round of abx wipes it out.
# ignoring the week 37 sample:
# pg$summary_stats
# Category    Number
# 1       Total     41005
# 2        Core     12437
# 3       Shell      8376
# 4       Cloud     20192

# Bacteroides uniformis also impacted by second round of antibiotics

# clostridium_Q symbiosium is present during first ABX, but then disappears mostly

# E. bolteae disturbance succession
#*C. bolteae* is a member of the normal gut microbiota but is an opportunistic pathogen that exploits compromised intestinal barriers [@doi:10.1186/s12864-016-3152-x].
#It is associated with disturbance succession and has increased gene expression during gut dysbiosis [@doi:10.1101/gr.138198.112; @doi:10.1038/s41586-019-1237-9].

# Parabacteroides merdae blooms during second course of ABX, but is wiped out by the third