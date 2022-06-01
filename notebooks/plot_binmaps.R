library(pagoo)
library(dplyr)
library(readr)
library(ggplot2)
library(aplot)

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

destfile <- "inputs/hmp2_metadata.csv"
url <- "https://ibdmdb.org/tunnel/products/HMP2/Metadata/hmp2_metadata.csv"
if (!file.exists(destfile)) {
  download.file(url, destfile, method="auto")
}
hmp_metadata <- read_csv(destfile)

# serology data_type contains information on abx types
# tmp <- hmp_metadata %>%
#   select(participant_id = "Participant ID", data_id = "External ID", data_type,
#          week_num, diagnosis, antibiotics = "Antibiotics",
#          fecalcal, metronidazole = "Flagyl (Metronidazole)",
#          cipro = "Cipro (Ciprofloxin)", rifaxamin = "Xifaxin (rifaxamin)",
#          levaquin = Levaquin, other_abx = "Other Antibiotic:") %>%
#   filter(participant_id =="H4017")

h4017 <- hmp_metadata %>%
  select(data_id = "External ID", data_type,
         week_num, diagnosis, antibiotics = "Antibiotics") %>%
  filter(data_type == "metagenomics") %>%
  filter(data_id %in% c('HSM67VF9', 'HSM67VFD', 'HSM67VFJ', 'HSM6XRQB',
                        'HSM6XRQI', 'HSM6XRQK', 'HSM6XRQM', 'HSM6XRQO',
                        'HSM7CYY7', 'HSM7CYY9', 'HSM7CYYB', 'HSM7CYYD'))



# gather results ----------------------------------------------------------

destfile <- "inputs/gtdb-rs202.taxonomy.v2.csv"
url <- "https://osf.io/p6z3w/download"
if (!file.exists(destfile)) {
  download.file(url, destfile, method="auto") 
}
gtdb_lineages <- read_csv(destfile)

gather_results <- Sys.glob("outputs/sample_gather/*genomic.csv") %>%
  set_names() %>%
  map_dfr(read_csv, col_types = c("dddddlllcccddddcccd"), .id = "sample") %>%
  mutate(sample = gsub("outputs/sample_gather/", "", sample)) %>%
  mutate(sample = gsub("_gather_gtdb-rs202-genomic.csv", "", sample)) %>%
  separate(col = name, into = c("accession"), remove = F, sep = " ", extra = "drop")

gather_results <- left_join(gather_results, gtdb_lineages, by = c("accession" = "ident")) %>%
  left_join(h4017, by = c("sample" = "data_id"))

# filter to sgc species
common_species <- gather_results %>%
  select(sample, species) %>%
  distinct() %>%
  group_by(species) %>%
  tally() %>%
  filter(n == 12)
 
sgc_species <- gather_results %>%
  filter(species %in% common_species$species) %>%
  group_by(species) %>%
  summarize(sum_f_unique_to_query = sum(f_unique_to_query)) %>%
  filter(sum_f_unique_to_query > 0.2)

# get per species, per metagenome fractional abundances of each organism
# also get base pairs?
gather_summarized <- gather_results %>%
  filter(species %in% sgc_species$species) %>%
  group_by(sample, species, week_num) %>%
  summarize(sum_f_unique_to_query = sum(f_unique_to_query),
            sum_intersect_bp = sum(intersect_bp))

# add in information about the number of bins per species per sample
View(gather_summarized)
bin_info <- Sys.glob("outputs/metabat2_gather_labelled_bins/*txt") %>%
  set_names() %>%
  map_dfr(read_tsv, col_names = "bin_name", .id = "species") %>%
  mutate(bin_name = gsub(".sig", "", basename(bin_name)),
         species =  gsub(".*-", "", basename(species)),
         species =  gsub(".txt", "", species),
         species =  gsub("_", " ", species), 
         species =  gsub("  ", "__", species)) %>%
  separate(bin_name, into = c("sample", "bin"), sep = "_") %>%
  group_by(sample, species) %>%
  tally() %>%
  select(sample, species, n_bins = n)

gather_summarized <- left_join(gather_summarized, bin_info) %>%
  mutate(n_bins = replace_na(n_bins, 0))

# function plot binmap --------------------------------------------------------
plot_binmap <- function(species_string, metadata = h4017, gather_summarized_species, label_x_axis = F) {
  #metadata <- h4017
  #gather_summarized_species = gather_summarized %>% filter(species == "s__Bacteroides fragilis")
  #species_string <- "GCF_003458955.1-s__Bacteroides_fragilis"
  species_string1 <- species_string
  species_string2 <- gsub("-", ".", species_string1)
  species_string3 <- gsub(".*s__", "", species_string2)
  species_string3 <- gsub("_", " ", species_string3)
  species_string4 <- gsub(".*\\.", "", species_string2)
  species_string4 <- gsub("(.*__.*?)_(.*?)", "\\1 \\2 ", species_string4)
  species_string4 <- gsub("  ", " ", species_string4)
  
  pg <- read_long_sketch_table_as_pagoo(paste0("outputs/nbhd_sketch_tables_species/", 
                                               species_string1, "_long.csv"), threshold = 0)
  pg_organisms <- gsub(paste0(".", species_string2), "", as.character(pg$organisms$org))
  metadata <- metadata %>%
    mutate(org = paste0(data_id, ".", species_string2)) %>%
    filter(data_id %in% pg_organisms)
  pg$add_metadata(map = "org", as.data.frame(metadata))
  
  tpm <- t(pg$pan_matrix)
  tpm[which(tpm > 0, arr.ind = TRUE)] <- 1L
  bm <- as.data.frame(tpm)
  or <- order(rowSums(bm), decreasing = TRUE)
  lvls <- rownames(bm)[or]
  bm$Cluster <- factor(rownames(bm), levels = lvls)
  bm <- reshape2::melt(bm, 'Cluster')
  bm$value <- factor(bm$value, levels = c(1, 0))
  bm <- left_join(bm, metadata, by = c("variable" = "org"))
  bm$value_abx <- paste0(bm$value, "_", bm$antibiotics)
  bm$value_abx <- factor(bm$value_abx, levels = c("0_No", "1_No", "0_Yes",  "1_Yes"))
  colnames(bm)[which(colnames(bm) == 'variable')] <- "Organism"
  
  binmap <- ggplot(bm, aes(Cluster, as.factor(week_num), fill=value_abx)) +
    geom_raster() +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(face = "italic", size = 8),
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7),
          legend.position = "none") +
    scale_fill_brewer(palette = "Paired", labels = c("0, No", "1, No", "0, Yes", "1, Yes")) +
    labs(x = "protein k-mer", y = "week number", fill = "presence\n& antibiotics",
         title = species_string3)
  print("bin map plotted")
   
  bp <- ggplot(gather_summarized_species, 
               aes(x = sum_intersect_bp, y = as.factor(week_num))) + 
    geom_col() +
    theme_minimal() +
    xlim(0, 139200000) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7)) +
    labs(x = "base pairs", y = "week number")
  print("bp plotted")
  
  frac <- ggplot(gather_summarized_species, 
                 aes(x = sum_f_unique_to_query, y = as.factor(week_num))) + 
    geom_col() +
    theme_minimal() +
    xlim(0, 0.31) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7)) +
    labs(x = "fraction of metagenome", y = "week number")
  print("frac plotted")
  
  n_bins <- ggplot(gather_summarized_species, 
                   aes(x = n_bins, y = as.factor(week_num))) + 
    geom_col() +
    theme_minimal() +
    xlim(0, 4) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7)) +
    labs(x = "number of bins", y = "week number")
  print("bins plotted")
    
  if(label_x_axis == F){
    binmap <- binmap + theme(axis.title.x = element_blank())
    bp <- bp + theme(axis.title.x = element_blank())
    frac <- frac + theme(axis.title.x = element_blank())
    n_bins <- n_bins + theme(axis.title.x = element_blank())
  }
  
  plt <- bp %>% aplot::insert_left(binmap) %>% insert_right(frac) %>% insert_right(n_bins)
  return(plt)
}

# apply function to different species 
list.files("outputs/nbhd_sketch_tables_species/")


bf <- plot_binmap(species_string = "GCF_003458955.1-s__Bacteroides_fragilis",
                  gather_summarized_species = gather_summarized %>% filter(species == "s__Bacteroides fragilis"))
bu <- plot_binmap(species_string = "GCF_009020325.1-s__Bacteroides_uniformis",
                  gather_summarized_species = gather_summarized %>% filter(species == "s__Bacteroides uniformis"))
eb <- plot_binmap(species_string = "GCF_003433765.1-s__Enterocloster_bolteae",
                  gather_summarized_species = gather_summarized %>% filter(species == "s__Enterocloster bolteae"))
pd <- plot_binmap(species_string = "GCA_000162535.1-s__Parabacteroides_distasonis",
                  gather_summarized_species = gather_summarized %>% filter(species == "s__Parabacteroides distasonis"))
pm <- plot_binmap(species_string = "GCF_003475305.1-s__Parabacteroides_merdae", label_x_axis = T,
                  gather_summarized_species = gather_summarized %>% filter(species == "s__Parabacteroides merdae"))
pv <- plot_binmap(species_string = "GCF_009025805.1-s__Phocaeicola_vulgatus", label_x_axis = T,
                  gather_summarized_species = gather_summarized %>% filter(species == "s__Phocaeicola vulgatus"))

png("tmp.png", height = 5, width = 10, res = 300, units = "in")
plot_list(bf, bu, eb, pd, pm, pv, nrow = 3, ncol = 2, heights = c(1, 1, 1.1))
dev.off()



# text --------------------------------------------------------------------

binmap <- ggplot(bm, aes(Cluster, as.factor(week_num), fill=value_abx)) +
  geom_raster() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(face = "italic", size = 8),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        legend.position = "bottom") +
  scale_fill_brewer(palette = "Paired", labels = c("0, No", "1, No", "0, Yes", "1, Yes")) +
  labs(x = "protein k-mer", y = "week number", fill = "presence & \n antibiotics",
       title = species_string3)

png("tmp.png", height = 6, width = 9, res = 300, units = "in")
plot_list(plot_list(bf, bu, eb, pd, pm, pv, nrow = 3, ncol = 2, heights = c(1, 1, 1.1)), 
          as_ggplot(get_legend(binmap)), nrow = 2, heights = c(7, 1))
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