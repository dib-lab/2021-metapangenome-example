library(dplyr)
library(readr)
library(purrr)
library(ggplot2)
library(tidyr)
library(aplot)

# read in metadata and results --------------------------------------------

destfile <- "inputs/gtdb-rs202.taxonomy.v2.csv"
url <- "https://osf.io/p6z3w/download"
if (!file.exists(destfile)) {
  download.file(url, destfile, method="auto") 
}
gtdb_lineages <- read_csv(destfile)

destfile <- "inputs/hmp2_metadata.csv"
url <- "https://ibdmdb.org/tunnel/products/HMP2/Metadata/hmp2_metadata.csv"
if (!file.exists(destfile)) {
  download.file(url, destfile, method="auto") 
}
hmp_metadata <- read_csv(destfile)
h4017 <- hmp_metadata %>%
  select(participant_id = "Participant ID", data_id = "External ID", data_type, 
         week_num, diagnosis, antibiotics = "Antibiotics") %>%
  filter(data_type == "metagenomics") %>%
  filter(data_id %in% c('HSM67VF9', 'HSM67VFD', 'HSM67VFJ', 'HSM6XRQB',
                        'HSM6XRQI', 'HSM6XRQK', 'HSM6XRQM', 'HSM6XRQO',
                        'HSM7CYY7', 'HSM7CYY9', 'HSM7CYYB', 'HSM7CYYD'))

#gather_results <- unlist(snakemake@input[['gather']]) %>%
gather_results <- Sys.glob("outputs/sample_gather/*genomic.csv") %>%
  set_names() %>%
  map_dfr(read_csv, col_types = c("dddddlllcccddddcccd"), .id = "sample") %>%
  mutate(sample = gsub("outputs/sample_gather/", "", sample)) %>%
  mutate(sample = gsub("_gather_gtdb-rs202-genomic.csv", "", sample)) %>%
  separate(col = name, into = c("accession"), remove = F, sep = " ", extra = "drop")

gather_results <- left_join(gather_results, gtdb_lineages, by = c("accession" = "ident")) %>%
  left_join(h4017, by = c("sample" = "data_id"))

# plot --------------------------------------------------------------------

abx_plt <- ggplot(h4017, aes(x = week_num, y = antibiotics)) +
  geom_point() +
  theme_minimal() +
  theme(axis.title.x = element_blank())+
  labs(x = "week number")

# library(ggthemes)
# common_phyla <- gather_results %>%
#   group_by(phylum) %>%
#   tally() %>%
#   filter(n >= 10)
# phylum_plt <- ggplot(gather_results %>%
#          mutate(phylum2 = ifelse(phylum %in% common_phyla$phylum, phylum, "other")),
#        aes(x = week_num, y = f_unique_to_query, fill = phylum2)) +
#   geom_col() +
#   theme_minimal() + 
#   ylim(0, 1) +
#   scale_fill_tableau(palette = "Tableau 20")
# 
# abx_plt %>% insert_bottom(phylum_plt)

# decide which genome to query with ---------------------------------------

common_species <- gather_results %>%
  select(sample, species) %>%
  distinct() %>%
  group_by(species) %>%
  tally() %>%
  filter(n == 12)

# common_species_plt <- ggplot(gather_results %>%
#                                filter(species %in% common_species$species),
#                              aes(x = week_num, y = f_unique_to_query)) +
#   geom_col() +
#   facet_wrap(~species) +
#   theme_minimal() 

sgc_species <- gather_results %>%
  filter(species %in% common_species$species) %>%
  group_by(species) %>%
  summarize(sum_f_unique_to_query = sum(f_unique_to_query)) %>%
  filter(sum_f_unique_to_query > 0.1)

gather_results2 <- gather_results %>%
  mutate(species2 = ifelse(species %in% sgc_species$species, species, "other")) %>%
  mutate(species2 = gsub("s__", "", species2))
gather_results2$species2 <- factor(gather_results2$species2, 
                                   levels =c('other', 
                                             'Bacteroides fragilis',
                                             'Bacteroides uniformis',
                                             'Clostridium_Q symbiosum',
                                             'Enterocloster bolteae',
                                             'Parabacteroides distasonis',
                                             'Parabacteroides merdae',
                                             'Phocaeicola vulgatus',    
                                             'Roseburia intestinalis', 
                                             'Ruminococcus_B gnavus'))
sgc_plt <- ggplot(gather_results2,
                     aes(x = week_num, y = f_unique_to_query, fill = species2)) +
  geom_col() +
  labs(x = "week number", y = "fraction of metagenome", fill = "species")+
  theme_minimal() +
  theme(legend.text = element_text(face = "italic")) +
  ylim(0, 1) + 
  scale_fill_brewer(palette = "Paired")

pdf("figures/common_species_breakdown.pdf", width = 6, height = 3)
#pdf(snakemake@output[['pdf']], width = 6, height = 3)
abx_plt %>% insert_bottom(sgc_plt)
dev.off()

# output sgc queries as gather file ---------------------------------------

# get rep accession
sgc_genomes <- gather_results %>%
  filter(species %in% sgc_species$species) %>%
  group_by(name, species) %>%
  summarize(sum_f_unique_to_query = sum(f_unique_to_query)) %>%
  ungroup() %>%
  group_by(species) %>%
  arrange(desc(sum_f_unique_to_query)) %>%
  slice(1)

gather_sgc_queries <- unlist(snakemake@input[['gather']])[6] %>%
#gather_sgc_queries <- Sys.glob("outputs/sample_gather/*genomic.csv")[6] %>%
  map_dfr(read_csv, col_types = c("dddddlllcccddddcccd")) %>%
  filter(name %in% sgc_genomes$name)

write_csv(gather_sgc_queries, snakemake@output[['gather_grist']])
