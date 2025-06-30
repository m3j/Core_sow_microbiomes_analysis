library(tidyverse)
library(phyloseq)

source("code/Functions/physeq_utils.R")

# get phyloseq object
load("data/physeq_objects_filt0.1.RData")

# aglommerate to species level
physeq_species <- tax_glom(ps, "Species")

# get prevalence
prevalence <- as.data.frame(get_prevalence(physeq_species)) %>% rownames_to_column("otu_id")

# get total abundance across the dataset
total_abundance <- as.data.frame(get_abundance(physeq_species))

# get taxa information from phyloseq object
physeq_tax <- as.data.frame(tax_table(ps)) %>% rownames_to_column("otu_id")

# create dataframe combining prevalence, abundance and taxa
all_otu_abu_prev <- prevalence %>%
    left_join(total_abundance, by = "otu_id") %>%
    left_join(physeq_tax, by = "otu_id") %>%
    rename("relative_abundance" = "abundance_%") %>%
    select(Domain, Phylum, Class, Order, Family, Genus, Species, otu_id, prevalence, "relative_abundance", abundance_reads) %>%
    arrange(., desc(relative_abundance))

# save as tsv
write_tsv(all_otu_abu_prev, "tables/table_S1.tsv")
