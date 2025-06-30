library(tidyverse)
library(phyloseq)

# get phyloseq object
load("data/physeq_objects_filt0.1.RData")

# get metadata from phyloseq object
metadata <- data.frame(sample_data(ps))

# get country counts distributed on performance categories
performance_distribution <- table(metadata$no_piglets_weaned_per_year_sow_not_for_gilts_index, metadata$Country)

# create a sum column and make rownames to a column
performance_dist_total <- performance_distribution %>%
    as.data.frame.matrix() %>%
    mutate(Total = rowSums(across(where(is.numeric)))) %>%
    rownames_to_column("performance")

# reset performance values
performance_dist_total$performance <- c("Good", "Medium", "Bad")

# save table
write_tsv(performance_dist_total, "tables/table_S6.tsv")
