library(tidyverse)
library(phyloseq)

source("code/Functions/physeq_utils.R")

# get phyloseq object
load("data/physeq_objects_filt0.1.RData")

# get the species distribution of a single otu given name at a certain rank across a specific grouping variable.
species_distributions <- function(physeq_obj, rank, name, grouping_var, post_fix = "", with_proportion = FALSE) {
    # change phyloseq counts to samplewise relative abundance in percentages
    physeq_obj_rel_abundance_perc <- phyloseq::transform_sample_counts(physeq_obj, function(x) ((x / sum(x)) * 100))

    # get the otu ids of the species in question
    species_otu <- as.data.frame(tax_table(physeq_obj_rel_abundance_perc)) %>%
        rownames_to_column("otu") %>%
        filter(.data[[rank]] == name) %>%
        pull("otu")

    # get the abundance of the otus
    species_counts <- data.frame(otu_table(physeq_obj_rel_abundance_perc))[species_otu, ]

    # rename sample id
    species_counts <- species_counts %>%
        rename_with(~ stringr::str_replace(., "X", "")) %>%
        rename_with(~ stringr::str_replace_all(., "\\.", "-"))

    # transpose the species abundances and make sample id a column and create a sum column per otu
    counts_transposed <- as.data.frame(t(species_counts)) %>%
        rownames_to_column("sampleid") %>%
        mutate(OTU_Sum = rowSums(select(., starts_with("OTU"))))

    # get metadata from phyloseq object and make a sample id column from index
    metadata <- data.frame(sample_data(physeq_obj)) %>% rownames_to_column("sampleid")

    # add grouping variable values to abundance data
    species_group <- metadata %>%
        select(all_of(grouping_var), sampleid) %>%
        left_join(counts_transposed, by = "sampleid")

    # get the mean value of otu abundance per group variable and rename including the post_fix
    species_dist_table <- species_group %>%
        group_by(.data[[grouping_var]]) %>%
        summarise(mean_abundance = mean(OTU_Sum)) %>%
        rename("{paste0('mean_abundance', post_fix)}" := mean_abundance)

    # if propotionality is selected the porportion of the abundance across the grouping variable will be added as a seperate coulmn
    if (with_proportion) {
        mean_abundance_col <- paste0("mean_abundance", post_fix)
        with_propotions <- species_dist_table %>% mutate(proportion = !!sym(mean_abundance_col) / sum(!!sym(mean_abundance_col)) * 100)
        return(with_proportion)
    }

    return(species_dist_table)
}

# will calculate species distribution of 4 predefined phyloseq objects grouping on country and return a combinded df
species_count_comparison <- function(species_name) {
    # get species distributions
    species_ges_good <- species_distributions(physeq_obj = physeq_ges_good, rank = "Species", name = species_name, grouping_var = "Country", post_fix = "_ges_good")

    species_lac_good <- species_distributions(physeq_obj = physeq_lac_good, rank = "Species", name = species_name, grouping_var = "Country", post_fix = "_lac_good")

    species_ges_bad <- species_distributions(physeq_obj = physeq_ges_bad, rank = "Species", name = species_name, grouping_var = "Country", post_fix = "_ges_bad")

    species_lac_bad <- species_distributions(physeq_obj = physeq_lac_bad, rank = "Species", name = species_name, grouping_var = "Country", post_fix = "_lac_bad")

    # combine into one dataframe with species as rows
    combined_counts <- species_ges_bad %>%
        left_join(species_ges_good, by = "Country") %>%
        left_join(species_lac_bad, by = "Country") %>%
        left_join(species_lac_good, by = "Country")

    return(combined_counts)
}

# agglommerate to species level
physeq_species <- phyloseq::tax_glom(ps, "Species")

# devide species phyloseq object on sampling stage
physeq_divided <- split_physeq_in_two(physeq_species, "sampling_stage", "L21", "G109")
physeq_lac <- physeq_divided[[1]]
physeq_ges <- physeq_divided[[2]]

# define phyloseq objects based sampling stage and performance
physeq_lac_good <- microViz::ps_filter(physeq_lac, no_piglets_weaned_per_year_sow_not_for_gilts_index2 == "Good>37.1")
physeq_ges_good <- microViz::ps_filter(physeq_ges, no_piglets_weaned_per_year_sow_not_for_gilts_index2 == "Good>37.1")
physeq_ges_bad <- microViz::ps_filter(physeq_ges, no_piglets_weaned_per_year_sow_not_for_gilts_index2 == "Bad<37.1")
physeq_lac_bad <- microViz::ps_filter(physeq_lac, no_piglets_weaned_per_year_sow_not_for_gilts_index2 == "Bad<37.1")

# get species abundances per predefined phyloseq objects
xerosis <- species_count_comparison(species_name = "Corynebacterium_xerosis")
suis <- species_count_comparison(species_name = "Streptococcus_suis")
claveliimonas <- species_count_comparison(species_name = "Claveliimonas_uncultured_bacterium")
tepidibacter <- species_count_comparison(species_name = "Tepidibacter_uncultured_Gram-positive")
pneumoniae <- species_count_comparison(species_name = "Streptococcus_pneumoniae")
coli <- species_count_comparison(species_name = "Escherichia_coli")
bacillus <- species_count_comparison(species_name = "Bacillus_uncultured_bacterium")

# create named list with abundances per species
species_tbls <- list(
    Corynebacterium_xerosis = xerosis,
    Streptococcus_suis = suis,
    Claveliimonas_uncultured_bacterium = claveliimonas,
    Tepidibacter_uncultured_Gram_positive = tepidibacter,
    Streptococcus_pneumoniae = pneumoniae,
    Escherichia_coli = coli,
    Bacillus_uncultured_bacterium = bacillus
)

# combine the lists into a dataframe
combind <- dplyr::bind_rows(species_tbls, .id = "species")

# reorder columns
combind_order <- combind %>%
    relocate(mean_abundance_ges_good, .before = 3) %>%
    relocate(mean_abundance_lac_good, .before = 5)

# save table
write_tsv(combind_order, "tables/table_S7.tsv")
