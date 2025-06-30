library(tidyverse)
library(phyloseq)
library(vegan)

source("code/Functions/stats_utils.R")
source("code/Functions/dist_utils.R")

# get phyloseq object
load("data/physeq_objects_filt0.1.RData")

# get robust aitchison distance matrix of phyloseq object
r_aitchison_dist <- robust_ait(ps)

# get metadata from phyloseq object
metadata <- data.frame(sample_data(ps))

# collect relevant features
relevant_cols <- metadata %>%
    dplyr::select(c(
        "GD109_temp_gestation_unit", "LD21_temp_lactation_unit", "parity", "no_piglets_alive",
        "no_stillborn_piglets", "no_piglets_weaned_per_year_sow_not_for_gilts", "sampling_stage",
        "Country", "farm_number", "probiotic_any_binary", "antibiotic_binary", "feed_combined", "runid"
    ))

# collect those that have a corresponding index
has_coresponding_index <- metadata %>%
    dplyr::select(ends_with("_index")) %>%
    colnames() %>%
    sub("_index", "", .)

# initiate empty dataframe for dispersion results
dispersion_results <- data.frame()

# loop through relevant features
for (x in seq_along(relevant_cols)) {
    feature <- colnames(relevant_cols[x])

    # if feature has a corresponding index use this column
    if (feature %in% has_coresponding_index) {
        feature_groups <- paste(feature, "_index", sep = "")
    } else {
        feature_groups <- feature
    }

    # calculate dispersion
    dispersion <- vegan::betadisper(r_aitchison_dist,
        group = as.vector(metadata[[feature_groups]]),
        type = "median",
        bias.adjust = TRUE
    )

    # calculate anova statistict on the dispersion
    anova_disp <- as.data.frame(anova(dispersion))

    # add the feature name to the results
    anova_disp$feature <- feature

    # add the results to the result dataframe
    dispersion_results <- rbind(anova_disp, dispersion_results)
}

# calculate adjusted p values with BH-fdr and save df - dispersion
# remove the residual rows
p_values_disp <- na.omit(dispersion_results[5])
# calculate the adjusted p values
p_adj_disp <- p.adjust(p_values_disp[, 1], "fdr")

# add adj. p-values to the dataframe
p_values_disp$p_adj <- p_adj_disp

# make rownames to column in both dataframes to be able to join on "index"
p_values_disp <- rownames_to_column(p_values_disp, "index")
dispersion_results$index <- rownames(dispersion_results)

# join results and adj. p-values
dispersion_results <- merge(dispersion_results, p_values_disp %>% select(index, p_adj),
    by = "index", all.x = TRUE
)

# get the number indicating the model and residual pairs
dispersion_results$pair_number <- as.numeric(gsub("[^0-9]", "", dispersion_results$index))

# order on pair number such that the model row is followed by the residual row
dispersion_results_sorted <- dispersion_results[order(dispersion_results$pair_number), ]

# remove paired number again
dispersion_results_clean <- dispersion_results_sorted %>% select(-c("pair_number"))

# adjust number of decimals
dispersion_results_clean$"Sum Sq" <- round(dispersion_results_clean$"Sum Sq", 2)
dispersion_results_clean$"Mean Sq" <- round(dispersion_results_clean$"Mean Sq", 2)
dispersion_results_clean$"F value" <- round(dispersion_results_clean$"F value", 3)
dispersion_results_clean$"Pr(>F)" <- round(dispersion_results_clean$"Pr(>F)", 4)
dispersion_results_clean$"p_adj" <- round(dispersion_results_clean$"p_adj", 4)

# move feature column up as the first
dispersion_results_clean <- dispersion_results_clean %>% relocate(feature, .before = 1)

# save table
write_tsv(as.data.frame(dispersion_results_clean), "tables/table_S5.tsv")
