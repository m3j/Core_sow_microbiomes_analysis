library(phyloseq)

source("code/Functions/stats_utils.R")

# get phyloseq object
load("data/physeq_objects_filt0.1.RData")

# calculate pariwise anova for countries
country_starta_sampling <- calculate_pairwise_anova("Country", ps, strata = "sampling_stage")

# collect comparisons results and save tsv.
pairwise_res2tsv("tabel_S2.tsv", country_starta_sampling, no_comparisons = 15)
