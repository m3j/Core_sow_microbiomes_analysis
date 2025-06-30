source("code/Functions/stats_utils.R")

# get phyloseq object
load("data/physeq_objects_filt0.1.RData")

# calculate pairwise anove in the three performance metrics
performance_starta_sampling <- calculate_pairwise_anova("no_piglets_weaned_per_year_sow_not_for_gilts_index", ps, strata = "sampling_stage")
alive_starta_sampling <- calculate_pairwise_anova("no_piglets_alive_index", ps, strata = "sampling_stage")
stillborn_starta_sampling <- calculate_pairwise_anova("no_stillborn_piglets_index", ps, strata = "sampling_stage")

# combine results and save
pairwise_res2tsv("tabel_S4.tsv", c(performance_starta_sampling, stillborn_starta_sampling, alive_starta_sampling), no_comparisons = 3)
