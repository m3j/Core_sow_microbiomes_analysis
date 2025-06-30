library(tidyverse)
library(phyloseq)
library(patchwork)

source("code/Functions/physeq_utils.R")
source("code/Functions/plot_utils.R")
source("code/Functions/differential_abundance.R")
source("code/Functions/PCoA_plot.R")
source("code/Functions/dist_utils.R")


# get phyloseq objetct
load("data/physeq_objects_filt0.1.RData")

# agglommerate to species
ps_species <- phyloseq::tax_glom(ps, "Species")

# split phyloseq to gestation and lactation
physeq_divided <- split_physeq_in_two(ps_species, feature = "sampling_stage", value1 = "L21", value2 = "G109")
physeq_lac <- physeq_divided[[1]]
physeq_ges <- physeq_divided[[2]]

# differentially abundance analysis on performance in lactation
productivity3_lac <- ancom_2groups(
    physeq_obj = physeq_lac,
    formula = "no_piglets_weaned_per_year_sow_not_for_gilts_index2",
    group = "no_piglets_weaned_per_year_sow_not_for_gilts_index2",
    positive_direction = "DA in high",
    negative_direction = "DA in low",
    rank = "Species"
)

# differentially abundance analysis on performance in gestation
productivity3_ges <- ancom_2groups(
    physeq_obj = physeq_ges,
    formula = "no_piglets_weaned_per_year_sow_not_for_gilts_index2",
    group = "no_piglets_weaned_per_year_sow_not_for_gilts_index2",
    positive_direction = "DA in high",
    negative_direction = "DA in low",
    rank = "Species"
)

# relative abundance plot of DA species in lactation
da_lac <- rabu_with_DA(productivity3_lac,
    lfc_col = "lfc_no_piglets_weaned_per_year_sow_not_for_gilts_index2Good>37.1",
    predictor = "no_piglets_weaned_per_year_sow_not_for_gilts_index2",
    val1 = "Good>37.1",
    val2 = "Bad<37.1",
    filt_threshold = 0.0005,
    colors = prod_col_2,
    physeq_obj = physeq_lac,
    legend_title = "No weaned piglets/sow/year"
)

# relative abundance plot of DA species in gestation
da_ges <- rabu_with_DA(productivity3_ges,
    lfc_col = "lfc_no_piglets_weaned_per_year_sow_not_for_gilts_index2Good>37.1",
    predictor = "no_piglets_weaned_per_year_sow_not_for_gilts_index2",
    val1 = "Good>37.1",
    al2 = "Bad<37.1",
    filt_threshold = 0.0005,
    colors = prod_col_2,
    physeq_obj = physeq_ges,
    legend_title = "No weaned piglets/sow/year"
)

# get robust aitchison distance matrix
robust_ait_full <- robust_ait(ps)

# plot pcoa of perfomance of the 4 groups of bad and good performers in gestation and lactation
stage_performance_pcoa <- PCoA_plot_generation(
    physeq_obj = ps,
    dist = robust_ait_full,
    fill_var = "enterotypes4",
    shape_vector = rep(c(21), 4),
    color_vector = stage_performance_col,
    label_vector = c(
        "Gestation_good" = "Gestation Good, n=18", "Gestation_bad" = "Gestation Bad, n=97",
        "Lactation_bad" = "Lactation Bad, n=94", "Lactation_good" = "Lactation Good, n=17"
    ),
    plt_title = "",
    legend_title = "Performance groups",
    legend_pos = c(0.87, 0.11)
)

# add tags and adjust sizes for rabu plots
stage_performance_pcoa <- stage_performance_pcoa + labs(tag = "C")
da_ges <- da_ges + labs(tag = "A") + theme(aspect.ratio = 0.5)
da_lac <- da_lac + labs(tag = "B") + theme(aspect.ratio = 0.5)

# combine into one plot
combined <- da_ges / da_lac / stage_performance_pcoa + plot_layout(guides = "collect")

# save plot
ggsave(
    file = "Figures/figure5.png",
    plot = combined,
    height = 800,
    width = 700,
    units = "mm"
)
