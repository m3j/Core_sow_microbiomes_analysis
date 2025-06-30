library(tidyverse)
library(phyloseq)

source("code/Functions/box_plot_utils.R")
source("code/Functions/dist_utils.R")
source("code/Functions/PCoA_plot.R")

# get phyloseq object
load("data/physeq_objects_filt0.1.RData")

# get robust aitchison distance matrices
ra_lac <- robust_ait(ps_lactation)


# create pcoa plot for temperature in lactation unit
temp_lac_pcoa <- PCoA_plot_generation(
    physeq_obj = ps_lactation,
    dist = ra_lac,
    fill_var = "LD21_temp_lactation_unit_index",
    shape_vector = rep(c(21), 3),
    color_vector = c("#1982c4", "#ffca3a", "#d70e17"),
    label_vector = c(
        "<=22" = "Low, <=22, n=48",
        ">25" = "High, >25 , n=78",
        "22 < x <= 25" = "Medium, 22 < x <= 25, n=66"
    ),
    plt_title = "",
    legend_title = "Temperature in lactation unit",
    legend_pos = c(0.84, 0.06)
)

temp_lac_pcoa <- temp_lac_pcoa +
    labs(tag = "A) Microbial composition") +
    theme(plot.tag.position = c(0.55, 0.98))


# get metadata from phyloseq object
metadata <- data.frame(sample_data(ps))

# create boxplot comparing performance across temperatures in lactation
productivity_temp_lac_plt <- box_plotting_compare(
    metadata %>% filter(!is.na(no_piglets_weaned_per_year_sow_not_for_gilts), sampling_stage == "L21"),
    x = "LD21_temp_lactation_unit_index",
    y = "no_piglets_weaned_per_year_sow_not_for_gilts",
    fill = "no_piglets_weaned_per_year_sow_not_for_gilts_index",
    x_label = "Temperature in lactation unit",
    y_label = "Performance",
    legend_title = "Weaned piglets/sow/year",
    compare_list = list(c("<=22", ">25"), c("<=22", "22 < x <= 25"), c(">25", "22 < x <= 25")),
    legend_pos = "right",
    color_vector = c("prod_col_3"),
    label_vector = c("<28.4" = "Bad", "28.4 <= x <= 37.1" = "Medium", ">37.1" = "Good")
)
productivity_temp_lac_plt <- productivity_temp_lac_plt +
    labs(tag = "B) Performance across temperature") +
    theme(plot.tag.position = c(0.55, 0.98))

temp_combined <- temp_lac_pcoa + productivity_temp_lac_plt

# save plot
ggsave(
    filename = "Figures/figure_S12.png",
    plot(temp_combined),
    height = 300,
    width = 550,
    units = "mm"
)
