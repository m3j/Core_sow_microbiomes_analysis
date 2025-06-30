rm(list = ls())

library(patchwork)
library(tidyverse)
source("code/Functions/box_plot_utils.R")
source("code/Functions/plot_utils.R")

# get phyloseq object
load("data/physeq_objects_filt0.1.RData")

# get metadata for gestation
metadata <- data.frame(sample_data(ps_gestation))

# create boxplot of performance data
performance_distribution_boxplot <- box_plotting(
    metadata = metadata %>% filter(!is.na("no_piglets_weaned_per_year_sow_not_for_gilts")),
    x = factor(0),
    y = "no_piglets_weaned_per_year_sow_not_for_gilts",
    fill = "no_piglets_weaned_per_year_sow_not_for_gilts_index",
    x_label = "",
    y_label = "Number of weaned piglets/sow/year",
    legend_pos = "right",
    col_vector = prod_col_3,
    label_vector = c("<28.4" = "Bad, n=21", ">37.1" = "Good, n=18", "28.4 <= x <= 37.1" = "Medium, n=76"),
    legend_title = "Number of piglets/sow/year"
)

# save plot
ggsave(
    filename = "Figures/Figure_S3.png",
    plot(performance_distribution_boxplot),
    height = 300,
    width = 250,
    units = "mm"
)
