library(patchwork)
library(tidyverse)
source("code/Functions/PCoA_plot.R")
source("code/Functions/dist_utils.R")
source("code/Functions/box_plot_utils.R")

# get phyloseq object
load("data/physeq_objects_filt0.1.RData")

# get meta data from gestation
metadata <- data.frame(sample_data(ps_gestation))

# get robust aitchison distance matrix from gestation samples
robust_ait_ges <- robust_ait(ps_gestation)

# create pcoa plot on feed form in gestation
feed_form_pcoa <- PCoA_plot_generation(
    physeq_obj = ps_gestation,
    dist = robust_ait_ges,
    fill_var = "LD21_feed_form",
    shape_vector = rep(c(21), 3),
    color_vector = c("#1982c4", "#e9c46a", "#676767"),
    label_vector = c(Liquid = "Liquid, n=36", Mash = "Mash, n=43", Pellets = "Pellets, n =118"),
    plt_title = "",
    legend_title = "Feed Form",
    legend_pos = c(0.87, 0.11)
)

# add tag, position and title
feed_form_pcoa <- feed_form_pcoa +
    labs(tag = "A) Microbial composition") +
    theme(plot.tag.position = c(0.55, 0.98))

# greate boxplot of feed form and perfromance in gestation
productivity_feed_form_ges_plt <- box_plotting_compare(
    metadata %>% filter(!is.na(no_piglets_weaned_per_year_sow_not_for_gilts), !is.na(GD109_feed_form)),
    x = "LD21_feed_form",
    y = "no_piglets_weaned_per_year_sow_not_for_gilts",
    fill = "no_piglets_weaned_per_year_sow_not_for_gilts_index",
    x_label = "Feed form in gestation",
    y_label = "Performance",
    legend_title = "Weaned piglets/sow/year",
    compare_list = list(c("Liquid", "Mash"), c("Liquid", "Pellets"), c("Mash", "Pellets")),
    legend_pos = "right",
    color_vector = c("prod_col_3"),
    label_vector = c("<28.4" = "Bad", "28.4 <= x <= 37.1" = "Medium", ">37.1" = "Good")
)

# add tag, position and title
productivity_feed_form_ges_plt <- productivity_feed_form_ges_plt +
    labs(tag = "B) Performance across feed forms") +
    theme(plot.tag.position = c(0.55, 0.98))

# combine into one plot
figure8 <- feed_form_pcoa + productivity_feed_form_ges_plt

# save plot
ggsave(
    file = "Figures/Figure_13S.png",
    plot = figure8,
    height = 300,
    width = 450,
    units = "mm"
)
