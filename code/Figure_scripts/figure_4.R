library(patchwork)
source("code/Functions/plot_utils.R")
source("code/Functions/dist_utils.R")
source("code/Functions/PCoA_plot.R")

# get phyloseq object
load("data/physeq_objects_filt0.1.RData")


#  get distance matrix for gestation and lactation in robust aitchison
ra_lac <- robust_ait(ps_lactation)
ra_ges <- robust_ait(ps_gestation)

# plot  pcoa for performance in gestation
weaned_per_year_ges_pcoa <- PCoA_plot_generation(
    physeq_obj = ps_gestation,
    dist = ra_ges,
    fill_var = "no_piglets_weaned_per_year_sow_not_for_gilts_index",
    shape_vector = rep(c(21), 3),
    color_vector = prod_col_3,
    label_vector = c("<28.4" = "Bad, n=21", ">37.1" = "Good, n=18", "28.4 <= x <= 37.1" = "Medium, n=76"),
    plt_title = "",
    legend_title = "No. weaned piglets/sow/year",
    legend_pos = c(0.83, 0.08)
)

# add tag and postition
weaned_per_year_ges_pcoa <- weaned_per_year_ges_pcoa + labs(tag = "A) Gestation") + theme(plot.tag.position = c(0.55, 0.98))

# plot  pcoa for performance in lactation
weaned_per_year_lac_pcoa <- PCoA_plot_generation(
    physeq_obj = ps_lactation,
    dist = ra_lac,
    fill_var = "no_piglets_weaned_per_year_sow_not_for_gilts_index",
    shape_vector = rep(c(21), 3),
    color_vector = prod_col_3,
    label_vector = c("<28.4" = "Bad, n=18", ">37.1" = "Good, n=17", "28.4 <= x <= 37.1" = "Medium, n=76"),
    plt_title = "",
    legend_title = "No. weaned piglets/sow/year",
    legend_pos = c(0.83, 0.08)
)

# add tag and postition
weaned_per_year_lac_pcoa <- weaned_per_year_lac_pcoa + labs(tag = "B) Lactation") + theme(plot.tag.position = c(0.55, 0.98))

# combine into one plot
performance_combined <- weaned_per_year_ges_pcoa + weaned_per_year_lac_pcoa + plot_annotation(
    title = "Perfomance's link to microbial composition",
    theme = theme(
        plot.title = element_text(size = 34)
    )
)

# save plot
ggsave(
    filename = "Figures/Figure4.png",
    plot(performance_combined),
    height = 300,
    width = 550,
    units = "mm"
)
