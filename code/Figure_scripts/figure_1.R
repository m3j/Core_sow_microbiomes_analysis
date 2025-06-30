library(patchwork)
source("code/Functions/PCoA_plot.R")
source("code/Functions/dist_utils.R")

# get phyloseq obj
load("data/physeq_objects_filt0.1.RData")

# get distance matric with robust aitchison
robust_ait_full <- robust_ait(ps)

# create pcoa plot for country
country_pcoa <- PCoA_plot_generation(
    physeq_obj = ps,
    dist = robust_ait_full,
    fill_var = "Country",
    shape_vector = rep(c(21), 6),
    color_vector = country_col,
    label_vector = c(Denmark = "Denmark, n=46, farms=2", Spain = "Spain, n=93,
                                    farms=4", Netherlands = "Netherlands, n=96, farms=4", France = "France, n=58,
                                    farms=2", UK = "UK, n=84, farms=3", Germany = "Germany, n=24, farms=1"),
    plt_title = "",
    legend_title = "Country",
    legend_pos = c(0.87, 0.11)
)

# add tag and title
country_pcoa <- country_pcoa + labs(tag = "A) Country") + theme(plot.tag.position = c(0.55, 0.98))

# create pcoa plot for sampling stage
sampling_pcoa <- PCoA_plot_generation(
    physeq_obj = ps,
    dist = robust_ait_full,
    fill_var = "sampling_stage",
    shape_vector = rep(c(21), 2),
    color_vector = stage_col,
    label_vector = c(L21 = "Lactation day 21, n=197", G109 = "Gestation day 109, n=204"),
    plt_title = "",
    legend_title = "Stage in reproductive cycle",
    legend_pos = c(0.87, 0.05)
)

# add tag and title
sampling_pcoa <- sampling_pcoa + labs(tag = "B) Stage in Reproductive Cycle") + theme(plot.tag.position = c(0.55, 0.98))

# combine into one plot
combined <- (country_pcoa + sampling_pcoa)

# save plot
ggsave(
    filename = "Figures/Figure1.png",
    plot(combined),
    height = 300,
    width = 550,
    units = "mm"
)
