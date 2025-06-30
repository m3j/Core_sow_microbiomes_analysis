library(patchwork)
source("code/Functions/PCoA_plot.R")
source("code/Functions/dist_utils.R")

# get phyloseq obj
load("data/physeq_objects_filt0.1.RData")

# get distance matric with robust aitchison
robust_ait_full <- robust_ait(ps)

# create pcoa plot for parity
parity_pcoa <- PCoA_plot_generation(
    physeq_obj = ps,
    dist = robust_ait_full,
    fill_var = "parity_index",
    shape_vector = rep(c(21), 5),
    color_vector = universal_colors_10[c(3, 4, 6, 9, 10)],
    label_vector = c("1" = "1, n=117", "2" = "2, n=66", "3" = "3, n=69", "4" = "4, n=61", ">5" = ">5, n=85"),
    plt_title = "",
    legend_title = "Parity",
    legend_pos = c(0.90, 0.10)
)

# add tag and title
parity_pcoa <- parity_pcoa +
    labs(tag = "A) Parity") +
    theme(plot.tag.position = c(0.55, 0.98))

# create pcoa plot for probiotics
probiotics_pcoa <- PCoA_plot_generation(
    physeq_obj = ps,
    dist = robust_ait_full,
    fill_var = "probiotic_any_binary",
    shape_vector = rep(c(21), 2),
    color_vector = probiotic_col,
    label_vector = c("0" = "No, n=229", "1" = "Yes, n=148"),
    plt_title = "",
    legend_title = "Probiotic usage",
    legend_pos = c(0.89, 0.05)
)

# add tag and title
probiotics_pcoa <- probiotics_pcoa +
    labs(tag = "B) Probiotic Usage") +
    theme(plot.tag.position = c(0.55, 0.98))

# combine into one plot
combined <- (parity_pcoa + probiotics_pcoa)

# save plot
ggsave(
    filename = "Figures/Figure_S8.png",
    plot(combined),
    height = 300,
    width = 500,
    units = "mm"
)
