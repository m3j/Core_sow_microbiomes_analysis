library(phyloseq)
library(microViz)
library(patchwork)

source("code/Functions/dist_utils.R")
source("code/Functions/PCoA_plot.R")
source("code/Functions/physeq_utils.R")

# get phyloseq object
load("data/physeq_objects_filt0.1.RData")

# get metadata from phyloseq object
metadata <- data.frame(sample_data(ps))

# split phyloseq object on diet_change
split_result_diet <- split_physeq_in_two(ps, feature = "diet_change", value1 = "no_change", value2 = "change")
ps_no_change <- split_result_diet$physeq_1
ps_change <- split_result_diet$physeq_2

# plot pcoa of with diet change with fill on sampling stage
pcoa_change <- PCoA_plot_generation(
    ps_change,
    dist = robust_ait(ps_change),
    fill_var = "sampling_stage",
    shape_vector = c(21, 21),
    color_vector = stage_col,
    label_vector = c(L21 = "Lactation day 21", G109 = "Gestation day 109"),
    plt_title = "",
    legend_title = "Sampling Stage"
)

# add tag and position
pcoa_change <- pcoa_change +
    labs(tag = "A) Sampling stage with major ingredients change") +
    theme(plot.tag.position = c(0.41, 0.99))

# plot pcoa of with no diet change with fill on sampling stage
pcoa_no_change <- PCoA_plot_generation(
    ps_no_change,
    dist = robust_ait(ps_no_change),
    fill_var = "sampling_stage",
    shape_vector = c(21, 21),
    color_vector = stage_col,
    label_vector = c(L21 = "Lactation day 21", G109 = "Gestation day 109"),
    plt_title = "",
    legend_title = "Sampling Stage"
)

# add tag and position
pcoa_no_change <- pcoa_no_change +
    labs(tag = " B) Sampling stage with no major ingredients change") +
    theme(plot.tag.position = c(0.34, 0.99))

# combine into one plot
combined_pc12 <- pcoa_change + pcoa_no_change + plot_annotation(
    # tag_levels = "A",
    theme = theme(
        plot.title = element_text(size = 34)
    )
) + plot_layout(guides = "collect")

# save plot
ggsave(
    filename = "Figures/figure_S7.png",
    plot(combined_pc12),
    height = 300,
    width = 550,
    units = "mm"
)
