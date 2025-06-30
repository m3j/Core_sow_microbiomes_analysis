library(patchwork)
library(vegan)
library(otuSummary)
source("code/Functions/box_plot_utils.R")

# get phyloseq object
load("data/physeq_objects_filt0.1.RData")

# get sample data and otu counts from phyloseq obj
metadata <- data.frame(sample_data(ps))
otu_counts <- data.frame(otu_table(ps))

# set seed and rarefy to the size of smallest sample
set.seed(238)
min_sample_size <- min(colSums(otu_counts))
otu_counts_rarefied <- vegan::rrarefy(t(otu_counts), sample = min_sample_size)

# calculate alpha diversity measures
alpha_all <- otuSummary::alphaDiversity(otu_counts_rarefied, siteInCol = FALSE)

# add observed and shannon diversity measures to the metadata
metadata$alpha_observed <- alpha_all$allBio$observed
metadata$alpha_shannon <- alpha_all$allBio$shannon

# plot shannon data
shannon_plt <- box_plotting_no_fill(
    metadata = metadata,
    x = "sampling_stage",
    y = "alpha_shannon",
    x_label = "",
    y_label = "Shannon's diversity index"
)

# plot observed
observed_plt <- box_plotting_no_fill(
    metadata = metadata,
    x = "sampling_stage",
    y = "alpha_observed",
    x_label = "",
    y_label = "Observed OTUs"
)

# combine plots
alpha_plots <- shannon_plt + observed_plt + plot_annotation(
    tag_levels = "A",
    theme = theme_alpha
)

# save plot
ggsave(
    file = "Figures/figure2.png",
    plot = alpha_plots,
    height = 300,
    width = 200,
    units = "mm"
)
