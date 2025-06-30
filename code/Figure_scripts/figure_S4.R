library(tidyverse)
library(phyloseq)
library(vegan)
library(otuSummary)

source("code/Functions/box_plot_utils.R")
source("code/Functions/plot_utils.R")

# get phyloseq object
load("data/physeq_objects_filt0.1.RData")

# get metadata and otu counts from phyloseq obj
metadata <- data.frame(sample_data(ps))
otu_counts <- data.frame(otu_table(ps))

# set seed and rarefy to count of lowest sample
set.seed(238)
min_sample_size <- min(colSums(otu_counts))
otu_counts_rarefied <- vegan::rrarefy(t(otu_counts), sample = min_sample_size)

# calculate alpha diversity measures
alpha_all <- otuSummary::alphaDiversity(otu_counts_rarefied, siteInCol = FALSE)

# add shannon metrics to metadata
metadata$alpha_shannon <- alpha_all$allBio$shannon

# create boxplot compareing shannong index across countries
country_shannon_plt <- box_plotting_compare(
    metadata,
    x = "country_code",
    y = "alpha_shannon",
    fill = "country_code",
    x_label = "Country",
    y_label = "Shannon",
    compare_list = list(
        c("UK", "DE"), c("UK", "DK"), c("UK", "ES"), c("UK", "FR"), c("UK", "NL"), c("NL", "DE"), c("NL", "DK"),
        c("NL", "ES"), c("NL", "FR"), c("DE", "FR"), c("DK", "FR"), c("ES", "FR"), c("ES", "DE"), c("ES", "DK"),
        c("DE", "DK")
    ),
    color_vector = c("universal_colors_10[c(1, 2, 4, 6, 10, 9)]"),
    legend_pos = "none"
)

# save plot
ggsave(
    file = "Figures/figure_S4.png",
    plot = country_shannon_plt,
    height = 300,
    width = 350,
    units = "mm"
)
