library(tidyverse)
library(ggplot2)
library(vegan)

# get phyloseq object
load("data/physeq_objects_filt0.1.RData")

# get otu counts
otu_counts_full <- otu_table(ps)

# transposed and turn into a dataframe
otu_counts_t <- as.data.frame(t(otu_counts_full))

# for each sample define number of species per 5000 reads, and create a dataframe per sample
data_rare_steps <- lapply(rownames(otu_counts_t), function(sample) {
    s <- otu_counts_t[sample, ]
    rarefied <- vegan::rarefy(s, sample = seq(1, sum(s), by = 5000))
    data.frame(Sample = sample, Effort = seq(1, sum(s), by = 5000), Richness = rarefied)
})

# combine each dataframe per sample into one big dataframe
rare_df <- bind_rows(data_rare_steps)

# add a richness column to plot
rare_df_plot <- rare_df %>%
    rowwise() %>%
    mutate(
        Richness = get(paste0("Richness.N", Effort), envir = as.environment(cur_data()))
    ) %>%
    ungroup()

# plot the rarefaction curves
rarefaction_plot <- ggplot(rare_df_plot, aes(x = Effort, y = Richness, color = Sample)) +
    geom_line(size = 1) +
    labs(
        title = "Rarefaction Curves",
        x = "Sampling Effort",
        y = "Species Richness"
    ) +
    theme_minimal() +
    theme(
        legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20)
    )

# save plot
ggsave(rarefaction_plot,
    file = "Figures/figure_S2.png",
    width = 400,
    height = 200,
    units = "mm"
)
