library(tidyverse)
library(phyloseq)
library(patchwork)
source("code/Functions/plot_utils.R")

make_density_plot <- function(data, x_val, x_lab) {
    ggplot(data, aes(x = {{ x_val }})) +
        stat_density(fill = prod_col_2[2], alpha = 0.5) +
        facet_wrap(~farm_number, ncol = 1) +
        labs(
            x = x_lab,
            y = "Density"
        ) +
        theme_minimal() +
        theme(strip.text = element_blank(), panel.grid.minor = element_blank()) +
        scale_y_continuous(position = "right")
}

make_bar_plot <- function(data, x_lab) {
    ggplot(data, aes(x = as.factor(parity))) +
        geom_bar(fill = universal_colors_10[3], alpha = 0.4) +
        facet_wrap(~farm_number, ncol = 1) +
        labs(
            x = x_lab,
            y = "No of Sows"
        ) +
        theme_classic() +
        theme(strip.text = element_blank()) +
        scale_x_discrete(position = "top") +
        guides(fill = FALSE)
}

# get phyloseq object
load("data/physeq_objects_filt0.1.RData")

# get metadata from lactation and gestation samples
metadata_lac <- data.frame(sample_data(ps_lactation))
metadata_ges <- data.frame(sample_data(ps_gestation))

# get number of sows per farm - gestation
ges_farm_count <- metadata_ges %>%
    group_by(farm_number) %>%
    count(sampling_stage) %>%
    select(farm_number, n) %>%
    rename("G109" = "n")

# get feed information - gestation
ges_feed <- metadata_ges %>%
    select(farm_number, GD019_feed_formulation_major_ingredients) %>%
    unique()

# combine sow per farm and feed information - gestation
ges_farm_feed <- ges_farm_count %>% left_join(ges_feed, by = "farm_number")

# get number of sows per farm - lactation
lac_farm_count <- metadata_lac %>%
    group_by(farm_number) %>%
    count(sampling_stage) %>%
    select(farm_number, n) %>%
    rename("L21" = "n")

# get feed information - lactation
lac_feed <- metadata_lac %>%
    select(farm_number, LD21_feed_formulation_major_ingredients) %>%
    unique()

# combine sow per farm and feed information - lactation
lac_farm_feed <- lac_farm_count %>%
    left_join(lac_feed, by = "farm_number") %>%
    data.frame()

# table with farm id and number of sows per timepoint
farmwise_info <- cbind(ges_farm_feed, lac_farm_feed %>% select("L21", "LD21_feed_formulation_major_ingredients")) %>% select(1, 2, 4, 3, 5)

# save table
write_tsv(farmwise_info, "Figures/figure_S1_table.tsv")

# get performance data create density plot
performance_data_per_farm <- metadata_ges %>%
    select(farm_number, no_piglets_weaned_per_year_sow_not_for_gilts)

# create density plot
farm_performance_density <- make_density_plot(data = performance_data_per_farm, x_val = no_piglets_weaned_per_year_sow_not_for_gilts, x_lab = "No piglets weaned/sow/year")

# get parity data
parity_ges <- metadata_ges %>%
    select(farm_number, parity)

# create bar plot
farm_parity <- make_bar_plot(data = parity_ges, x_lab = "Parity")

# combine plots into one figure
combined_fig <- farm_parity | farm_performance_density

# save fig
ggsave(combined_fig,
    file = "Figures/figure_S1.png",
    height = 400,
    width = 250,
    units = "mm"
)
