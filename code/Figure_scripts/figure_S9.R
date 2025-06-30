source("code/Functions/plot_utils.R")
source("code/Functions/heatmap.R")

# get phyloseq object
load("data/physeq_objects_filt0.1.RData")

# get metadata from phyloseq object
metadata <- data.frame(sample_data(ps))

# set annotation bar colors
country_cols <- country_col
names(country_cols) <- unique(metadata$Country)

farm_cols <- c(
    "#9e0142", "#d70e17", "#ff595e", "#ff924c",
    "#efb366", "#e9c46a", "#c5ca30", "#8ab17d", "#2a9d8f",
    "#287271", "#264653", "#36949d", "#1982c4", "#4267ac", "#565aa0", "#6a4c93"
)
names(farm_cols) <- unique(metadata$farm_number)

# make named list for colors
col_list <- list(farm_number = farm_cols, Country = country_cols)

# create heatmap for top 40 in gestation
species_top40_gestation <- abundance_heatmap_generation(
    physeq_obj = ps_gestation,
    tax_level = "Species",
    top_taxa = 40,
    annotation_vec = c("Country", "farm_number"),
    file_name = "figure_S9A.svg",
    height = 250,
    width = 550,
    plot_title = "Top 40 abundant species in gestation",
    annotation_color_list = col_list
)

# create heatmap of top 40 of lactation
species_top40_lactation <- abundance_heatmap_generation(
    physeq_obj = ps_lactation,
    tax_level = "Species",
    top_taxa = 40,
    annotation_vec = c("Country", "farm_number"),
    file_name = "figure_S9B.svg",
    height = 250,
    width = 550,
    plot_title = "Top 40 abundant species in lactation",
    annotation_color_list = col_list
)
