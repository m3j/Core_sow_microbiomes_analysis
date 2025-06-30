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

stage_cols <- stage_col
names(stage_cols) <- unique(metadata$sampling_stage)

# create named list with colours for annotation bars
col_list <- list(farm_number = farm_cols, Country = country_cols, sampling_stage = stage_cols)

# create heatmap with built in save function
abundance_heatmap_generation(
    physeq_obj = ps,
    tax_level = "Species",
    top_taxa = 20,
    annotation_vec = c("Country", "farm_number", "sampling_stage"),
    file_name = "figure_S6.png",
    height = 250,
    width = 900,
    plot_title = "Top 20 abundant species",
    annotation_color_list = col_list
)
