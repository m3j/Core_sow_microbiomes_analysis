library(phyloseq)

source("code/Functions/physeq_utils.R")
source("code/Functions/plot_utils.R")
source("code/Functions/differential_abundance.R")

# get phyloseq object
load("data/physeq_objects_filt0.1.RData")

# aglommerate in species level
ps_species <- phyloseq::tax_glom(ps, "Species")

# run differential analysis on sampling stage taking Country into account
sampling_country <- ancom_2groups(
    physeq_obj = ps_species,
    formula = "Country + sampling_stage",
    group = "sampling_stage",
    positive_direction = "DA in lactation",
    negative_direction = "DA in gestation",
    rank = "Species"
)

# plot relative abunance of DA species
DA_rabu_sampling_stage <- rabu_with_DA(
    DA_res = sampling_country,
    lfc_col = "lfc_sampling_stageL21",
    predictor = "sampling_stage",
    val1 = "L21",
    val = "G109",
    filt_threshold = 0.005,
    colors = stage_col,
    physeq_obj = ps_species,
    legend_title = "Sampling stage"
)

# save plot
ggsave(
    file = "Figures/figure_S10.png",
    plot = DA_rabu_sampling_stage,
    height = 600,
    width = 700,
    units = "mm"
)
