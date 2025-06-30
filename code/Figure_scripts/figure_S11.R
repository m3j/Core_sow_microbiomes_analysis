library(phyloseq)

source("code/Functions/physeq_utils.R")
source("code/Functions/differential_abundance.R")

# get phyloseq object
load("data/physeq_objects_filt0.1.RData")

# aglommerate to genus level
ps_genus <- phyloseq::tax_glom(ps, "Genus")

# differential abundance analysis on sampling stage taking country into account
sampling_country_genus <- ancom_2groups(
    physeq_obj = ps_genus,
    formula = "Country + sampling_stage",
    group = "sampling_stage",
    positive_direction = "DA in lactation",
    negative_direction = "DA in gestation",
    rank = "Genus"
)

# create relative abundance plot
DA_rabu_sampling_stage <- rabu_with_DA(
    DA_res = sampling_country_genus,
    lfc_col = "lfc_sampling_stageL21",
    predictor = "sampling_stage",
    val1 = "L21",
    val = "G109",
    filt_threshold = 0.005,
    colors = stage_col,
    physeq_obj = ps_genus,
    legend_title = "Sampling stage",
    type = "genus"
)

# save plot
ggsave(
    file = "Figures/figure_S11.png",
    plot = DA_rabu_sampling_stage,
    height = 500,
    width = 600,
    units = "mm"
)
