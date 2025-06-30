library(microViz)
library(phyloseq)
library(ggVennDiagram)

source("code/Functions/plot_utils.R")

# get phyloseq object
load("data/physeq_objects_filt0.1.RData")

# define species in each country
dk_species <- microViz::ps_filter(ps, Country == "Denmark") %>%
    phyloseq::tax_table(.) %>%
    data.frame(.) %>%
    pull("Species")
ge_species <- microViz::ps_filter(ps, Country == "Germany") %>%
    phyloseq::tax_table(.) %>%
    data.frame(.) %>%
    pull("Species")
fr_species <- microViz::ps_filter(ps, Country == "France") %>%
    phyloseq::tax_table(.) %>%
    data.frame(.) %>%
    pull("Species")
nl_species <- microViz::ps_filter(ps, Country == "Netherlands") %>%
    phyloseq::tax_table(.) %>%
    data.frame(.) %>%
    pull("Species")
uk_species <- microViz::ps_filter(ps, Country == "UK") %>%
    phyloseq::tax_table(.) %>%
    data.frame(.) %>%
    pull("Species")
es_species <- microViz::ps_filter(ps, Country == "Spain") %>%
    phyloseq::tax_table(.) %>%
    data.frame(.) %>%
    pull("Species")

# create species and names list
species_sets <- list(dk_species, ge_species, fr_species, nl_species, uk_species, es_species)
country_names <- c("Denmark", "Germany", "France", "Netherlands", "UK", "Spain")

# create venn diagram
species_overlap_country <- ggVennDiagram(species_sets,
    category.names = country_names,
    label = "count",
    set_color = country_col,
    label_alpha = 0,
    label_color = "#827f7f"
)

# save plot
ggsave(filename = "Figures/figure_S5.png", species_overlap_country)
