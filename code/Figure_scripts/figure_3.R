library(tidyverse)
source("code/Functions/plot_utils.R")
source("code/Functions/core_microbiome.R")

# get phyloseq object
load("data/physeq_objects_filt0.1.RData")

# define core species in gestation/lactation
core_lactation_species <- define_core_microbiome(ps_lactation, "Species")
core_gestation_species <- define_core_microbiome(ps_gestation, "Species")

# define sampling stage specific species
gestation_only_species <- anti_join(core_gestation_species, core_lactation_species, by = "Species")
lactation_only_species <- anti_join(core_lactation_species, core_gestation_species, by = "Species")

# define shared species
shared <- inner_join(core_gestation_species, core_lactation_species, by = "Species")

# select only relevant columns
shared <- shared %>% select("Family.x", "Genus.x", "Species", "otu_id.x", "prevalence_proportion.x", "prevalence_proportion.y", "abundance_%.x", "abundance_%.y")

# define column names and rename
column_names <- c("Family", "Genus", "Species", "OTU id", "Prevalence G109", "Prevalence L21", "Relative abundance G109", "Relative abundance L21")
colnames(shared) <- column_names

# create tsv files
write_tsv(lactation_only_species, "Figures/table_fig3_lactation.tsv")
write_tsv(gestation_only_species, "Figures/table_fig3_gestation.tsv")
write_tsv(shared, "Figures/table_fig4_shared.tsv")

# get species as a list for each core
l_species <- core_lactation_species %>% pull(Species)
g_species <- core_gestation_species %>% pull(Species)

# plot venn diagram
plot_venn_diagram_2sets(l_species, "Lactation", g_species, "Gestation", out_path = "Figures/figure3.png")
