# Packages needed
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

devtools::install_github("jbisanz/qiime2R")
devtools::install_github("vegandevs/vegan")
devtools::install_github("FrederickHuangLin/ANCOMBC", ref = "bugfix")
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
devtools::install_github("david-barnett/microViz@0.12.0")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    devtools::install_version("BiocManager", version = "1.30.26")
}
BiocManager::install("Maaslin2")

devtools::install_version("ape", version = "5.9-1")
devtools::install_version("compositions", version = "2.0-8")
devtools::install_version("dplyr", version = "1.1.4")
devtools::install_version("ggdist", version = "3.3.3")
devtools::install_version("ggplot2", version = "3.5.2")
devtools::install_version("ggpubr", version = "0.6.0")
devtools::install_version("ggVennDiagram", version = "1.5.4")
devtools::install_version("hash", version = "2.2.6.3")
devtools::install_version("lmerTest", version = "3.1-3")
devtools::install_version("metagenomeSeq", version = "1.46.0")
devtools::install_version("microViz", version = "0.12.4")
devtools::install_version("otuSummary", version = "0.1.2")
devtools::install_version("pairwiseAdonis", version = "0.4.1")
devtools::install_version("patchwork", version = "1.3.0")
devtools::install_version("pheatmap", version = "1.0.13")
devtools::install_version("phyloseq", version = "1.48.0")
devtools::install_version("qiime2R", version = "0.99.6")
devtools::install_version("RColorBrewer", version = "1.1-3")
devtools::install_version("readr", version = "2.1.5")
devtools::install_version("scales", version = "1.4.0")
devtools::install_version("stringr", version = "1.5.1")
devtools::install_version("tibble", version = "3.3.0")
devtools::install_version("tidyr", version = "1.3.1")
devtools::install_version("tidyverse", version = "2.0.0")
devtools::install_version("vegan", version = "2.6-10")
devtools::install_version("VennDiagram", version = "1.7.3")
