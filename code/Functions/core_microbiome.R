library(phyloseq)
library(VennDiagram)

# source color parameters
source("code/Functions/physeq_utils.R")

define_core_microbiome <- function(ps_obj, rank) {
    # defined at abundance 0.5% and prevelance 90%

    # agglommerate at desired level
    ps_obj_glom <- tax_glom(ps_obj, rank)

    # get count table
    otu_count <- as.data.frame(otu_table(ps_obj_glom))

    # define total reads and read number at cutoff
    total_reads <- sum(otu_count)
    abundance_cutoff_0.5perc. <- total_reads * 0.5 / 100

    # calculate abundance and filter based on cutoff
    core_abundant_otus <- get_abundance(ps_obj_glom) %>%
        filter(abundance_reads >= abundance_cutoff_0.5perc.)

    # get otu id of abundant taxa
    abundant_otus <- core_abundant_otus %>% pull(otu_id)

    # get the counts for the abundant otus
    abundant_counts <- otu_count[abundant_otus, ]

    # calculate prevalence proportions
    prevalence_proportion <- apply(abundant_counts, 1, function(x) round(mean(x > 0) * 100, 2)) %>%
        sort(decreasing = TRUE)

    # filter on prevelance and join with abundant otus
    core_otus <- as.data.frame(prevalence_proportion) %>%
        filter(prevalence_proportion > 90) %>%
        rownames_to_column("otu_id") %>%
        left_join(core_abundant_otus, by = "otu_id")

    # get taxa information
    tax <- as.data.frame(tax_table(ps_obj_glom)) %>% rownames_to_column("otu_id")

    # combine core otus with taxa info
    core_otus_with_tax <- core_otus %>%
        left_join(tax, by = "otu_id") %>%
        select(c(otu_id, prevalence_proportion, "abundance_%", Family, Genus, Species))

    return(core_otus_with_tax)
}


# plot venn diagram for two sets
plot_venn_diagram_2sets <- function(set1, title1, set2, title2, out_path, color_list = stage_col) {
    x <- setNames(list(set1, set2), c(title1, title2))
    venn.otu <- venn.diagram(
        x = x,
        filename = out_path,
        fill = rev(color_list),
        alpha = 0.5,
        cex = 4,
        cat.cex = 1.5,
        cat.pos = c(0, 180),
        cat.dist = 0.05,
        output = TRUE
    )
}
