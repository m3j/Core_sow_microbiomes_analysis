library(tidyverse)
library(phyloseq)
library(pheatmap)
library(compositions)

source("code/Functions/plot_utils.R")
source("code/Functions/physeq_utils.R")


abundance_heatmap_generation <- function(
    physeq_obj, tax_level = NA, top_taxa,
    annotation_vec, file_name, width, height,
    plot_title, annotation_color_list) {
    # get metadata from phyloseq object
    metadata <- data.frame(sample_data(physeq_obj))

    # transform to get relative abundance
    relative_physeq_obj <- phyloseq::transform_sample_counts(physeq_obj, function(x) x / sum(x))

    if (!is.na(tax_level)) {
        # if tax_level is provided agglomorate phyloseq object to that level
        relative_physeq_obj_glom <- phyloseq::tax_glom(relative_physeq_obj, taxrank = tax_level)
    }

    # get otu id of the most abundant species
    top <- get_abundance(relative_physeq_obj_glom)[1:top_taxa, ] %>% pull(otu_id)

    # create a phyloseq object of only the most abundant species
    abundance_top <- phyloseq::prune_taxa(top, relative_physeq_obj_glom)

    # if tax_level is not defined set to Species
    if (is.na(tax_level)) {
        tax_level <- "Species"
    }

    # set rownames to taxon based on tax_level
    row_names <- as.data.frame(tax_table(abundance_top)) %>%
        dplyr::select(all_of(tax_level))

    # centered log ratio transformation of data
    clr_otu <- compositions::clr(otu_table(abundance_top))

    # set rownames to that of the tax_level
    row.names(clr_otu) <- row_names[[tax_level]]

    # make names italic
    italic_names <- lapply(
        rownames(clr_otu),
        function(x) bquote(italic(.(x)))
    )

    #  Call heatmap function
    abundance_heatmap <- pheatmap(
        mat = clr_otu,
        show_colnames = TRUE,
        show_rownames = TRUE,
        border_color = NA,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        cutree_cols = 4,
        cutree_rows = 5,
        clustering_method = "complete",
        annotation_col = metadata %>% dplyr::select(all_of(annotation_vec)),
        fontsize = 10,
        fontsize_col = 7,
        cellwidth = 5,
        cellheight = 12,
        main = plot_title,
        annotation_colors = annotation_color_list,
        labels_row = as.expression(italic_names),
        fontfamily = "helvetica"
    )

    # save figure
    ggsave(
        file = paste0("Figures/", file_name),
        plot = abundance_heatmap,
        height = height,
        width = width,
        units = "mm"
    )

    return(abundance_heatmap)
}
