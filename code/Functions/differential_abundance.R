library(ANCOMBC)
library(tidyverse)
library(phyloseq)
library(ggplot2)

source("code/Functions/rabuplot.R")
source("code/Functions/physeq_utils.R")

# set seed before running ANCOMBC2 for reproducibility
set.seed(123)

# will clean ancombcc results
clean_ancombc_data <- function(group, ancom_res, positive_direction, negative_direction, physeq_obj, rank) {
    # define column names dynamically
    lfc_col <- colnames(dplyr::select(ancom_res, contains(paste0("lfc_", group))))
    diff_col <- colnames(dplyr::select(ancom_res, contains(paste0("diff_", group))))
    ss_col <- colnames(dplyr::select(ancom_res, contains(paste0("passed_ss_", group))))

    # print a check for the direction of the differntially abundance analysis
    print(paste0("The positive direction is: ", lfc_col, ". Is that the same as ", positive_direction, "?"))

    # prepare data frame for rabuplot
    res_df <-
        ancom_res %>%
        dplyr::select(taxon, contains(group)) %>%
        # filter LFC above abs(1)
        filter(abs(.[[lfc_col]]) >= 1) %>%
        # arrange according to lfc values
        arrange(desc(.[[lfc_col]]))

    # factorize taxons (numerical order)
    res_df$taxon_order <- factor(res_df$taxon, levels = res_df$taxon)

    # change from otu no. to taxa names
    res_df <- otu_id2tax_names_col(res_df, "taxon", physeq_obj, rank)

    return(res_df)
}

# run ancombc and clean data
ancom_2groups <- function(physeq_obj, formula, group, positive_direction, negative_direction, rank) {
    # run ancombc
    ancombc_output <- ancombc2(
        data = physeq_obj, fix_formula = formula,
        group = group,
        neg_lb = TRUE,
        struc_zero = TRUE,
        verbose = TRUE,
    )

    # extract the results
    ancom_res <- ancombc_output$res

    # define the column indicating differnetial species
    diff_col <- colnames(dplyr::select(ancom_res, contains(paste0("diff_", group))))

    # handle the case of no or only one differntial abundant species.
    if (!any(ancom_res[[diff_col]] == TRUE)) {
        print("no differential abundant species detected")
        ancom_res <- otu_id2tax_names_col(otu_df = ancom_res, otu_col = "taxon", physeq_obj = physeq_obj, rank = "Species")
        return(NULL)
    }
    if (sum(ancom_res[[diff_col]]) == 1) {
        print("only one differentail abundant species/genus - see TSV result file")
        ancom_res <- otu_id2tax_names_col(otu_df = ancom_res, otu_col = "taxon", physeq_obj = physeq_obj, rank = "Species")
        return(ancom_res)
    }

    # clean ancombc data
    ancom_plot_data <- clean_ancombc_data(group = group, ancom_res = ancom_res, positive_direction = positive_direction, negative_direction = negative_direction, physeq_obj = physeq_obj, rank = rank)

    return(ancom_plot_data)
}


# Filter ancom results based on median abundance thresholds in either group
individual_filt <- function(ancom_res, median_threshold, physeq_obj, physeq_split_var, value1, value2) {
    # get differentially abundant OTU ids
    DA_otus <- ancom_res$taxon

    # split the physeeq object based on the variable differnetiated between.
    physeq_list <- split_physeq_in_two(physeq_obj = physeq_obj, feature = physeq_split_var, value1 = value1, value2 = value2)

    # extract the first physeq obj
    physeq_1 <- physeq_list[[1]]

    # filt on median_threshold
    physeq_1_ra_0.5 <- filt_on_median(physeq_1, median_threshold = median_threshold)

    if (nrow(physeq_1_ra_0.5) > 0) {
        # select otus that pass the threshold and is also differentially abundnant(DA_otus)
        physeq_1_ra_0.5_DA <- physeq_1_ra_0.5[physeq_1_ra_0.5$taxon %in% DA_otus, ]
        # extract the otu ids
        otus_1 <- physeq_1_ra_0.5_DA %>% pull(taxon)
    } else {
        # if no otus pass the filter set to NULL
        otu_1 <- NULL
    }

    # extract the second physeq obj
    physeq_2 <- physeq_list[[2]]

    # filt on median_threshold
    physeq_2_ra_0.5 <- filt_on_median(physeq_2, median_threshold = median_threshold)

    if (nrow(physeq_2_ra_0.5) > 0) {
        # select otus that pass the threshold and is also differentially abundnant(DA_otus)
        physeq_2_ra_0.5_DA <- physeq_2_ra_0.5[physeq_2_ra_0.5$taxon %in% DA_otus, ]
        # extract the otu ids
        otus_2 <- physeq_2_ra_0.5_DA %>% pull(taxon)
    } else {
        # if no otus pass the filter set to NULL
        otu_2 <- NULL
    }

    # combine otus that pass the filter in either physeq object and take unique otus.
    DA_otus <- unique(c(otus_1, otus_2))

    return(DA_otus)
}


# prepare data and plot Rabuplot
rabu_with_DA <- function(
    DA_res, lfc_col, predictor, val1, val2, filt_threshold, colors,
    physeq_obj, legend_title, type = "otu", return_list = FALSE) {
    # handle empty case
    if (is.null(DA_res)) {
        print(paste0("no differnetially abundant otus comparing ", predictor))
        return()
    }

    # define the end of the lfc_col string name
    lfc_end <- substr(lfc_col, nchar(lfc_col) - 3, nchar(lfc_col))

    # filter on differentially abundant species and passed sensetivity test
    DA_filtered <- DA_res %>%
        filter(if_any(matches(paste0("^diff.*", lfc_end, "$")), ~.)) %>%
        filter(if_any(matches(paste0("^passed_ss.*", lfc_end, "$")), ~.))

    # run individual filt on median
    DA_otus_filtered <- individual_filt(
        ancom_res = DA_filtered,
        median_threshold = filt_threshold,
        physeq_obj = physeq_obj,
        physeq_split_var = predictor,
        value1 = val1,
        value2 = val2
    )

    # handle the case when no otus pass the individual filtering
    if (!length(DA_otus_filtered) > 0) {
        print("no otus passed the individual filtering")
        return()
    }

    # subset results to selected columns
    lfc_vals_df <- DA_res %>% select(all_of(c(lfc_col, "taxon", "tax_names")))

    # subset to otus passed filtering
    lfc_vals <- lfc_vals_df[lfc_vals_df$taxon %in% DA_otus_filtered, ]

    # define plot title
    title <- "Differentially abundant species"

    # incase type is different from otu change otu to tax
    if (type != "otu") {
        DA_otus_filtered <- lfc_vals$tax_names
        lfc_vals <- lfc_vals %>%
            mutate(taxon = tax_names)
        # reset title to match genera
        title <- "Differentially abundant genera"
    }

    # if return list is passed as parameter the list will be returned instead of the plot
    if (return_list) {
        return(data.frame(lfc_vals)[, c(1, 3)])
    }

    # call the rabuplot function
    rabu_plot <- rabuplot_modified(
        phylo_ob = physeq_obj,
        predictor = predictor,
        type = type,
        list_taxa = DA_otus_filtered,
        no_other_type = TRUE,
        colors = colors,
        otu_id = FALSE,
        predefined_stats = lfc_vals,
        lfc_col = lfc_col,
        legend_title = legend_title,
        main = title
    )
}
