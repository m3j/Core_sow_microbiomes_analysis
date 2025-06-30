library(vegan)
library(phyloseq)
library(microViz)
library(pairwiseAdonis)
library(readr)

source("code/Functions/dist_utils.R")

# calculate permaniva and return results as a string
perm_anova_res_string <- function(physeq_obj, fill_var, dist) {
    # get metadata
    metadata <- data.frame(sample_data(physeq_obj))

    # create formula
    formula <- reformulate(termlabels = as.vector(fill_var), response = "dist")

    # run permanova
    permanova <- adonis2(formula,
        data = metadata,
        permutations = 999,
        na.action = na.omit
    )
    # extract the resutls and create a string
    perm_res <- paste0(
        "R2 = ", round(permanova[1, 3], 4),
        "\nF = ", round(permanova[1, 4], 4),
        "\np = ", permanova[1, 5]
    )

    return(perm_res)
}

# calculate permanove and return result object
get_permanova <- function(physeq_obj, fill_var, dist) {
    metadata <- data.frame(sample_data(physeq_obj))
    formula <- reformulate(termlabels = as.vector(fill_var), response = "dist")
    permanova <- adonis2(formula,
        data = metadata,
        permutations = 999,
        na.action = na.omit
    )

    return(permanova)
}

# calculate pairwise anova of robust aitchison distances and return result object
calculate_pairwise_anova <- function(feature, physeq_obj, strata = NULL) {
    # get metadata
    metadata <- data.frame(sample_data(physeq_obj))

    # remove na if any in the feature column
    if (any(is.na(metadata[feature]))) {
        metadata_no_na <- metadata[!is.na(metadata[feature]), ]
        physeq_no_na <- microViz::ps_filter(physeq_obj, !is.na(metadata[feature]))
        # get robust aitchison distances
        robust_ait_dist <- robust_ait(physeq_no_na)
    } else {
        # get robust aitchison distances
        robust_ait_dist <- robust_ait(physeq_obj)
        metadata_no_na <- metadata
    }

    # create formula
    formula <- reformulate(termlabels = as.vector(feature), response = "robust_ait_dist")

    # call pairwise adonis without starta
    if (is.null(strata)) {
        result <- pairwise.adonis2(
            formula,
            data = metadata_no_na,
            p.adjust.m = "holm",
            nperm = 9999
        )
        # call pairwose adonis with strata
    } else {
        result <- pairwise.adonis2(
            formula,
            strata = strata,
            data = metadata_no_na,
            p.adjust.m = "holm",
            nperm = 9999
        )
    }

    return(result)
}

# Create and save table of pairwise results which can be a list if results
pairwise_res2tsv <- function(file_title, pairwise_res, no_comparisons) {
    # create empty dataframe
    pairwise_res_list <- data.frame()

    # loop through results
    for (i in seq_along(pairwise_res)) {
        # for the line showing the comparison formula we wish to get the feature only and add it to the df
        if (i == 1 || i == no_comparisons + 2 || i == 2 * no_comparisons + 3) {
            res_feature <- data.frame(pairwise_res[[i]])
            pairwise_res_list <- dplyr::bind_rows(pairwise_res_list, res_feature)
            next
        }
        # get the model and residuals
        res <- pairwise_res[[i]][1:2, ]

        # get the name of compared variables
        name <- names(pairwise_res[i])[[1]]

        # add the name to the result
        res_name <- dplyr::bind_rows(data.frame(name), res)

        # round the different values
        res_name$SumOfSqs <- round(res_name$SumOfSqs, 0)
        res_name$R2 <- round(res_name$R2, 3)
        res_name$F <- round(res_name$F, 3)
        res_name$"Pr(>F)" <- round(res_name$"Pr(>F)", 4)

        # add the res_name df to the combined res_list
        pairwise_res_list <- dplyr::bind_rows(pairwise_res_list, res_name)
    }
    # save the combined results
    write_tsv(pairwise_res_list, paste0("tables/", file_title))
}
