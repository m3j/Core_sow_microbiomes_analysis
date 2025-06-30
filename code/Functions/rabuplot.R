library(Maaslin2)
library(devtools)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(lmerTest)
library(metagenomeSeq)
library(phyloseq)
library(tidyr)

source("code/Functions/physeq_utils.R")

# Rabuplot modified from https://github.com/jstokholm/rabuplot
# added changes are commented:
# predefined_stats is a df containing columns:
# taxon(or otu_ID if type = "otu"), tax_names, lfc column or p-values based on usecase - can also contain adj_p for adjusted p-values.
# otu_id: bool indicating if the plot should contain the OTU_ids rather than tax names.
rabuplot_modified <- function(
    phylo_ob, predictor = "none", type = "genus", relative_abun = TRUE,
    id = NULL, xlabs = "Relative abundance (%)", ylabs = "Average relative abundance",
    main = "Relative abundance plot", violin = TRUE, violin_scale = "width",
    legend_title = predictor, N_taxa = NULL, By_median = TRUE,
    no_other_type = FALSE, legend_names = NULL, Time = "Time",
    Timepoint = NULL, Strata = NULL, Strata_val = "1", no_legends = FALSE,
    no_names = FALSE, italic_names = TRUE, Only_sig = FALSE,
    log = TRUE, log_max = 100, stat_out = FALSE, p_val = TRUE,
    p_stars = FALSE, stats = "non-parametric", p_adjust = FALSE,
    p_adjust_method = "fdr", p_adjust_full = FALSE, colors = NULL,
    color_by = NULL, order = TRUE, reverse = FALSE, list_taxa = NULL,
    select_taxa = NULL, select_type = "genus", bar_chart = FALSE,
    bar_chart_stacked = FALSE, facet_wrap = NULL, facet_label = NULL,
    facet_n = TRUE, percent = FALSE, order_by = "Time", order_val = NULL,
    text_angle_x = 0, rank = NULL, select_rank = NULL, otu_id = FALSE,
    predefined_stats = NULL, lfc_col = NULL) {
    if (!is.null(list_taxa) & is.null(N_taxa)) {
        N_taxa <- length(list_taxa)
    }
    if (is.null(N_taxa) & is.null(list_taxa)) {
        N_taxa <- 15
    }
    if (!is.null(rank)) {
        type <- rank
    }
    if (!is.null(select_rank)) {
        select_type <- select_rank
    }
    options(dplyr.summarise.inform = FALSE)
    if (bar_chart_stacked == TRUE) {
        bar_chart <- TRUE
        p_val <- FALSE
    }
    if (predictor == "none") {
        sample_data(phylo_ob)$none <- "All samples"
        p_val <- FALSE
        if (bar_chart_stacked == FALSE & is.null(color_by)) {
            no_legends <- TRUE
        }
    }
    phylo_ob <- prune_samples(sample_sums(phylo_ob) > 0, phylo_ob)
    otu_mat <- as(otu_table(phylo_ob), "matrix")
    if (taxa_are_rows(phylo_ob)) {
        otu_mat <- t(otu_mat)
    }
    if (!is.null(facet_wrap)) {
        index <- !is.na(get_variable(phylo_ob, predictor)) &
            !is.na(get_variable(phylo_ob, facet_wrap))
    } else {
        index <- !is.na(get_variable(phylo_ob, predictor))
    }
    if (length(unique(index)) != 1) {
        message(paste(length(which(index == F)), "samples have been removed from full dataset (predictor/facet_wrap NAs)"))
    }
    otu_mat <- otu_mat[index, ]
    otu_mat <- otu_mat[, colSums(otu_mat) > 0]
    OTU_index <- colnames(otu_mat)
    tax <- as(tax_table(phylo_ob), "matrix") %>% data.frame(stringsAsFactors = FALSE)
    tax <- tax[rownames(tax) %in% OTU_index, , drop = FALSE]
    tax[is.na(tax)] <- "unclassified"
    tax[tax == ""] <- "unclassified"
    names(tax) <- tolower(names(tax))
    type <- tolower(type)
    if (!is.null(select_type)) {
        select_type <- tolower(select_type)
    }
    tax$OTU <- rownames(tax)
    samp <- data.frame(sample_data(phylo_ob), stringsAsFactors = TRUE)
    samp <- samp[index, ]
    if (is.null(facet_wrap)) {
        samp$wrap <- ""
    }
    if (!is.null(facet_wrap)) {
        samp$wrap <- samp[, facet_wrap]
    }
    if (!is.null(Timepoint)) {
        index <- rownames(samp[(samp[, Time] == Timepoint), ])
        otu_mat <- otu_mat[rownames(otu_mat) %in% index, ]
        otu_mat <- otu_mat[, colSums(otu_mat) > 0]
        OTU_index <- colnames(otu_mat)
        tax <- tax[rownames(tax) %in% OTU_index, ]
        samp <- samp[rownames(samp) %in% index, ]
    }
    if (type == "otu") {
        tax <- rownames_to_column(tax, "otu")
    }
    list <- as.character(tax[, type])
    unique_tax <- unique(list)
    abund <- as.data.frame(matrix(rep(0, (length(unique_tax) *
        nrow(otu_mat))), ncol = length(unique_tax)))
    row.names(abund) <- row.names(otu_mat)
    names(abund) <- unique_tax
    for (i in names(abund)) {
        if (is.array(otu_mat[, list == i])) {
            abund[, i] <- rowSums(otu_mat[, list == i])
        } else {
            abund[, i] <- otu_mat[, list == i]
        }
    }
    abund_org <- abund
    if (relative_abun == TRUE) {
        abund <- apply(abund, 1, function(x) x / sum(x)) %>%
            t() %>%
            as.data.frame()
    }
    abund_all <- abund
    if (is.null(list_taxa) & !is.null(select_taxa)) {
        list_taxa <- NULL
        for (i in 1:length(select_taxa)) {
            list_taxa <- c(list_taxa, (as.character(unique(tax[grep(select_taxa[[i]],
                tax[, select_type],
                ignore.case = TRUE
            ), type]))))
        }
    }
    if (!is.null(list_taxa)) {
        if (is.null(N_taxa)) {
            N_taxa <- length(list_taxa)
        }
        abund <- abund[, colnames(abund) %in% list_taxa, drop = FALSE]
        unique_tax <- names(abund)
    }
    if (length(abund) > 1) {
        index <- !is.na(rownames(samp))
        if (!is.null(order_val)) {
            index <- samp[, order_by] == order_val
        }
        abund <- abund[, order(-colSums(abund[index, ]))]
        if (By_median) {
            abund <- abund[, order(-apply(
                abund[index, ], 2,
                median
            ))]
        }
        if ("unclassified" %in% unique_tax) {
            abund <- abund[c(
                setdiff(names(abund), "unclassified"),
                "unclassified"
            )]
        }
        if (N_taxa < length(unique_tax)) {
            abund <- abund[-(length(unique_tax) - (length(unique_tax) -
                N_taxa) + 1):-length(unique_tax)]
        }
        if (no_other_type == FALSE) {
            abund[, paste("Other", type)] <- rowSums(abund_all[
                ,
                !names(abund_all) %in% names(abund)
            ])
        }
    }
    # if predifined stats are included
    if (!is.null(predefined_stats)) {
        # make names fit removing "_"
        predefined_stats$tax_names <- gsub("_", " ", predefined_stats$tax_names)
        # create pval tibble with taxon lfc_col renamed to "variable" and 'pval' to fit original pval df
        pval <- as_tibble(predefined_stats) %>%
            select(taxon, lfc_col) %>%
            mutate(wrap = "") %>%
            mutate(p_adjust = "")
        colnames(pval)[1:2] <- c("variable", "pval")
    } else {
        index <- !is.na(rownames(samp))
        if (!is.null(Strata)) {
            index <- samp[, Strata] == Strata_val
        }
        samp2 <- samp %>% filter(index)
        if (p_val == TRUE & (bar_chart == FALSE | (bar_chart == TRUE &
            bar_chart_stacked == FALSE))) {
            if (p_adjust_full == TRUE | stats == "mgs_feature") {
                abund2 <- abund_org %>% filter(index)
                if (relative_abun == TRUE & stats != "mgs_feature") {
                    abund2 <- apply(abund2, 1, function(x) x / sum(x)) %>%
                        t() %>%
                        as.data.frame()
                }
            } else {
                abund2 <- abund %>% filter(index)
            }
            if (stats == "mgs_feature" & length(levels(factor(samp2[
                ,
                predictor
            ]))) > 2) {
                stats <- "non-parametric"
                message("MGS not available for >2 predictors, switching to non-parametric")
            }
            if (stats == "mixed" & is.null(id)) {
                stats <- "non-parametric"
                message("No id variable for mixed model, switching to non-parametric")
            }
            if (stats == "mixed") {
                message("Mixed model statistics")
                pred <- samp2[, predictor]
                id <- samp2[, id]
                if (!is.null(facet_wrap)) {
                    message(paste0("Using ", facet_wrap, " as mixed model group variable"))
                    facet <- samp2[, "wrap"]
                    pval <- cbind(abund2, pred, id, facet) %>%
                        as_tibble() %>%
                        gather(variable, value, -"pred", -"id", -"facet") %>%
                        group_by(variable) %>%
                        summarize(pval = lmerTest::lmer(value ~
                            pred + factor(facet) + (1 | id)) %>% anova() %>%
                            filter(row_number() == 1) %>% .$"Pr(>F)", .groups = "drop")
                } else {
                    message(paste0("No group variable defined for mixed model in facet_wrap"))
                    pval <- cbind(abund2, pred, id) %>%
                        as_tibble() %>%
                        gather(variable, value, -"pred", -"id") %>%
                        group_by(variable) %>%
                        summarize(
                            pval = lmerTest::lmer(value ~
                                pred + (1 | id)) %>% anova() %>% .$"Pr(>F)",
                            .groups = "drop"
                        )
                }
                pval <- pval %>% mutate(wrap = "Mixed", p_adjust = p.adjust(
                    pval,
                    p_adjust_method
                ))
                pval$wrap <- factor(pval$wrap, levels = c(levels(factor(samp2[
                    ,
                    facet_wrap
                ])), "Mixed"))
            } else {
                pval <- data.frame()
                for (i in 1:length(unique(samp2$wrap))) {
                    index <- samp2$wrap == unique(samp2$wrap)[[i]]
                    abund3 <- abund2 %>% filter(index)
                    pred <- samp2[index, predictor]
                    if (stats == "mgs_feature") {
                        mgs <- metagenomeSeq::newMRexperiment(counts = t(abund3))
                        mgsp <- metagenomeSeq::cumNormStat(mgs)
                        mgs <- metagenomeSeq::cumNorm(mgs, mgsp)
                        mod <- model.matrix(~ as.numeric(pred == unique(pred)[1]))
                        if (length(unique(samp2$wrap)) > 1) {
                            message(paste0(
                                "MGS FeatureModel for facet_wrap = ",
                                unique(samp2$wrap)[[i]]
                            ))
                        } else {
                            message("MGS FeatureModel")
                        }
                        mgsfit <- metagenomeSeq::fitFeatureModel(
                            obj = mgs,
                            mod = mod
                        )
                        pval_tmp <- data.frame(
                            variable = mgsfit$taxa,
                            pval = mgsfit$pvalues
                        )
                    }
                    if (stats == "non-parametric") {
                        if (i == 1) {
                            message("Non-parametric statistics")
                        }
                        pval_tmp <- cbind(abund3, pred) %>%
                            as_tibble() %>%
                            gather(variable, value, -"pred") %>%
                            group_by(variable) %>%
                            summarize(
                                pval = kruskal.test(value ~ pred)$p.value,
                                .groups = "drop"
                            )
                    }
                    if (stats == "parametric") {
                        if (i == 1) {
                            message("Parametric statistics")
                        }
                        pval_tmp <- cbind(abund3, pred) %>%
                            as_tibble() %>%
                            gather(variable, value, -"pred") %>%
                            group_by(variable) %>%
                            summarize(
                                pval = oneway.test(value ~ pred)$p.value,
                                .groups = "drop"
                            )
                    }
                    if (stats == "maaslin2") {
                        if (i == 1) {
                            message("Maaslin2 statistics")
                        }
                        test_set <- setNames(abund3, paste0("X", seq_along(abund3)))
                        capture.output(
                            {
                                fit_data <- Maaslin2::Maaslin2(test_set,
                                    data.frame(pred) %>% data.frame(row.names = rownames(test_set)),
                                    output = tempdir(), fixed_effects = "pred",
                                    normalization = "NONE", transform = "LOG",
                                    min_abundance = 0, min_prevalence = 0,
                                    min_variance = 0, plot_scatter = FALSE,
                                    plot_heatmap = FALSE, standardize = FALSE
                                )$results
                            },
                            file = nullfile()
                        )
                        pval_tmp <- data.frame(
                            variable = names(abund3),
                            pval = fit_data[
                                match(names(test_set), fit_data$feature),
                                "pval"
                            ]
                        )
                    }
                    pval_tmpe <<- pval_tmp
                    pval_tmp <- pval_tmp %>% mutate(
                        wrap = unique(samp2$wrap)[[i]],
                        p_adjust = p.adjust(pval, p_adjust_method)
                    )
                    pval <- rbind(pval, pval_tmp)
                }
            }
            if (p_adjust) {
                message(paste(
                    p_adjust_method, "correction applied for",
                    length(unique(pval$variable)), "taxa"
                ))
            }
        }
    }
    # if input type is otu but do not to wish to use otu id we get the corresponding taxa
    if (type == "otu" && otu_id == FALSE) {
        # get tax names using code/physeq_utils::otu2tax()
        col_names_tmp <- sapply(colnames(abund), function(x) otu2tax(physeq_obj = phylo_ob, rank = "species", otu = x))
        col_names <- as.data.frame(col_names_tmp) %>% pull(col_names_tmp)
        # set column names to taxa instead of otu id
        colnames(abund) <- col_names

        # get the tax names corresponding to pval$variable (OTUid)
        pval_variable_names <- sapply(pval$variable, function(x) otu2tax(physeq_obj = phylo_ob, rank = "species", otu = x))
        var_tax <- as.data.frame(pval_variable_names) %>% pull(pval_variable_names)
        # change $variable to tax names(instead of OTU id)
        pval$variable <- var_tax
    }
    bacteria <- rev(names(abund))
    subset <- cbind(samp[!names(samp) %in% bacteria], abund)
    subset$predictor2 <- as.factor(subset[, predictor])
    subset$ID <- rownames(subset)
    if (!is.null(Strata)) {
        subset[, Strata] <- as.factor(subset[, Strata])
    }
    if (!is.null(facet_wrap)) {
        subset$wrap <- as.factor(subset[, facet_wrap])
        if (!is.null(Strata)) {
            molten <- subset[, c(
                "ID", paste(bacteria), "predictor2",
                Strata, "wrap"
            )] %>% gather(
                variable, value,
                -"predictor2", -"ID", -all_of(Strata), -"wrap"
            )
        } else {
            molten <- subset[, c(
                "ID", paste(bacteria), "predictor2",
                "wrap"
            )] %>% gather(
                variable, value, -"predictor2",
                -"ID", -"wrap"
            )
        }
    }
    if (is.null(facet_wrap)) {
        if (!is.null(Strata)) {
            molten <- subset[, c(
                "ID", paste(bacteria), "predictor2",
                Strata
            )] %>% gather(
                variable, value, -"predictor2",
                -"ID", -all_of(Strata)
            )
        } else {
            molten <- subset[, c("ID", paste(bacteria), "predictor2")] %>%
                gather(variable, value, -"predictor2", -"ID")
        }
    }
    if (!is.null(color_by)) {
        molten[molten$variable != paste("Other", type), "colvar"] <- molten %>%
            dplyr::filter(variable != paste("Other", type)) %>%
            .[, "variable"] %>%
            match(tax[, type]) %>%
            tax[
                .,
                color_by
            ] %>%
            as.character()
        molten[molten$variable == paste("Other", type), "colvar"] <- paste(
            "Other",
            color_by
        ) %>% as.character()
    }
    molten$variable <- gsub("_", " ", molten$variable)
    # if predefined stats is provided set levels and factorize tax_names
    if (order && !is.null(predefined_stats)) {
        ordered <- levels(factor(predefined_stats$tax_names, levels = unique(predefined_stats$tax_names)))
    }
    molten <- molten %>% arrange(factor(variable, levels = ordered))
    if (order && is.null(predefined_stats)) {
        ordered <- unique(molten$variable)
    }
    if (!order) {
        ordered <- sort(unique(molten$variable))
    }
    molten$variable <- factor(molten$variable, levels = ordered)
    if (is.null(color_by)) {
        molten$colvar <- molten$variable
    }
    if (!is.null(Strata)) {
        molten <- molten[which(molten[, Strata] == Strata_val), ]
    }
    if (is.null(colors)) {
        cols <- c(
            brewer.pal(8, "Set1"), brewer.pal(7, "Dark2"),
            brewer.pal(7, "Set2"), brewer.pal(12, "Set3"), brewer.pal(
                7,
                "Accent"
            ), brewer.pal(12, "Paired"), "gray"
        )
        cols <- cols[1:length(levels(factor(molten$predictor2)))]
    }
    if (!is.null(colors)) {
        cols <- colors
    }
    if (bar_chart == TRUE & bar_chart_stacked == FALSE & is.null(legend_names)) {
        legend_names <- as.character(levels(factor(molten$predictor2)))
    }
    if (is.null(legend_names)) {
        legend_names <- as.character(levels(factor(molten$predictor2)))
    }
    ordered2 <- rev(unique(molten$colvar))
    if (reverse) {
        if (bar_chart == FALSE) {
            molten$predictor2 <- factor(molten$predictor2, levels = rev(levels(molten$predictor2)))
            legend_names <- rev(legend_names)
            cols <- rev(cols)
        }
        if (bar_chart == TRUE) {
            molten$colvar <- factor(molten$colvar, levels = rev(levels(factor(molten$colvar))))
            molten$variable <- factor(molten$variable, levels = rev(levels(factor(molten$variable))))
            cols <- rev(cols)
            ordered2 <- rev(ordered2)
        }
    }
    if (bar_chart) {
        log <- FALSE
        cols <- c(
            brewer.pal(8, "Set1"), brewer.pal(7, "Dark2"),
            brewer.pal(7, "Set2"), brewer.pal(12, "Set3"), brewer.pal(
                7,
                "Accent"
            ), brewer.pal(12, "Paired"), "gray"
        )
        if (is.null(color_by) & bar_chart_stacked == FALSE) {
            cols <- cols[1:length(levels(factor(molten$predictor2)))]
        } else {
            cols <- cols[c(1:length(levels(factor(molten$colvar))) -
                1, length(cols))]
        }
        if (!is.null(colors)) {
            cols <- colors
        }
        if (is.null(color_by) & reverse == FALSE) {
            cols <- rev(cols)
        }
        if (!is.null(color_by) & reverse == TRUE) {
            cols <- rev(cols)
        }
        if (is.null(facet_wrap)) {
            molten$wrap <- ""
        }
        molten_mean <- molten %>%
            dplyr::group_by(
                variable, predictor2,
                wrap, colvar
            ) %>%
            dplyr::summarize(value = mean(value))
        molten_mean$colvar <- factor(molten_mean$colvar, levels = ordered2)
    }

    if (p_val == TRUE & ((bar_chart == TRUE & bar_chart_stacked ==
        FALSE) | bar_chart == FALSE) & is.null(color_by)) {
        if (is.null(facet_wrap)) {
            molten$wrap <- ""
        }
        if (!is.null(facet_wrap)) {
            pval <- data.frame(pval = pval[gsub("_", " ", pval$variable) %in%
                ordered, ]$pval, p_adjust = pval[gsub(
                "_", " ",
                pval$variable
            ) %in% ordered, ]$p_adjust, variable = gsub(
                "_",
                " ", pval[gsub("_", " ", pval$variable) %in%
                    ordered, ]$variable
            ), wrap = pval[gsub(
                "_",
                " ", pval$variable
            ) %in% ordered, ]$wrap)
        } else {
            pval <- data.frame(pval = pval[gsub("_", " ", pval$variable) %in%
                ordered, ]$pval, p_adjust = pval[gsub(
                "_", " ",
                pval$variable
            ) %in% ordered, ]$p_adjust, variable = gsub(
                "_",
                " ", pval[gsub("_", " ", pval$variable) %in%
                    ordered, ]$variable
            ))
            # if (length(pval$variable) - length(ordered) < 0) {
            #   pval <- pval[match(pval$variable, ordered[length(pval$variable) -
            #     length(ordered)]), ]
            # }
        }

        pval$predictor2 <- molten$predictor2[1]
        pval$pval <- ifelse(is.na(pval$pval), 1, pval$pval)
        pval$p_adjust <- ifelse(is.na(pval$p_adjust), 1, pval$p_adjust)
        if (Only_sig) {
            index <- pval[pval$pval < 0.05, "variable"]
            molten <- molten[molten$variable %in% index, ]
            pval <- pval[pval$pval < 0.05, ]
        }
        if (stat_out) {
            median_iqr <<- molten %>%
                dplyr::group_by(
                    variable,
                    predictor2
                ) %>%
                dplyr::summarize(
                    N = length(value),
                    median = median(value) * 100, Q1 = quantile(
                        value,
                        1 / 4
                    ) * 100, Q3 = quantile(value, 3 / 4) * 100,
                    IQR = IQR(value)
                ) %>%
                as.data.frame()
            pval_out <<- pval
            mean_sd <<- molten %>%
                dplyr::group_by(
                    variable,
                    predictor2
                ) %>%
                dplyr::summarize(
                    N = length(value),
                    mean = mean(value) * 100, sd = sd(value) * 100
                ) %>%
                as.data.frame()
        }
    }
    if (bar_chart == FALSE) {
        if (ncol(tax) >= 6) {
            molten$value <- molten$value + 1e-06
        } else {
            molten$value <- molten$value + 0.001
        }
        ordered <- levels(factor(molten$colvar))
        p <- ggplot(molten, aes(x = variable, y = value, fill = predictor2)) +
            {
                if (violin) {
                    geom_violin(
                        scale = violin_scale, width = 0.65,
                        position = position_dodge(width = 0.9), size = 1,
                        color = "#00000000"
                    )
                } else {
                    geom_boxplot(
                        width = 0.55, position = position_dodge(width = 0.8),
                        size = 0.3, outlier.size = 0, outlier.color = "grey"
                    )
                }
            } +
            {
                if (violin) {
                    stat_summary(
                        fun = median, fun.min = min, fun.max = max,
                        geom = "point", size = 0.8, color = "black",
                        position = position_dodge(width = 0.9)
                    )
                } else {
                    stat_summary(
                        fun = median, fun.min = min, fun.max = max,
                        geom = "point", size = 0.8, color = "#00000000",
                        position = position_dodge(width = 0.9)
                    )
                }
            } +
            theme_bw() +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), legend.key = element_blank(),
                axis.title = element_text(size = 28), legend.text = element_text(size = 26),
                axis.text = element_text(size = 26), strip.text = element_text(size = 26),
                legend.key.size = unit(0.5, "cm"), text = element_text(size = 26)
            ) +
            coord_flip() +
            xlab(NULL) +
            ylab(xlabs) +
            ggtitle(main)
        if (length(unique(molten$variable)) > 1) {
            p <- p + geom_vline(xintercept = seq(1.5, length(unique(molten$variable)) -
                0.5, 1), lwd = 0.2, colour = "grey")
        }
        p <- p + scale_fill_manual(values = cols, labels = legend_names) +
            guides(fill = guide_legend(
                title = legend_title,
                reverse = TRUE, override.aes = list(
                    linetype = 0,
                    shape = 16, color = rev(cols), size = 5, bg = "white"
                )
            ))
    }
    if (bar_chart == TRUE) {
        if (bar_chart_stacked == TRUE) {
            p <- ggplot(molten_mean, aes(
                x = factor(predictor2,
                    levels = legend_names, labels = legend_names
                ),
                y = value, fill = variable
            )) +
                theme_bw() +
                geom_bar(stat = "identity") +
                theme_bw() +
                theme(
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), legend.key = element_blank(),
                    axis.title = element_text(size = 20), legend.text = element_text(size = 16),
                    axis.text = element_text(size = 16), strip.text = element_text(size = 16),
                    legend.key.size = unit(0.5, "cm"), text = element_text(size = 16)
                ) +
                xlab(NULL) +
                ylab(ylabs) +
                ggtitle(main) +
                scale_fill_manual(
                    values = cols,
                    labels = ordered
                ) +
                guides(fill = guide_legend(title = NULL))
        }
        if (bar_chart_stacked == FALSE) {
            if (!is.null(color_by)) {
                p <- ggplot(molten_mean, aes(
                    x = variable, y = value,
                    fill = colvar, group = wrap
                )) +
                    geom_bar(
                        stat = "identity",
                        position = position_dodge(width = 0.95)
                    ) +
                    scale_fill_manual(values = cols, labels = ordered2) +
                    guides(fill = guide_legend(title = color_by))
            } else {
                p <- ggplot(molten_mean, aes(
                    x = variable, y = value,
                    fill = predictor2
                )) +
                    geom_bar(
                        stat = "identity",
                        position = position_dodge(width = 0.95)
                    ) +
                    scale_fill_manual(values = cols, labels = legend_names) +
                    guides(fill = guide_legend(title = legend_title))
            }
            if (length(unique(molten_mean$variable)) > 1) {
                p <- p + geom_vline(xintercept = seq(1.5, length(unique(molten_mean$variable)) -
                    0.5, 1), lwd = 0.2, colour = "grey")
            }
            p <- p + theme_bw() + theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), legend.key = element_blank(),
                axis.title = element_text(size = 20), legend.text = element_text(size = 20),
                axis.text = element_text(size = 20), strip.text = element_text(size = 20),
                legend.key.size = unit(0.5, "cm"), text = element_text(size = 20)
            ) +
                xlab(NULL) + ylab(ylabs) + ggtitle(main) + theme(strip.background = element_blank()) +
                coord_flip()
        }
    }
    if (!is.null(facet_wrap)) {
        if (is.null(facet_label)) {
            label_names <- levels(factor(samp[, facet_wrap]))
        }
        if (!is.null(facet_label)) {
            label_names <- facet_label
        }
        if (facet_n == TRUE) {
            label_names <- samp2 %>%
                dplyr::group_by(get(facet_wrap)) %>%
                dplyr::summarise(n = n()) %>%
                dplyr::mutate(pasted_label = paste0(levels(factor(samp2[
                    ,
                    facet_wrap
                ])), ", n = ", n))
            label_names <- as.character(label_names$pasted_label)
        }
        names(label_names) <- levels(factor(samp2[, facet_wrap]))
        if (stats == "mixed") {
            label_names <- c(label_names, "Mixed")
            names(label_names) <- c(
                levels(factor(samp2[, facet_wrap])),
                "Mixed"
            )
        }
        p <- p + facet_grid(~wrap,
            labeller = labeller(wrap = label_names),
            scales = "free", space = "free"
        ) + theme(strip.background = element_blank())
        if (bar_chart == FALSE) {
            p$layers[4:5] <- NULL
        }
    }
    if (italic_names == TRUE & (bar_chart == FALSE | (bar_chart ==
        TRUE & bar_chart_stacked == FALSE))) {
        if (type == "genus" | type == "family" | type == "species" | type == "otu") {
            if (no_other_type == F) {
                p <- suppressWarnings(p + theme(axis.text.y = element_text(face = c(
                    "plain",
                    rep("italic", length(unique(molten$variable)) -
                        1)
                ))))
            } else {
                p <- p + theme(axis.text.y = element_text(face = "italic"))
            }
        }
    }
    if (!is.null(color_by)) {
        if (color_by == "genus" | color_by == "family" | color_by ==
            "species") {
            p <- p + theme(legend.text = element_text(face = "italic"))
        }
        if (color_by == type & bar_chart_stacked == FALSE) {
            p <- p + theme(legend.position = "none")
        }
    }
    if (p_val == TRUE) {
        if (log == FALSE) {
            if (bar_chart == TRUE) {
                pval$y <- max(molten_mean$value) * 1.1
            } else {
                pval$y <- max(molten$value) * 1.15
            }
        } else {
            pval$y <- ifelse(log_max == 100, 10, ifelse(log_max ==
                10, 0.126, 0.0126))
        }
        if (p_adjust == TRUE) {
            if (log == FALSE & bar_chart == FALSE) {
                pval$y_adjust <- 1.22
            }
            if (log == FALSE & bar_chart == TRUE) {
                pval$y_adjust <- max(molten_mean$value) * 1.25
            }
            if (log == TRUE) {
                pval$y_adjust <- ifelse(log_max == 100, 105,
                    ifelse(log_max == 10, 1.26, 0.126)
                )
            }
        }
    }
    if (log == TRUE) {
        if (p_val == FALSE) {
            if (log_max == 100) {
                p <- p + scale_y_log10(breaks = c(
                    1e-06, 0.001,
                    0.01, 0.1, 1
                ), labels = c(
                    "0%", "0.1%", "1%",
                    "10%", "100%"
                ))
            }
            if (log_max == 10) {
                p <- p + scale_y_log10(
                    limits = c(0.001, 0.13),
                    breaks = c(0.001, 0.01, 0.05, 0.1), labels = c(
                        "0%",
                        "1%", "5%", "10%"
                    )
                )
            }
            if (log_max == 1) {
                p <- p + scale_y_log10(
                    limits = c(0.001, 0.013),
                    breaks = c(0.001, 0.01), labels = c("0%", "1%")
                )
            }
        }
        if (p_val == TRUE) {
            if (p_adjust) {
                if (log_max == 100) {
                    p <- p + scale_y_log10(breaks = c(
                        1e-06, 0.001,
                        0.01, 0.1, 1, 7, 70
                    ), labels = c(
                        "0%", "0.1%",
                        "1%", "10%", "100%", "P-value", "q-value"
                    ))
                }
                if (log_max == 10) {
                    p <- p + scale_y_log10(breaks = c(
                        0.001, 0.01,
                        0.05, 0.1, 0.126, 1.26
                    ), labels = c(
                        "0%",
                        "1%", "5%", "10%", "P-value", "q-value"
                    ))
                }
                if (log_max == 1) {
                    p <- p + scale_y_log10(breaks = c(
                        0.001, 0.01,
                        0.0126, 0.126
                    ), labels = c(
                        "0%", "1%", "P-value",
                        "q-value"
                    ))
                }
            } else {
                if (log_max == 100) {
                    p <- p + scale_y_log10(breaks = c(
                        1e-06, 0.001,
                        0.01, 0.1, 1, 7
                    ), labels = c(
                        "0%", "0.1%",
                        "1%", "10%", "100%", "LFC-value"
                    ))
                }
                if (log_max == 10) {
                    p <- p + scale_y_log10(breaks = c(
                        0.001, 0.01,
                        0.05, 0.1, 0.126
                    ), labels = c(
                        "0%", "1%",
                        "5%", "10%", "LCF-value"
                    ))
                }
                if (log_max == 1) {
                    p <- p + scale_y_log10(breaks = c(
                        0.001, 0.01,
                        0.0126
                    ), labels = c("0%", "1%", "LFC-value"))
                }
            }
        }
    }
    if (log == FALSE) {
        if (p_val == FALSE) {
            p <- p + scale_y_continuous(breaks = c(
                0, 0.25, 0.5,
                0.75, 1
            ), labels = c(
                "0%", "25%", "50%", "75%",
                "100%"
            ))
        }
        if (p_val == TRUE) {
            if (p_adjust == TRUE) {
                if (max(molten_mean$value) <= 1 & max(molten_mean$value) >=
                    0.75) {
                    p <- p + scale_y_continuous(breaks = c(
                        0, 0.25,
                        0.5, 0.75, 1, max(molten_mean$value) * 1.07,
                        max(molten_mean$value) * 1.25
                    ), labels = c(
                        "0%",
                        "25%", "50%", "75%", "100%", "P-value", "q-value"
                    ))
                }
                if (max(molten_mean$value) < 0.75 & max(molten_mean$value) >=
                    0.5) {
                    p <- p + scale_y_continuous(breaks = c(
                        0, 0.25,
                        0.5, max(molten_mean$value) * 1.07, max(molten_mean$value) *
                            1.25
                    ), labels = c(
                        "0%", "25%", "50%", "P-value",
                        "q-value"
                    ))
                }
                if (max(molten_mean$value) < 0.5 & max(molten_mean$value) >=
                    0.25) {
                    p <- p + scale_y_continuous(breaks = c(
                        0, 0.1,
                        0.2, 0.3, 0.4, max(molten_mean$value) * 1.07,
                        max(molten_mean$value) * 1.25
                    ), labels = c(
                        "0%",
                        "10%", "20%", "30%", "40%", "P-value", "q-value"
                    ))
                }
                if (max(molten_mean$value) < 0.25) {
                    p <- p + scale_y_continuous(breaks = c(
                        0, 0.05,
                        0.1, 0.15, 0.2, max(molten_mean$value) *
                            1.07, max(molten_mean$value) * 1.25
                    ), labels = c(
                        "0%",
                        "5%", "10%", "15%", "20%", "P-value", "q-value"
                    ))
                }
            }
            if (p_adjust == FALSE) {
                if (max(molten_mean$value) <= 1 & max(molten_mean$value) >=
                    0.75) {
                    p <- p + scale_y_continuous(
                        breaks = c(
                            0, 0.25,
                            0.5, 0.75, 1, max(molten_mean$value) * 1.07
                        ),
                        labels = c(
                            "0%", "25%", "50%", "75%", "100%",
                            "P-value"
                        )
                    )
                }
                if (max(molten_mean$value) < 0.75 & max(molten_mean$value) >=
                    0.5) {
                    p <- p + scale_y_continuous(breaks = c(
                        0, 0.25,
                        0.5, max(molten_mean$value) * 1.07
                    ), labels = c(
                        "0%",
                        "25%", "50%", "P-value"
                    ))
                }
                if (max(molten_mean$value) < 0.5 & max(molten_mean$value) >=
                    0.25) {
                    p <- p + scale_y_continuous(
                        breaks = c(
                            0, 0.1,
                            0.2, 0.3, 0.4, max(molten_mean$value) * 1.07
                        ),
                        labels = c(
                            "0%", "10%", "20%", "30%", "40%",
                            "P-value"
                        )
                    )
                }
                if (max(molten_mean$value) < 0.25) {
                    p <- p + scale_y_continuous(breaks = c(
                        0, 0.05,
                        0.1, 0.15, 0.2, max(molten_mean$value) *
                            1.07
                    ), labels = c(
                        "0%", "5%", "10%", "15%",
                        "20%", "P-value"
                    ))
                }
            }
        }
    }
    p <- p + theme(
        plot.background = element_blank(), panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)
    )
    p <- p + theme(axis.text.x = element_text(
        angle = text_angle_x,
        vjust = ifelse(text_angle_x < 0 & text_angle_x > -90,
            0, ifelse(text_angle_x >= 0 & text_angle_x < 90,
                1, 0.5
            )
        ), hjust = ifelse(text_angle_x == 0 |
            text_angle_x == 180 | text_angle_x == -180 | text_angle_x ==
            180, 0.5, ifelse((text_angle_x < 0 & text_angle_x >=
            -90) | text_angle_x >= 270, 0, 1))
    ))
    if (bar_chart == TRUE & bar_chart_stacked == FALSE & percent ==
        TRUE) {
        p <- p + geom_text(aes(label = paste0(sprintf(
            "%.2f",
            value * 100
        ), "%")), hjust = -0.12, position = position_dodge(width = 0.95)) +
            scale_y_continuous(limits = c(0, max(molten_mean$value) +
                0.2), labels = scales::percent)
    }
    if (no_legends) {
        p <- p + theme(legend.position = "none")
    }
    if (no_names) {
        p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }
    stars.pval <- function(p.value) {
        unclass(symnum(p.value, corr = FALSE, na = FALSE, cutpoints = c(
            0,
            0.001, 0.01, 0.05, 1
        ), symbols = c(
            "***", "**", "*",
            "NS"
        )))
    }
    if (p_stars == TRUE & p_val == TRUE) {
        p <- p + geom_text(data = pval, aes(
            x = variable, y = y,
            label = paste(stars.pval(pval))
        ), size = 3, hjust = 1)
    }
    if (p_stars == FALSE & p_val == TRUE & (bar_chart == FALSE |
        (bar_chart == TRUE & bar_chart_stacked == FALSE))) {
        # if predifined stats are provided use them as p without filtering
        if (!is.null(predefined_stats)) {
            p <- p + geom_text(data = pval, aes(
                x = variable, y = y,
                label = round(pval, 3)
            ), size = 9, hjust = 1) # size of LFC text
        } else {
            p <- p + geom_text(
                data = pval, aes(
                    x = variable, y = y,
                    label = ifelse(pval < 0.05, paste(format.pval(pval,
                        1, 0.001,
                        nsmall = 3
                    )), "")
                ), size = 5, hjust = 1,
                fontface = "bold"
            )
            p <- p + geom_text(data = pval, aes(
                x = variable, y = y,
                label = ifelse(pval >= 0.05, paste(format.pval(pval,
                    1, 0.001,
                    nsmall = 3
                )), "")
            ), size = 5, hjust = 1)
        }
        if (p_adjust) {
            p <- p + geom_text(data = pval, aes(
                x = variable,
                y = y_adjust, label = ifelse(p_adjust < 0.05,
                    paste(format.pval(p_adjust, 1, 0.001, nsmall = 3)),
                    ""
                )
            ), size = 5, hjust = 1, fontface = "bold")
            p <- p + geom_text(data = pval, aes(
                x = variable,
                y = y_adjust, label = ifelse(p_adjust >= 0.05,
                    paste(format.pval(p_adjust, 1, 0.001, nsmall = 3)),
                    ""
                )
            ), size = 5, hjust = 1)
        }
        if (stats == "mixed" & !is.null(facet_wrap)) {
            if (bar_chart == FALSE) {
                p <- p + expand_limits(y = 2)
            }
            if (bar_chart == TRUE) {
                p <- p + expand_limits(y = max(molten_mean2$value))
            }
        }
    }
    p
}
