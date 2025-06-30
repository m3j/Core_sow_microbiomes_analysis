library(ggpubr)
library(ggplot2)

## Set theme for alpha diversity plots
theme_alpha <- theme(
    axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1,
        size = 22,
        color = "grey30"
    ),
    strip.text.x = element_text(size = 22, color = "gray30", face = "bold"),
    strip.background = element_rect(color = "gray20"),
    panel.grid = element_line(color = "gray80"),
    panel.grid.major.x = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    axis.ticks.x = element_line(color = "grey30"),
    axis.ticks.y = element_line(color = "grey30"),
    axis.title.x = element_text(size = 22, color = "black"),
    axis.title.y = element_text(size = 22, color = "black"),
    axis.text.y = element_text(size = 22, color = "grey30"),
    legend.title = element_text(size = 22, color = "black"),
    legend.text = element_text(size = 17),
    legend.spacing.y = unit(0.05, "cm"),
    legend.key.size = unit(0.5, "cm"),
    plot.tag = element_text(size = 22, face = "bold"),
    plot.title = element_text(size = 34)
)


# function to plot a transparent boxplot with filled points and half violin plot showing the distribution on the right side of the boxplot. Comparing means with a standard t-test.
box_plotting <- function(metadata, x, y, fill, x_label, y_label, legend_pos, col_vector, label_vector, legend_title = fill) {
    plt <-
        metadata %>%
        ggplot(aes_string(x = x, y = y)) +
        ggdist::stat_halfeye(
            adjust = .5,
            width = .6,
            .width = 0,
            justification = -.3,
            point_colour = NA
        ) +
        geom_point(aes_string(fill = fill),
            shape = 21,
            size = 2,
            alpha = .8,
            position = position_jitter(
                seed = 1, width = .1
            )
        ) +
        geom_boxplot(
            alpha = 0.5,
            width = .25,
            outlier.shape = NA
        ) +
        ggpubr::stat_compare_means(vjust = -2, method = "t.test", size = 6) +
        labs(
            x = x_label,
            y = y_label
        ) +
        scale_color_brewer(palette = "Dark2") +
        theme_alpha +
        theme(legend.position = legend_pos)

    if (is.factor(metadata[[fill]]) || is.character(metadata[[fill]])) {
        # Discrete scale
        plt <- plt + scale_fill_discrete(
            type = col_vector,
            name = paste0(legend_title, ":"),
            labels = label_vector
        )
    } else {
        # Continuous scale
        plt <- plt + scale_fill_gradient(
            type = col_vector,
            name = paste0(legend_title, ":"),
            labels = label_vector
        )
    }

    return(plt)
}

# function to plot a transparent boxplot with filled points with pairwise comparison of means with more than two groups with a standard t-test. including half violin distribution on the right of each boxplot
box_plotting_compare <- function(metadata, x, y, fill, x_label, y_label, legend_pos, compare_list, color_vector, legend_title = fill, label_vector = NA) {
    plt <-
        metadata %>%
        ggplot(aes_string(x = x, y = y)) +
        ggdist::stat_halfeye(
            adjust = .5,
            width = .6,
            .width = 0,
            justification = -.3,
            point_colour = NA
        ) +
        geom_point(aes_string(fill = fill),
            shape = 21,
            size = 2,
            alpha = .8,
            position = position_jitter(
                seed = 1, width = .1
            )
        ) +
        geom_boxplot(
            alpha = 0.5,
            width = .25,
            outlier.shape = NA
        ) +
        ggpubr::stat_compare_means(
            comparisons = compare_list, label = "p.format",
            size = 6, method = "t.test", exact = FALSE, , hide.ns = TRUE,
        ) +
        labs(
            x = x_label,
            y = y_label
        ) +
        scale_color_brewer(palette = "Dark2") +
        scale_color_manual(values = eval(parse(text = color_vector))) +
        theme_alpha +
        theme(legend.position = legend_pos)

    if (is.factor(metadata[[fill]]) || is.character(metadata[[fill]])) {
        # Discrete scale
        if (all(is.na(label_vector))) {
            plt <- plt + scale_fill_discrete(
                type = eval(parse(text = color_vector)),
                name = paste0(legend_title, ":")
            )
        } else {
            plt <- plt + scale_fill_discrete(
                type = eval(parse(text = color_vector)),
                name = paste0(legend_title, ":"),
                labels = label_vector
            )
        }
    } else {
        # Continuous scale
        plt <- plt + scale_fill_gradient(name = paste0(legend_title, ":"))
    }

    return(plt)
}

# function to plot a transparent boxplot with unfilled points and half violin plot showing the distribution on the right side of the boxplot. Comparing means with a standard t-test.
box_plotting_no_fill <- function(metadata, x, y, x_label, y_label) {
    plt <-
        metadata %>%
        ggplot(aes_string(x = x, y = y)) +
        ggdist::stat_halfeye(
            adjust = .5,
            width = .6,
            .width = 0,
            justification = -.3,
            point_colour = NA
        ) +
        geom_point(
            shape = 21,
            size = 2,
            alpha = .8,
            position = position_jitter(
                seed = 1, width = .1
            )
        ) +
        geom_boxplot(
            alpha = 0.5,
            width = .25,
            outlier.shape = NA
        ) +
        ggpubr::stat_compare_means(vjust = -1.2, method = "t.test", size = 6) +
        labs(
            x = x_label,
            y = y_label
        ) +
        theme_alpha

    return(plt)
}
